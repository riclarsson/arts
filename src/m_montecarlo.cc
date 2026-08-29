#include <arts_constants.h>
#include <atm_path.h>
#include <mc_antenna.h>
#include <path_point.h>
#include <physics_funcs.h>
#include <scattering_species.h>
#include <workspace.h>

#include <algorithm>
#include <chrono>
#include <cmath>

namespace {

struct ScalarOptics {
  Numeric extinction{};
  Numeric absorption{};
  Numeric temperature{};
};

ScalarOptics scalar_optics(const Workspace&                ws,
                           Numeric                         frequency,
                           const PropagationPathPoint&     point,
                           const AtmPoint&                 atm,
                           const ArrayOfScatteringSpecies& scattering_species,
                           const Agenda&                   spectral_propmat_agenda) {
  const AscendingGrid freq_grid{frequency};
  PropmatVector gas;
  StokvecVector src;
  PropmatMatrix gas_jac;
  StokvecMatrix src_jac;
  spectral_propmat_agendaExecute(ws,
                                 gas,
                                 src,
                                 gas_jac,
                                 src_jac,
                                 freq_grid,
                                 Vector3{0.0, 0.0, 0.0},
                                 JacobianTargets{},
                                 {},
                                 point,
                                 atm,
                                 spectral_propmat_agenda);

  Numeric particle_ext = 0.0;
  Numeric particle_abs = 0.0;
  if (not scattering_species.species.empty()) {
    auto za = std::make_shared<scattering::ZenithAngleGrid>(
        scattering::IrregularZenithAngleGrid(Vector{0.0, 180.0}));
    const auto bulk = scattering_species.get_bulk_scattering_properties_tro_gridded(atm, freq_grid, za);
    particle_ext    = bulk.extinction_matrix[0, 0, 0];
    particle_abs    = bulk.absorption_vector[0, 0, 0];
  }

  return {.extinction = gas[0].A() + particle_ext,
          .absorption = gas[0].A() + particle_abs,
          .temperature = atm.temperature};
}

PropagationPathPoint interpolate(const PropagationPathPoint& a,
                                 const PropagationPathPoint& b,
                                 Numeric                     x) {
  PropagationPathPoint p = a;
  p.pos_type              = PathPositionType::atm;
  p.los_type              = PathPositionType::atm;
  p.pos                   = (1.0 - x) * a.pos + x * b.pos;
  p.los                   = (1.0 - x) * a.los + x * b.los;
  p.nreal                 = (1.0 - x) * a.nreal + x * b.nreal;
  p.ngroup                = (1.0 - x) * a.ngroup + x * b.ngroup;
  return p;
}

Stokvec boundary_radiance(const Workspace&            ws,
                          Numeric                     frequency,
                          const PropagationPathPoint& point,
                          const SurfaceField&         surf_field,
                          const SubsurfaceField&      subsurf_field,
                          const Agenda&               spectral_rad_space_agenda,
                          const Agenda&               spectral_rad_surface_agenda) {
  StokvecVector value;
  StokvecMatrix jac;
  const AscendingGrid grid{frequency};
  if (point.los_type == PathPositionType::space) {
    spectral_rad_space_agendaExecute(ws, value, jac, grid, JacobianTargets{}, point, spectral_rad_space_agenda);
  } else if (point.los_type == PathPositionType::surface) {
    spectral_rad_surface_agendaExecute(
        ws, value, jac, grid, JacobianTargets{}, point, surf_field, subsurf_field, spectral_rad_surface_agenda);
  } else {
    ARTS_USER_ERROR("Monte Carlo path ended without reaching space or the surface")
  }
  return value[0];
}

Vector2 isotropic_los(RandomNumberGenerator<>& rng) {
  auto uniform = rng.get<std::uniform_real_distribution>(0.0, 1.0);
  return {Conversion::rad2deg(std::acos(1.0 - 2.0 * uniform())), 360.0 * uniform() - 180.0};
}

struct Collision {
  bool                     happened{false};
  PropagationPathPoint     point{};
  AtmPoint                 atm{};
  ScalarOptics             optics{};
  PropagationPathPoint     end{};
};

Collision next_collision(const Workspace&                ws,
                         RandomNumberGenerator<>&         rng,
                         Numeric                         frequency,
                         const Vector3&                  pos,
                         const Vector2&                  los,
                         const AtmField&                 atm_field,
                         const SurfaceField&             surf_field,
                         const ArrayOfScatteringSpecies& scattering_species,
                         const Agenda&                   ray_path_observer_agenda,
                         const Agenda&                   spectral_propmat_agenda) {
  ArrayOfPropagationPathPoint path;
  ray_path_observer_agendaExecute(ws, path, pos, los, ray_path_observer_agenda);
  ARTS_USER_ERROR_IF(path.empty(), "The ray-path agenda returned an empty path")
  if (path.size() == 1) return {.end = path.back()};

  const auto atm       = forward_atm_path(path, atm_field);
  const auto distances = path::distance(path, surf_field.ellipsoid);
  auto uniform         = rng.get<std::uniform_real_distribution>(0.0, 1.0);
  const Numeric target = -std::log(std::max(uniform(), std::numeric_limits<Numeric>::min()));
  Numeric tau          = 0.0;
  auto previous = scalar_optics(ws, frequency, path[0], atm[0], scattering_species, spectral_propmat_agenda);

  for (Size i = 1; i < path.size(); ++i) {
    const auto current = scalar_optics(ws, frequency, path[i], atm[i], scattering_species, spectral_propmat_agenda);
    const Numeric dtau = 0.5 * (previous.extinction + current.extinction) * distances[i - 1];
    if (dtau > 0.0 and tau + dtau >= target) {
      const Numeric x = std::clamp((target - tau) / dtau, 0.0, 1.0);
      auto point      = interpolate(path[i - 1], path[i], x);
      auto atm_point  = atm_field.at(point.pos);
      auto optics     = scalar_optics(ws, frequency, point, atm_point, scattering_species, spectral_propmat_agenda);
      return {.happened = true, .point = point, .atm = std::move(atm_point), .optics = optics, .end = path.back()};
    }
    tau += std::max(0.0, dtau);
    previous = current;
  }
  return {.end = path.back()};
}

Numeric phase_weight(const ArrayOfScatteringSpecies& species,
                     const AtmPoint&                 atm,
                     Numeric                         frequency,
                     const Vector2&                  old_los,
                     const Vector2&                  new_los,
                     Numeric                         scattering_coefficient) {
  const auto unit = [](const Vector2& los) {
    const Numeric za = Conversion::deg2rad(los[0]);
    const Numeric aa = Conversion::deg2rad(los[1]);
    return Vector3{std::sin(za) * std::sin(aa), std::sin(za) * std::cos(aa), std::cos(za)};
  };
  const Vector3 a = unit(old_los), b = unit(new_los);
  const Numeric angle = Conversion::rad2deg(std::acos(std::clamp(dot(a, b), -1.0, 1.0)));
  auto za = std::make_shared<scattering::ZenithAngleGrid>(
      scattering::IrregularZenithAngleGrid(Vector{angle}));
  const auto bulk = species.get_bulk_scattering_properties_tro_gridded(atm, AscendingGrid{frequency}, za);
  ARTS_USER_ERROR_IF(not bulk.phase_matrix, "MCGeneral requires phase matrices from all scattering species")
  return 4.0 * Constant::pi * (*bulk.phase_matrix)[0, 0, 0, 0] / scattering_coefficient;
}

}  // namespace

void MCGeneral(const Workspace&                ws,
               Stokvec&                        mc_spectral_rad,
               Stokvec&                        mc_error,
               Index&                          mc_iteration_count,
               const AtmField&                 atm_field,
               const SurfaceField&             surf_field,
               const SubsurfaceField&          subsurf_field,
               const ArrayOfScatteringSpecies& scattering_species,
               const MCAntenna&                mc_antenna,
               const Agenda&                   ray_path_observer_agenda,
               const Agenda&                   spectral_propmat_agenda,
               const Agenda&                   spectral_rad_space_agenda,
               const Agenda&                   spectral_rad_surface_agenda,
               const Numeric&                  frequency,
               const Vector3&                  sensor_pos,
               const Vector2&                  sensor_los,
               const Index&                    mc_seed,
               const Index&                    mc_min_iter,
               const Index&                    mc_max_iter,
               const Index&                    mc_max_scatorder,
               const Numeric&                  mc_std_err,
               const Numeric&                  mc_max_time) try {
  ARTS_TIME_REPORT
  ARTS_USER_ERROR_IF(frequency <= 0.0, "frequency must be positive")
  ARTS_USER_ERROR_IF(mc_min_iter < 1 or mc_max_iter < mc_min_iter, "Require 1 <= mc_min_iter <= mc_max_iter")
  ARTS_USER_ERROR_IF(mc_max_scatorder < 0, "mc_max_scatorder cannot be negative")
  ARTS_USER_ERROR_IF(mc_std_err < 0.0 or mc_max_time < 0.0, "Stopping thresholds cannot be negative")

  RandomNumberGenerator<> rng(mc_seed);
  rng.force_seed(mc_seed);
  const Matrix33 ant_to_enu = rotmat_enu(sensor_los);
  Stokvec sum{0.0}, sumsq{0.0};
  const auto started = std::chrono::steady_clock::now();
  mc_iteration_count = 0;

  for (; mc_iteration_count < mc_max_iter; ++mc_iteration_count) {
    const auto [initial_los, ray_rotation] = mc_antenna.draw_los(rng, ant_to_enu, sensor_los);
    Vector3 pos = sensor_pos;
    Vector2 los = initial_los;
    Numeric weight = 1.0;
    Stokvec sample{0.0};

    for (Index order = 0;;) {
      const auto collision = next_collision(ws,
                                            rng,
                                            frequency,
                                            pos,
                                            los,
                                            atm_field,
                                            surf_field,
                                            scattering_species,
                                            ray_path_observer_agenda,
                                            spectral_propmat_agenda);
      if (not collision.happened) {
        sample = weight * boundary_radiance(ws,
                                            frequency,
                                            collision.end,
                                            surf_field,
                                            subsurf_field,
                                            spectral_rad_space_agenda,
                                            spectral_rad_surface_agenda);
        break;
      }

      auto uniform = rng.get<std::uniform_real_distribution>(0.0, 1.0);
      const Numeric extinction = collision.optics.extinction;
      const Numeric scattering = std::max(0.0, extinction - collision.optics.absorption);
      if (scattering <= 0.0 or uniform() >= scattering / extinction) {
        sample[0] = weight * planck(frequency, collision.optics.temperature);
        break;
      }
      if (++order > mc_max_scatorder) break;

      const Vector2 new_los = isotropic_los(rng);
      weight *= phase_weight(scattering_species, collision.atm, frequency, los, new_los, scattering);
      pos = collision.point.pos;
      los = new_los;
    }

    sum += sample;
    for (Index s = 0; s < 4; ++s) sumsq[s] += sample[s] * sample[s];
    const Index n = mc_iteration_count + 1;
    if (n < mc_min_iter) continue;
    const Numeric nn = static_cast<Numeric>(n);

    Numeric largest_error = 0.0;
    for (Index s = 0; s < 4; ++s) {
      const Numeric mean = sum[s] / nn;
      largest_error = std::max(largest_error, std::sqrt(std::max(0.0, sumsq[s] / nn - mean * mean) / nn));
    }
    const Numeric elapsed = std::chrono::duration<Numeric>(std::chrono::steady_clock::now() - started).count();
    if ((mc_std_err > 0.0 and largest_error <= mc_std_err) or (mc_max_time > 0.0 and elapsed >= mc_max_time)) {
      ++mc_iteration_count;
      break;
    }
  }

  ARTS_USER_ERROR_IF(mc_iteration_count == 0, "MCGeneral produced no samples")
  const Numeric n = static_cast<Numeric>(mc_iteration_count);
  for (Index s = 0; s < 4; ++s) {
    mc_spectral_rad[s] = sum[s] / n;
    const Numeric variance = std::max(0.0, sumsq[s] / n - mc_spectral_rad[s] * mc_spectral_rad[s]);
    mc_error[s] = std::sqrt(variance / n);
  }
}
ARTS_METHOD_ERROR_CATCH

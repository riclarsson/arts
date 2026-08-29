#include <arts_constants.h>
#include <atm_path.h>
#include <jacobian.h>
#include <path_point.h>
#include <rtepack.h>
#include <scattering_species.h>
#include <sensor_meta_info.h>
#include <time_report.h>
#include <workspace.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <string_view>
#include <vector>

namespace {

Muelmat as_muelmat(const Matrix& x) {
  ARTS_USER_ERROR_IF(x.nrows() != 4 or x.ncols() != 4, "Expected a 4 by 4 phase matrix, got {:B,}", x.shape())
  Muelmat out{0.0};
  for (Index i = 0; i < 4; ++i)
    for (Index j = 0; j < 4; ++j) out[i, j] = x[i, j];
  return out;
}

Numeric liebe93_k2(const Numeric frequency, const Numeric temperature) {
  ARTS_USER_ERROR_IF(
      temperature < Constant::temperature_at_0c - 40.0 or temperature > Constant::temperature_at_0c + 100.0,
      "Automatic radar k2 requires {} K <= ze_tref <= {} K; got {} K",
      Constant::temperature_at_0c - 40.0,
      Constant::temperature_at_0c + 100.0,
      temperature)
  ARTS_USER_ERROR_IF(frequency < 10e9 or frequency > 1000e9,
                     "Automatic radar k2 uses Liebe-93 and requires 10 GHz <= frequency <= 1000 GHz; got {} Hz",
                     frequency)

  const Numeric theta = 1.0 - 300.0 / temperature;
  const Numeric e0    = 77.66 - 103.3 * theta;
  const Numeric e1    = 0.0671 * e0;
  const Numeric f1    = 20.2 + 146.0 * theta + 316.0 * theta * theta;
  const Numeric e2    = 3.52;
  const Numeric f2    = 39.8 * f1;
  const Complex ifGHz{0.0, frequency / 1e9};
  const Complex n  = std::sqrt(e2 + (e1 - e2) / (1.0 - ifGHz / f2) + (e0 - e1) / (1.0 - ifGHz / f1));
  const Complex n2 = n * n;
  return std::norm((n2 - 1.0) / (n2 + 2.0));
}

Numeric ze_factor(const Numeric frequency, const Numeric ze_tref, const Numeric k2) {
  const Numeric dielectric = k2 < 0.0 ? liebe93_k2(frequency, ze_tref) : k2;
  ARTS_USER_ERROR_IF(dielectric <= 0.0, "Radar k2 must be positive, or negative to request Liebe-93; got {}", k2)
  const Numeric wavelength = Constant::speed_of_light / frequency;
  return 4e18 * std::pow(wavelength, 4) / (std::pow(Constant::pi, 4) * dielectric);
}

struct OpticalProfile {
  std::vector<PropmatVector> extinction;
  std::vector<PropmatVector> gas_extinction;
  std::vector<PropmatVector> particle_extinction;
  std::vector<MuelmatVector> backscatter;
  std::vector<PropmatMatrix> extinction_jac;
  std::vector<MuelmatMatrix> backscatter_jac;
};

OpticalProfile optical_profile(const Workspace&                   ws,
                               const AscendingGrid&               frequency,
                               const ArrayOfPropagationPathPoint& ray_path,
                               const AtmField&                    atm_field,
                               const ArrayOfScatteringSpecies&    scattering_species,
                               const Agenda&                      spectral_propmat_agenda,
                               const JacobianTargets&             jac_targets,
                               const Numeric                      pext_scaling) {
  const Size np = ray_path.size();
  const Size nf = frequency.size();
  const Size nq = jac_targets.target_count();

  OpticalProfile out;
  out.extinction.resize(np, PropmatVector(nf));
  out.gas_extinction.resize(np, PropmatVector(nf));
  out.particle_extinction.resize(np, PropmatVector(nf));
  out.backscatter.resize(np, MuelmatVector(nf, Muelmat{0.0}));
  out.extinction_jac.resize(np, PropmatMatrix(nq, nf));
  out.backscatter_jac.resize(np, MuelmatMatrix(nq, nf, Muelmat{0.0}));

  const ArrayOfAtmPoint atm_path = forward_atm_path(ray_path, atm_field);
  const auto            za_grid =
      std::make_shared<scattering::ZenithAngleGrid>(scattering::IrregularZenithAngleGrid(Vector{180.0}));

  for (Size ip = 0; ip < np; ++ip) {
    PropmatVector gas;
    StokvecVector source;
    PropmatMatrix gas_jac;
    StokvecMatrix source_jac;
    spectral_propmat_agendaExecute(ws,
                                   gas,
                                   source,
                                   gas_jac,
                                   source_jac,
                                   frequency,
                                   Vector3{0.0, 0.0, 0.0},
                                   jac_targets,
                                   {},
                                   ray_path[ip],
                                   atm_path[ip],
                                   spectral_propmat_agenda);
    ARTS_USER_ERROR_IF(gas.size() != nf, "spectral_propmat_agenda returned {} frequencies, expected {}", gas.size(), nf)

    const auto bulk = scattering_species.get_bulk_scattering_properties_tro_gridded(atm_path[ip], frequency, za_grid);
    ARTS_USER_ERROR_IF(not bulk.phase_matrix.has_value(),
                       "Radar single scattering requires phase matrices from all scattering species")

    const auto& phase = *bulk.phase_matrix;
    for (Size iv = 0; iv < nf; ++iv) {
      out.gas_extinction[ip][iv]       = gas[iv];
      out.particle_extinction[ip][iv]  = Propmat{bulk.extinction_matrix[0, iv, 0]};
      out.extinction[ip][iv]           = gas[iv];
      out.extinction[ip][iv]          += pext_scaling * out.particle_extinction[ip][iv];
      out.backscatter[ip][iv]          = as_muelmat(scattering::expand_phase_matrix(phase[0, iv, 0, joker]));
    }
    for (const auto& target : jac_targets.atm) {
      const auto dbulk = scattering_species.get_bulk_scattering_properties_tro_gridded_derivative(
          atm_path[ip], frequency, za_grid, target.type);
      ARTS_USER_ERROR_IF(not dbulk.phase_matrix, "Radar derivatives require scattering phase matrices")
      for (Size iv = 0; iv < nf; ++iv) {
        out.extinction_jac[ip][target.target_pos, iv] = gas_jac[target.target_pos, iv];
        out.extinction_jac[ip][target.target_pos, iv] +=
            pext_scaling * Propmat{dbulk.extinction_matrix[0, iv, 0]};
        out.backscatter_jac[ip][target.target_pos, iv] =
            as_muelmat(scattering::expand_phase_matrix((*dbulk.phase_matrix)[0, iv, 0, joker]));
      }
    }
  }

  return out;
}

enum class RadarRangeMode { Altitude, Distance, RoundTripTime };

RadarRangeMode range_mode(const String& mode, const Matrix& limits) {
  if (mode == "Altitude") return RadarRangeMode::Altitude;
  if (mode == "Distance") return RadarRangeMode::Distance;
  if (mode == "RoundTripTime") return RadarRangeMode::RoundTripTime;
  if (mode == "Legacy") {
    Numeric largest = -std::numeric_limits<Numeric>::infinity();
    for (Index i = 0; i < limits.nrows(); ++i)
      for (Index j = 0; j < limits.ncols(); ++j) largest = std::max(largest, limits[i, j]);
    return largest > 1.0 ? RadarRangeMode::Altitude : RadarRangeMode::RoundTripTime;
  }
  ARTS_USER_ERROR("Unknown radar range_mode \"{}\". Expected Altitude, Distance, RoundTripTime, or Legacy.", mode)
}

Vector range_coordinate(const ArrayOfPropagationPathPoint& ray_path,
                        const Vector3&                     observer,
                        const Vector2&                     ellipsoid,
                        const RadarRangeMode               mode) {
  const Size np = ray_path.size();
  Vector     out(np);
  if (mode == RadarRangeMode::Altitude) {
    for (Size ip = 0; ip < np; ++ip) out[ip] = ray_path[ip].pos[0];
    return out;
  }

  const Numeric initial = path::distance(observer, ray_path.front().pos, ellipsoid);
  if (mode == RadarRangeMode::Distance) {
    out[0]               = initial;
    const Vector segment = path::distance(ray_path, ellipsoid);
    for (Size ip = 1; ip < np; ++ip) out[ip] = out[ip - 1] + segment[ip - 1];
    return out;
  }

  out[0]               = 2.0 * initial / Constant::speed_of_light;
  const Vector segment = path::distance(ray_path, ellipsoid);
  for (Size ip = 1; ip < np; ++ip) {
    out[ip] =
        out[ip - 1] + segment[ip - 1] * (ray_path[ip - 1].ngroup + ray_path[ip].ngroup) / Constant::speed_of_light;
  }
  return out;
}

template <class T> std::optional<T> bin_average(const Vector& x, const std::vector<T>& y, Numeric lo, Numeric hi) {
  ARTS_USER_ERROR_IF(x.size() != y.size(), "Internal radar range/profile size mismatch")
  if (x.size() < 2) return std::nullopt;

  const bool    increasing = x.back() > x.front();
  const Numeric xmin       = std::min(x.front(), x.back());
  const Numeric xmax       = std::max(x.front(), x.back());
  lo                       = std::max(lo, xmin);
  hi                       = std::min(hi, xmax);
  if (not(lo < hi)) return std::nullopt;

  auto interpolate = [&](const Numeric q) {
    Size i = 0;
    if (increasing) {
      while (i + 1 < x.size() and x[i + 1] < q) ++i;
    } else {
      while (i + 1 < x.size() and x[i + 1] > q) ++i;
    }
    if (i + 1 >= x.size()) return y.back();
    const Numeric w = (q - x[i]) / (x[i + 1] - x[i]);
    return (1.0 - w) * y[i] + w * y[i + 1];
  };

  std::vector<Numeric> q{lo};
  q.reserve(x.size() + 2);
  for (auto v : x)
    if (v > lo and v < hi) q.push_back(v);
  q.push_back(hi);
  std::ranges::sort(q);

  T integral{0.0};
  T past = interpolate(q.front());
  for (Size i = 1; i < q.size(); ++i) {
    const T now  = interpolate(q[i]);
    integral    += 0.5 * (q[i] - q[i - 1]) * (past + now);
    past         = now;
  }
  return integral * (1.0 / (hi - lo));
}

struct PreparedSimulation {
  const AscendingGrid*        frequency{};
  const SensorPosLosVector*   poslos{};
  Size                        iposlos{};
  ArrayOfPropagationPathPoint ray_path{};
};

struct RadarResult {
  Vector measurement;
  Matrix jacobian;
  Matrix auxiliary;
};

RadarResult radar_forward(const Workspace&                       ws,
                          const ArrayOfSensorObsel&              sensor,
                          const Matrix&                          limits,
                          const std::vector<PreparedSimulation>& simulations,
                          const AtmField&                        atm_field,
                          const SurfaceField&                    surf_field,
                          const ArrayOfScatteringSpecies&        scattering_species,
                          const Agenda&                          spectral_propmat_agenda,
                          const Stokvec&                         transmitted,
                          const RadarRangeMode                   coordinate_mode,
                          const String&                          unit,
                          const Numeric                          ze_tref,
                          const Numeric                          k2,
                          const Numeric                          dbze_min,
                          const Numeric                          pext_scaling,
                          const ArrayOfString&                   aux_vars,
                          const JacobianTargets&                 jac_targets) {
  const Size        ny = sensor.size();
  const Size        nq = jac_targets.target_count();
  RadarResult       out{Vector(ny, 0.0), Matrix(ny, jac_targets.x_size(), 0.0), Matrix(aux_vars.size(), ny, 0.0)};
  std::vector<bool> seen(ny, false);
  Vector            scalar_weight_sum(ny, 0.0);

  for (const auto& simulation : simulations) {
    const auto& frequency = *simulation.frequency;
    const auto& poslos    = *simulation.poslos;
    const auto& ray_path  = simulation.ray_path;
    const Size  iposlos   = simulation.iposlos;
    if (ray_path.size() < 2) continue;

    const OpticalProfile optical = optical_profile(ws,
                                                   frequency,
                                                   ray_path,
                                                   atm_field,
                                                   scattering_species,
                                                   spectral_propmat_agenda,
                                                   jac_targets,
                                                   pext_scaling);
    const Size np = ray_path.size();
    const Size nf = frequency.size();

    std::vector<StokvecVector> attenuated(np, StokvecVector(nf));
    std::vector<StokvecVector> unattenuated(np, StokvecVector(nf));
    std::vector<StokvecVector> gas_ext(np, StokvecVector(nf));
    std::vector<StokvecVector> particle_ext(np, StokvecVector(nf));
    std::vector<Muelmat>       cumulative(nf, Muelmat::id());
    std::vector<MuelmatMatrix> cumulative_jac(nf, MuelmatMatrix(np, nq, Muelmat{0.0}));
    std::vector<std::vector<StokvecVector>> attenuated_jac(
        np * nq, std::vector<StokvecVector>(np, StokvecVector(nf, Stokvec{0.0})));

    for (Size ip = 0; ip < np; ++ip) {
      if (ip > 0) {
        const Numeric ds = path::distance(ray_path[ip - 1].pos, ray_path[ip].pos, surf_field.ellipsoid);
        for (Size iv = 0; iv < nf; ++iv) {
          const Muelmat previous = cumulative[iv];
          const rtepack::tran segment(optical.extinction[ip - 1][iv], optical.extinction[ip][iv], ds);
          const Muelmat transmission = segment();
          for (Size jp = 0; jp < np; ++jp)
            for (Size iq = 0; iq < nq; ++iq)
              cumulative_jac[iv][jp, iq] = transmission * cumulative_jac[iv][jp, iq];
          for (Size iq = 0; iq < nq; ++iq) {
            cumulative_jac[iv][ip - 1, iq] +=
                segment.deriv(transmission,
                              optical.extinction[ip - 1][iv],
                              optical.extinction[ip][iv],
                              optical.extinction_jac[ip - 1][iq, iv],
                              ds,
                              0.0) * previous;
            cumulative_jac[iv][ip, iq] +=
                segment.deriv(transmission,
                              optical.extinction[ip - 1][iv],
                              optical.extinction[ip][iv],
                              optical.extinction_jac[ip][iq, iv],
                              ds,
                              0.0) * previous;
          }
          cumulative[iv] = transmission * previous;
        }
      }
      for (Size iv = 0; iv < nf; ++iv) {
        unattenuated[ip][iv] = optical.backscatter[ip][iv] * transmitted;
        attenuated[ip][iv]   = cumulative[iv] * (optical.backscatter[ip][iv] * (cumulative[iv] * transmitted));
        gas_ext[ip][iv]      = Stokvec{optical.gas_extinction[ip][iv].A(), 0.0, 0.0, 0.0};
        particle_ext[ip][iv] = Stokvec{pext_scaling * optical.particle_extinction[ip][iv].A(), 0.0, 0.0, 0.0};
        const Stokvec incoming  = cumulative[iv] * transmitted;
        const Stokvec scattered = optical.backscatter[ip][iv] * incoming;
        for (Size jp = 0; jp < np; ++jp) {
          for (Size iq = 0; iq < nq; ++iq) {
            const Muelmat& dtrans = cumulative_jac[iv][jp, iq];
            Stokvec derivative = dtrans * scattered +
                                 cumulative[iv] * (optical.backscatter[ip][iv] * (dtrans * transmitted));
            if (jp == ip)
              derivative += cumulative[iv] * (optical.backscatter_jac[ip][iq, iv] * incoming);
            attenuated_jac[jp * nq + iq][ip][iv] = derivative;
          }
        }
      }
    }

    const Vector  coordinate = range_coordinate(ray_path, poslos[iposlos].pos, surf_field.ellipsoid, coordinate_mode);
    const Numeric background = ray_path.back().los_type == PathPositionType::space
                                   ? 0.0
                                   : (ray_path.back().los_type == PathPositionType::surface ? 1.0 : 2.0);

    for (Size iy = 0; iy < ny; ++iy) {
      const auto& obsel = sensor[iy];
      if (not obsel.same_freqs(frequency)) continue;

      StokvecVector bin_signal(nf, Stokvec{0.0});
      StokvecVector bin_backscatter(nf, Stokvec{0.0});
      StokvecVector bin_gas(nf, Stokvec{0.0});
      StokvecVector bin_particle(nf, Stokvec{0.0});
      std::vector<StokvecVector> bin_jac(np * nq, StokvecVector(nf, Stokvec{0.0}));
      bool          inside = false;
      for (Size iv = 0; iv < nf; ++iv) {
        std::vector<Stokvec> a(np), b(np), g(np), p(np);
        for (Size ip = 0; ip < np; ++ip) {
          a[ip] = attenuated[ip][iv];
          b[ip] = unattenuated[ip][iv];
          g[ip] = gas_ext[ip][iv];
          p[ip] = particle_ext[ip][iv];
        }
        if (auto avg = bin_average(coordinate, a, limits[iy, 0], limits[iy, 1])) {
          inside              = true;
          bin_signal[iv]      = *avg;
          bin_backscatter[iv] = *bin_average(coordinate, b, limits[iy, 0], limits[iy, 1]);
          bin_gas[iv]         = *bin_average(coordinate, g, limits[iy, 0], limits[iy, 1]);
          bin_particle[iv]    = *bin_average(coordinate, p, limits[iy, 0], limits[iy, 1]);

          const Numeric factor  = unit == "1" ? 1.0 : ze_factor(frequency[iv], ze_tref, k2);
          bin_signal[iv]       *= factor;
          bin_backscatter[iv]  *= factor;
          for (Size jp = 0; jp < np; ++jp) {
            for (Size iq = 0; iq < nq; ++iq) {
              std::vector<Stokvec> derivative(np);
              for (Size kp = 0; kp < np; ++kp)
                derivative[kp] = attenuated_jac[jp * nq + iq][kp][iv];
              bin_jac[jp * nq + iq][iv] =
                  factor * *bin_average(coordinate, derivative, limits[iy, 0], limits[iy, 1]);
            }
          }
        }
      }
      if (not inside) {
        if (not seen[iy]) {
          out.measurement[iy] = std::numeric_limits<Numeric>::quiet_NaN();
          for (Size ia = 0; ia < aux_vars.size(); ++ia)
            out.auxiliary[ia, iy] = std::numeric_limits<Numeric>::quiet_NaN();
        }
        continue;
      }

      if (not seen[iy] or std::isnan(out.measurement[iy])) {
        out.measurement[iy] = 0.0;
        for (Size ia = 0; ia < aux_vars.size(); ++ia) out.auxiliary[ia, iy] = 0.0;
      }
      const Numeric signal  = obsel.sumup(bin_signal, iposlos);
      out.measurement[iy]  += signal;
      for (const auto& target : jac_targets.atm) {
        const auto& field = atm_field[target.type];
        for (Size jp = 0; jp < np; ++jp) {
          const Numeric local = obsel.sumup(bin_jac[jp * nq + target.target_pos], iposlos);
          for (const auto& [index, weight] : field.flat_weight(ray_path[jp].pos))
            out.jacobian[iy, target.x_start + index] += weight * local;
        }
      }

      const Numeric scalar_weight  = obsel.sumup(Stokvec{1.0, 0.0, 0.0, 0.0}, iposlos);
      scalar_weight_sum[iy]       += scalar_weight;
      for (Size ia = 0; ia < aux_vars.size(); ++ia) {
        if (aux_vars[ia] == "Radiative background") {
          out.auxiliary[ia, iy] += scalar_weight * background;
        } else if (aux_vars[ia] == "Backscattering") {
          out.auxiliary[ia, iy] += obsel.sumup(bin_backscatter, iposlos);
        } else if (aux_vars[ia] == "Abs species extinction") {
          out.auxiliary[ia, iy] += obsel.sumup(bin_gas, iposlos);
        } else if (aux_vars[ia] == "Particle extinction") {
          out.auxiliary[ia, iy] += obsel.sumup(bin_particle, iposlos);
        }
      }
      seen[iy] = true;
    }
  }

  for (Size ia = 0; ia < aux_vars.size(); ++ia) {
    if (aux_vars[ia] == "Backscattering") continue;
    for (Size iy = 0; iy < ny; ++iy)
      if (seen[iy] and scalar_weight_sum[iy] != 0.0) out.auxiliary[ia, iy] /= scalar_weight_sum[iy];
  }

  if (unit == "dBZe") {
    const Numeric ze_min = std::pow(10.0, dbze_min / 10.0);
    for (Size iy = 0; iy < ny; ++iy) {
      if (std::isnan(out.measurement[iy])) continue;
      if (out.measurement[iy] <= ze_min) {
        out.measurement[iy] = dbze_min;
        out.jacobian[iy, joker] = 0.0;
      } else {
        const Numeric factor = 10.0 / (std::log(10.0) * out.measurement[iy]);
        out.jacobian[iy, joker] *= factor;
        out.measurement[iy] = 10.0 * std::log10(out.measurement[iy]);
      }
    }
    for (Size ia = 0; ia < aux_vars.size(); ++ia) {
      if (aux_vars[ia] != "Backscattering") continue;
      for (Size iy = 0; iy < ny; ++iy) {
        auto& value = out.auxiliary[ia, iy];
        if (not std::isnan(value)) value = value <= ze_min ? dbze_min : 10.0 * std::log10(value);
      }
    }
  }

  return out;
}

}  // namespace

void measurement_sensorAddSimpleRadar(ArrayOfSensorObsel&    measurement_sensor,
                                      ArrayOfSensorMetaInfo& measurement_sensor_meta,
                                      Matrix&                radar_range_limits,
                                      const AscendingGrid&   freq_grid,
                                      const Vector3&         pos,
                                      const Vector2&         los,
                                      const Stokvec&         pol,
                                      const AscendingGrid&   range_bins) try {
  ARTS_TIME_REPORT
  ARTS_USER_ERROR_IF(freq_grid.empty(), "Radar frequency grid cannot be empty")
  ARTS_USER_ERROR_IF(range_bins.size() < 2, "Radar range_bins needs at least two edges")

  const Size old_n = measurement_sensor.size();
  ARTS_USER_ERROR_IF(radar_range_limits.nrows() != 0 and
                         (radar_range_limits.ncols() != 2 or static_cast<Size>(radar_range_limits.nrows()) != old_n),
                     "radar_range_limits must be empty or have shape [{}, 2], got {:B,}",
                     old_n,
                     radar_range_limits.shape())

  const Size nbins = range_bins.size() - 1;
  const Size nnew  = freq_grid.size() * nbins;
  Matrix     new_limits(old_n + nnew, 2);
  if (old_n) new_limits[Range(0, old_n), joker] = radar_range_limits;

  auto poslos = std::make_shared<const SensorPosLosVector>(SensorPosLosVector{SensorPosLos{pos, los}});
  measurement_sensor.reserve(old_n + nnew);
  measurement_sensor_meta.reserve(measurement_sensor_meta.size() + freq_grid.size());

  Size row = old_n;
  for (Size iv = 0; iv < freq_grid.size(); ++iv) {
    auto                        one_frequency = std::make_shared<const AscendingGrid>(AscendingGrid{freq_grid[iv]});
    sensor::SparseStokvecMatrix weights(1, 1);
    weights[0, 0] = pol;

    for (Size ib = 0; ib < nbins; ++ib) {
      measurement_sensor.emplace_back(one_frequency, poslos, weights);
      new_limits[row, 0] = range_bins[ib];
      new_limits[row, 1] = range_bins[ib + 1];
      ++row;
    }

    SortedGriddedField1 meta;
    meta.data_name     = std::format("radar-{}-Hz", freq_grid[iv]);
    meta.gridname<0>() = "range";
    Vector centers(nbins);
    for (Size ib = 0; ib < nbins; ++ib) centers[ib] = std::midpoint(range_bins[ib], range_bins[ib + 1]);
    meta.grid<0>() = AscendingGrid{std::move(centers)};
    meta.data.resize(nbins);
    meta.data = 0.0;
    measurement_sensor_meta.push_back(SensorMetaInfo{.data = std::move(meta)});
  }
  radar_range_limits = std::move(new_limits);
}
ARTS_METHOD_ERROR_CATCH

void measurement_vecFromRadarSingleScattering(const Workspace&                ws,
                                              Vector&                         measurement_vec,
                                              Matrix&                         measurement_jac,
                                              Matrix&                         radar_aux,
                                              const ArrayOfSensorObsel&       measurement_sensor,
                                              const Matrix&                   radar_range_limits,
                                              const JacobianTargets&          jac_targets,
                                              const AtmField&                 atm_field,
                                              const SurfaceField&             surf_field,
                                              const ArrayOfScatteringSpecies& scattering_species,
                                              const Agenda&                   ray_path_observer_agenda,
                                              const Agenda&                   spectral_propmat_agenda,
                                              const Stokvec&                  transmitted_stokes,
                                              const String&                   range_mode_string,
                                              const String&                   unit,
                                              const Numeric&                  ze_tref,
                                              const Numeric&                  k2,
                                              const Numeric&                  dbze_min,
                                              const Numeric&                  pext_scaling,
                                              const ArrayOfString&            aux_vars) try {
  ARTS_TIME_REPORT
  const Size ny = measurement_sensor.size();
  ARTS_USER_ERROR_IF(ny == 0, "measurement_sensor cannot be empty for a radar calculation")
  ARTS_USER_ERROR_IF(radar_range_limits.nrows() != static_cast<Index>(ny) or radar_range_limits.ncols() != 2,
                     "radar_range_limits must have shape [{}, 2], got {:B,}",
                     ny,
                     radar_range_limits.shape())
  for (Size iy = 0; iy < ny; ++iy) {
    measurement_sensor[iy].check();
    ARTS_USER_ERROR_IF(not(radar_range_limits[iy, 0] < radar_range_limits[iy, 1]),
                       "Radar range limits in row {} are not strictly increasing: [{}, {}]",
                       iy,
                       radar_range_limits[iy, 0],
                       radar_range_limits[iy, 1])
  }
  ARTS_USER_ERROR_IF(transmitted_stokes.I() != 1.0,
                     "The transmitted Stokes I component must equal one, got {}",
                     transmitted_stokes.I())
  ARTS_USER_ERROR_IF(unit != "1" and unit != "Ze" and unit != "dBZe",
                     "Radar unit must be \"1\", \"Ze\", or \"dBZe\"; got \"{}\"",
                     unit)
  ARTS_USER_ERROR_IF(
      pext_scaling < 0.0 or pext_scaling > 2.0, "pext_scaling must be between 0 and 2; got {}", pext_scaling)
  ARTS_USER_ERROR_IF(scattering_species.species.empty(), "Radar single scattering requires scattering species")
  ARTS_USER_ERROR_IF(not jac_targets.surf.empty() or not jac_targets.subsurf.empty() or not jac_targets.line.empty() or
                         not jac_targets.sensor.empty() or not jac_targets.error.empty(),
                     "The radar forward model currently supports atmospheric Jacobian targets only")
  for (const auto& name : aux_vars) {
    ARTS_USER_ERROR_IF(name != "Radiative background" and name != "Backscattering" and
                           name != "Abs species extinction" and name != "Particle extinction",
                       "Unknown radar auxiliary variable \"{}\"",
                       name)
  }

  const RadarRangeMode            coordinate_mode = range_mode(range_mode_string, radar_range_limits);
  const SensorSimulations         collected       = collect_simulations(measurement_sensor);
  std::vector<PreparedSimulation> simulations;
  simulations.reserve(collected.size());
  for (const auto& simulation : collected) {
    const auto& [frequency, poslos, iposlos] = simulation;
    PreparedSimulation prepared{.frequency = &frequency, .poslos = &poslos, .iposlos = iposlos};
    ray_path_observer_agendaExecute(
        ws, prepared.ray_path, poslos[iposlos].pos, poslos[iposlos].los, ray_path_observer_agenda);
    simulations.push_back(std::move(prepared));
  }

  const auto calculate = [&](const AtmField& field) {
    return radar_forward(ws,
                         measurement_sensor,
                         radar_range_limits,
                         simulations,
                         field,
                         surf_field,
                         scattering_species,
                         spectral_propmat_agenda,
                         transmitted_stokes,
                         coordinate_mode,
                         unit,
                         ze_tref,
                         k2,
                         dbze_min,
                         pext_scaling,
                         aux_vars,
                         jac_targets);
  };

  RadarResult base = calculate(atm_field);
  measurement_vec  = std::move(base.measurement);
  measurement_jac  = std::move(base.jacobian);
  radar_aux        = std::move(base.auxiliary);

  const Size nx = jac_targets.x_size();
  if (nx == 0) return;

  Vector state(nx, 0.0);
  for (const auto& target : jac_targets.atm) target.update_state(state, atm_field);
  for (const auto& target : jac_targets.atm) target.update_jac(measurement_jac, state, atm_field);
}
ARTS_METHOD_ERROR_CATCH

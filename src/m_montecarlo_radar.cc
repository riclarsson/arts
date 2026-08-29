#include <arts_constants.h>
#include <atm_path.h>
#include <mc_antenna.h>
#include <path_point.h>
#include <rtepack.h>
#include <scattering_species.h>
#include <workspace.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

Muelmat as_muelmat(const Matrix& x) {
  ARTS_USER_ERROR_IF(x.nrows() != 4 or x.ncols() != 4, "Expected a 4 by 4 phase matrix, got {:B,}", x.shape())
  Muelmat out{0.0};
  for (Index i = 0; i < 4; ++i)
    for (Index j = 0; j < 4; ++j) out[i, j] = x[i, j];
  return out;
}

Index range_bin(const AscendingGrid& edges, const Numeric range) {
  const auto it = std::upper_bound(edges.begin(), edges.end(), range);
  if (it == edges.begin() or it == edges.end()) return -1;
  return static_cast<Index>(it - edges.begin() - 1);
}

Numeric ze_factor(const Numeric frequency, const Numeric k2) {
  const Numeric wavelength = Constant::speed_of_light / frequency;
  return 4e18 * std::pow(wavelength, 4) / (std::pow(Constant::pi, 4) * k2);
}

struct OpticalPoint {
  Propmat extinction{};
  Muelmat backscatter{0.0};
};

OpticalPoint optical_point(const Workspace&                ws,
                           const Numeric                   frequency,
                           const PropagationPathPoint&     ray_point,
                           const AtmPoint&                 atm_point,
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
                                 ray_point,
                                 atm_point,
                                 spectral_propmat_agenda);

  auto za_grid    = std::make_shared<scattering::ZenithAngleGrid>(scattering::IrregularZenithAngleGrid(Vector{180.0}));
  const auto bulk = scattering_species.get_bulk_scattering_properties_tro_gridded(atm_point, freq_grid, za_grid);

  ARTS_USER_ERROR_IF(not bulk.phase_matrix.has_value(), "MCRadar requires phase matrices from all scattering species")

  const auto& phase   = *bulk.phase_matrix;
  const auto  compact = phase[0, 0, 0, joker];

  OpticalPoint out;
  out.extinction   = gas[0];
  out.extinction  += Propmat{bulk.extinction_matrix[0, 0, 0]};
  out.backscatter  = as_muelmat(scattering::expand_phase_matrix(compact));
  return out;
}

}  // namespace

void MCRadar(const Workspace&                ws,
             Matrix&                         radar_signal,
             Matrix&                         radar_error,
             const AtmField&                 atm_field,
             const SurfaceField&             surf_field,
             const ArrayOfScatteringSpecies& scattering_species,
             const MCAntenna&                mc_antenna,
             const Agenda&                   ray_path_observer_agenda,
             const Agenda&                   spectral_propmat_agenda,
             const Numeric&                  frequency,
             const Vector3&                  sensor_pos,
             const Vector2&                  sensor_los,
             const Stokvec&                  mc_y_tx,
             const AscendingGrid&            range_bins,
             const Index&                    mc_seed,
             const Index&                    mc_max_iter,
             const Index&                    mc_max_scatorder,
             const Numeric&                  k2,
             const String&                   unit) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(frequency <= 0.0, "frequency must be positive")
  ARTS_USER_ERROR_IF(range_bins.size() < 2, "range_bins must contain at least two edges")
  ARTS_USER_ERROR_IF(range_bins.front() < 0.0, "range_bins cannot contain negative distances")
  ARTS_USER_ERROR_IF(mc_max_iter <= 0, "mc_max_iter must be positive")
  ARTS_USER_ERROR_IF(mc_max_scatorder != 1,
                     "This first ARTS3 MCRadar implementation supports mc_max_scatorder=1; got {}",
                     mc_max_scatorder)
  ARTS_USER_ERROR_IF(scattering_species.species.empty(), "MCRadar requires at least one scattering species")
  ARTS_USER_ERROR_IF(unit != "1" and unit != "Ze", "unit must be either \"1\" or \"Ze\"")
  ARTS_USER_ERROR_IF(k2 <= 0.0, "k2 must be positive")

  const Index nbins = static_cast<Index>(range_bins.size()) - 1;
  radar_signal.resize(nbins, 4);
  radar_error.resize(nbins, 4);
  radar_signal = 0.0;
  radar_error  = 0.0;
  Matrix squared_sum(nbins, 4, 0.0);

  RandomNumberGenerator<> rng(mc_seed);
  // The general RNG avoids reusing seeds process-wide.  A forward model with
  // an explicit seed instead promises reproducibility across repeated calls.
  rng.force_seed(mc_seed);
  const Matrix33 ant_to_enu = rotmat_enu(sensor_los);
  Matrix33       enu_to_ant;
  for (Index i = 0; i < 3; ++i)
    for (Index j = 0; j < 3; ++j) enu_to_ant[i, j] = ant_to_enu[j, i];

  for (Index iter = 0; iter < mc_max_iter; ++iter) {
    Matrix photon(nbins, 4, 0.0);
    const auto [sampled_los, ray_rotation] = mc_antenna.draw_los(rng, ant_to_enu, sensor_los);
    const Numeric antenna_weight           = mc_antenna.return_los(ray_rotation, enu_to_ant);
    const Muelmat tx_rotation              = rotmat_stokes(1.0, 1.0, ant_to_enu, ray_rotation);
    const Stokvec transmitted              = tx_rotation * mc_y_tx;

    ArrayOfPropagationPathPoint ray_path;
    ray_path_observer_agendaExecute(ws, ray_path, sensor_pos, sampled_los, ray_path_observer_agenda);
    if (ray_path.size() < 2) continue;

    const ArrayOfAtmPoint atm_path  = forward_atm_path(ray_path, atm_field);
    const Vector          distances = path::distance(ray_path, surf_field.ellipsoid);
    Muelmat               outbound  = Muelmat::id();
    Numeric               travelled = 0.0;
    OpticalPoint          previous =
        optical_point(ws, frequency, ray_path.front(), atm_path.front(), scattering_species, spectral_propmat_agenda);

    for (Size ip = 1; ip < ray_path.size(); ++ip) {
      const Numeric ds = distances[ip - 1];
      if (ds <= 0.0) continue;

      const OpticalPoint current =
          optical_point(ws, frequency, ray_path[ip], atm_path[ip], scattering_species, spectral_propmat_agenda);
      const Muelmat segment                = rtepack::tran(previous.extinction, current.extinction, ds)();
      const Muelmat midpoint_transmission  = rtepack::tran(previous.extinction, current.extinction, 0.5 * ds)();
      const Muelmat to_midpoint            = midpoint_transmission * outbound;
      travelled                           += 0.5 * ds;

      const Index ibin = range_bin(range_bins, travelled);
      if (ibin >= 0) {
        const Muelmat backscatter  = 0.5 * (previous.backscatter + current.backscatter);
        const Stokvec contribution = antenna_weight * ds * (to_midpoint * (backscatter * (to_midpoint * transmitted)));
        for (Index s = 0; s < 4; ++s) photon[ibin, s] += contribution[s];
      }

      travelled += 0.5 * ds;
      if (travelled >= range_bins.back()) break;
      outbound = segment * outbound;
      previous = current;
    }

    for (Index b = 0; b < nbins; ++b) {
      const Numeric width = range_bins[b + 1] - range_bins[b];
      for (Index s = 0; s < 4; ++s) {
        const Numeric value  = photon[b, s] / width;
        radar_signal[b, s]  += value;
        squared_sum[b, s]   += value * value;
      }
    }
  }

  const Numeric factor = unit == "Ze" ? ze_factor(frequency, k2) / (2.0 * Constant::pi) : 1.0;
  for (Index b = 0; b < nbins; ++b) {
    for (Index s = 0; s < 4; ++s) {
      const Numeric mean     = radar_signal[b, s] / static_cast<Numeric>(mc_max_iter);
      const Numeric variance = std::max(0.0, squared_sum[b, s] / static_cast<Numeric>(mc_max_iter) - mean * mean);
      radar_signal[b, s]     = factor * mean;
      radar_error[b, s]      = factor * std::sqrt(variance / static_cast<Numeric>(mc_max_iter));
    }
  }
}
ARTS_METHOD_ERROR_CATCH

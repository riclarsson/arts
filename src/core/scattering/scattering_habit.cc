#include "scattering_habit.h"

#include "configtypes.h"
#include "interpolation.h"

namespace scattering {

ScatteringHabit::ScatteringHabit(const ParticleHabit& particle_habit_,
                                 const PSD&           psd_,
                                 Numeric              mass_size_rel_a_,
                                 Numeric              mass_size_rel_b_)
    : particle_habit(particle_habit_), mass_size_rel_a(mass_size_rel_a_), mass_size_rel_b(mass_size_rel_b_), psd(psd_) {
  if ((mass_size_rel_a < 0.0) || (mass_size_rel_b < 0.0)) {
    auto size_param = std::visit([](const auto& psd) { return psd.get_size_parameter(); }, psd);
    auto [sizes, mass_size_rel_a_, mass_size_rel_b_] = particle_habit.get_size_mass_info(size_param);
    mass_size_rel_a                                  = mass_size_rel_a_;
    mass_size_rel_b                                  = mass_size_rel_b_;
  }
}

namespace detail {

Numeric max_relative_difference(const Vector& v1, const Vector& v2) {
  if (v1.size() != v2.size()) { return 1.0; }
  Numeric max_diff = 0.0;

  for (size_t i = 0; i < v1.size(); ++i) {
    Numeric denom = std::max(std::abs(v1[i]), std::abs(v2[i]));
    if (denom > std::numeric_limits<Numeric>::epsilon()) {  // Avoid division by zero
      Numeric rel_diff = std::abs(v1[i] - v2[i]) / denom;
      max_diff         = std::max(max_diff, rel_diff);
    }
  }
  return max_diff;
}

}  // namespace detail

////////////////////////////////////////////////////////////////////////////////
// Bulk properties TRO gridded
////////////////////////////////////////////////////////////////////////////////

/** Calculate bulk scattering properties for an amtospheric point
 *
 * This funciont evaluates the PSD at the given point, interpolates the scattering data along the temperature
 * dimension, and sum up the particle scattering data to calculate the bulk scattering properties.
 *
 * @param point: The AtmPoint for which to calculate the bulk properties.
 * @param f_grid: The frequencies for which to calculate the bulk scattering properties.
 * @param f_tol: Maximum relative difference between frequency vectors to trigger re-interpolation of
 *     along the frequency grid.
 * @return A struct containing the calculated bulk-scattering properties.
 */
BulkScatteringPropertiesTROGridded ScatteringHabit::get_bulk_scattering_properties_tro_gridded(const AtmPoint& point,
                                                                                               const Vector&   f_grid,
                                                                                               const Numeric) const {
  auto  sizes       = particle_habit.get_sizes(SizeParameter::DVeq);
  Index n_particles = sizes.size();
  auto  pnd         = std::visit(
      [&point, &sizes, this](const auto& psd) { return psd.evaluate(point, sizes, mass_size_rel_a, mass_size_rel_b); },
      psd);

  if (!particle_habit.grids.has_value()) {
    ARTS_USER_ERROR("Particle habit must be brought on a shared grid before buld properties can be computed.")
  }

  auto    grids  = particle_habit.grids.value();
  GridPos interp = find_interp_weights(*grids.t_grid, point[AtmKey::t]);

  Index n_freqs = grids.f_grid->size();
  Index n_angs  = grid_size(*grids.za_scat_grid);

  Tensor4 phase_matrix(n_freqs, n_angs, 4, 4);
  Tensor3 extinction_matrix(n_freqs, 4, 4);
  Matrix  absorption_vector(n_freqs, 4);

  for (Index part_ind = 0; part_ind < n_particles; ++part_ind) {
    try {
      auto ssd =
          std::get<SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>>(particle_habit[part_ind]);
      for (Index f_ind = 0; f_ind < n_freqs; ++f_ind) {
        for (Index ang_ind = 0; ang_ind < n_angs; ++ang_ind) {
          phase_matrix[f_ind, ang_ind, 0, 0] +=
              pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, ang_ind, 0];
          phase_matrix[f_ind, ang_ind, 0, 0] +=
              pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, ang_ind, 0];
          phase_matrix[f_ind, ang_ind, 0, 1] +=
              pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, ang_ind, 1];
          phase_matrix[f_ind, ang_ind, 0, 1] +=
              pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, ang_ind, 1];
          phase_matrix[f_ind, ang_ind, 1, 0] +=
              pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, ang_ind, 1];
          phase_matrix[f_ind, ang_ind, 1, 0] +=
              pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, ang_ind, 1];
          phase_matrix[f_ind, ang_ind, 1, 1] +=
              pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, ang_ind, 2];
          phase_matrix[f_ind, ang_ind, 1, 1] +=
              pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, ang_ind, 2];
          phase_matrix[f_ind, ang_ind, 2, 2] +=
              pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, ang_ind, 3];
          phase_matrix[f_ind, ang_ind, 2, 2] +=
              pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, ang_ind, 3];
          phase_matrix[f_ind, ang_ind, 2, 3] +=
              pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, ang_ind, 4];
          phase_matrix[f_ind, ang_ind, 2, 3] +=
              pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, ang_ind, 4];
          phase_matrix[f_ind, ang_ind, 3, 2] +=
              pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, ang_ind, 4];
          phase_matrix[f_ind, ang_ind, 3, 2] +=
              pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, ang_ind, 4];
          phase_matrix[f_ind, ang_ind, 3, 3] +=
              pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, ang_ind, 5];
          phase_matrix[f_ind, ang_ind, 3, 3] +=
              pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, ang_ind, 5];
        }
        for (Index stokes_ind = 0; stokes_ind < 4; ++stokes_ind) {
          extinction_matrix[f_ind, stokes_ind, stokes_ind] +=
              pnd[part_ind] * interp.fd[1] * ssd.extinction_matrix[interp.idx, f_ind, stokes_ind];
          extinction_matrix[f_ind, stokes_ind, stokes_ind] +=
              pnd[part_ind] * interp.fd[0] * ssd.extinction_matrix[interp.idx + 1, f_ind, stokes_ind];
          absorption_vector[f_ind, stokes_ind] +=
              pnd[part_ind] * interp.fd[1] * ssd.absorption_vector[interp.idx, f_ind, stokes_ind];
          absorption_vector[f_ind, stokes_ind] +=
              pnd[part_ind] * interp.fd[0] * ssd.absorption_vector[interp.idx + 1, f_ind, stokes_ind];
        }
      }
    } catch (const std::bad_variant_access& e) {
      ARTS_USER_ERROR("Scattering habit must be in TRO gridded format to extract bulk scattering properties.");
    }
  }

  auto f_diff = detail::max_relative_difference(f_grid, *grids.f_grid);

  // If frequency grids are the same, return calculated bulk properties.
  if (f_diff < 1e-3) { return BulkScatteringPropertiesTROGridded(phase_matrix, extinction_matrix, absorption_vector); }

  Index n_freqs_new = f_grid.size();

  // Otherwise perform frequency interpolation.
  ArrayOfGridPos interp_weights;
  gridpos(interp_weights, *grids.f_grid, f_grid);
  Tensor4 phase_matrix_new(n_freqs_new, n_angs, 4, 4);
  Tensor3 extinction_matrix_new(n_freqs_new, 4, 4);
  Matrix  absorption_vector_new(n_freqs_new, 4);

  for (Size f_ind = 0; f_ind < f_grid.size(); ++f_ind) {
    auto weights = interp_weights[f_ind];
    for (Index za_scat_ind = 0; za_scat_ind < n_angs; ++za_scat_ind) {
      phase_matrix_new[f_ind, za_scat_ind, 0, 0] += weights.fd[1] * phase_matrix[weights.idx, za_scat_ind, 0, 0];
      phase_matrix_new[f_ind, za_scat_ind, 0, 0] += weights.fd[0] * phase_matrix[weights.idx + 1, za_scat_ind, 0, 0];
      phase_matrix_new[f_ind, za_scat_ind, 0, 1] += weights.fd[1] * phase_matrix[weights.idx, za_scat_ind, 0, 1];
      phase_matrix_new[f_ind, za_scat_ind, 0, 1] += weights.fd[0] * phase_matrix[weights.idx + 1, za_scat_ind, 0, 1];
      phase_matrix_new[f_ind, za_scat_ind, 1, 0] += weights.fd[1] * phase_matrix[weights.idx, za_scat_ind, 1, 0];
      phase_matrix_new[f_ind, za_scat_ind, 1, 0] += weights.fd[0] * phase_matrix[weights.idx + 1, za_scat_ind, 1, 0];
      phase_matrix_new[f_ind, za_scat_ind, 1, 1] += weights.fd[1] * phase_matrix[weights.idx, za_scat_ind, 1, 1];
      phase_matrix_new[f_ind, za_scat_ind, 1, 1] += weights.fd[0] * phase_matrix[weights.idx + 1, za_scat_ind, 1, 1];
      phase_matrix_new[f_ind, za_scat_ind, 2, 2] += weights.fd[1] * phase_matrix[weights.idx, za_scat_ind, 2, 2];
      phase_matrix_new[f_ind, za_scat_ind, 2, 2] += weights.fd[0] * phase_matrix[weights.idx + 1, za_scat_ind, 2, 2];
      phase_matrix_new[f_ind, za_scat_ind, 2, 3] += weights.fd[1] * phase_matrix[weights.idx, za_scat_ind, 2, 3];
      phase_matrix_new[f_ind, za_scat_ind, 2, 3] += weights.fd[0] * phase_matrix[weights.idx + 1, za_scat_ind, 2, 3];
      phase_matrix_new[f_ind, za_scat_ind, 3, 2] += weights.fd[1] * phase_matrix[weights.idx, za_scat_ind, 3, 2];
      phase_matrix_new[f_ind, za_scat_ind, 3, 2] += weights.fd[0] * phase_matrix[weights.idx + 1, za_scat_ind, 3, 2];
      phase_matrix_new[f_ind, za_scat_ind, 3, 3] += weights.fd[1] * phase_matrix[weights.idx, za_scat_ind, 3, 3];
      phase_matrix_new[f_ind, za_scat_ind, 3, 3] += weights.fd[0] * phase_matrix[weights.idx + 1, za_scat_ind, 3, 3];
    }

    for (Index stokes_ind = 0; stokes_ind < 4; ++stokes_ind) {
      extinction_matrix_new[f_ind, stokes_ind, stokes_ind] +=
          weights.fd[1] * extinction_matrix[weights.idx, stokes_ind, stokes_ind];
      extinction_matrix_new[f_ind, stokes_ind, stokes_ind] +=
          weights.fd[0] * extinction_matrix[weights.idx + 1, stokes_ind, stokes_ind];
      absorption_vector_new[f_ind, stokes_ind] += weights.fd[1] * absorption_vector[weights.idx, stokes_ind];
      absorption_vector_new[f_ind, stokes_ind] += weights.fd[0] * absorption_vector[weights.idx + 1, stokes_ind];
    }
  }
  return BulkScatteringPropertiesTROGridded(phase_matrix_new, extinction_matrix_new, absorption_vector_new);
}

BulkScatteringProperties<Format::TRO, Representation::Gridded>
ScatteringHabit::get_bulk_scattering_properties_tro_gridded(const AtmPoint&                  point,
                                                            const Vector&                    f_grid,
                                                            std::shared_ptr<ZenithAngleGrid> za_scat_grid) const {
  const auto sizes = particle_habit.get_sizes(std::visit([](const auto& p) { return p.get_size_parameter(); }, psd));
  const auto pnd   = std::visit(
      [&point, &sizes, this](const auto& p) { return p.evaluate(point, sizes, mass_size_rel_a, mass_size_rel_b); },
      psd);
  ARTS_USER_ERROR_IF(pnd.size() != static_cast<Size>(particle_habit.size()), "PSD and particle-habit sizes differ.")

  auto grids = ScatteringDataGrids(
      std::make_shared<Vector>(Vector{point.temperature}), std::make_shared<Vector>(f_grid), std::move(za_scat_grid));
  using Bulk = BulkScatteringProperties<Format::TRO, Representation::Gridded>;
  std::optional<Bulk> result;
  for (Index i = 0; i < particle_habit.size(); ++i) {
    auto data = std::visit([&grids](const auto& ssd) { return ssd_to_tro_gridded(grids, ssd); }, particle_habit[i]);
    if (pnd[i] == 0.0) {
      if (data.phase_matrix) std::fill_n(data.phase_matrix->data_handle(), data.phase_matrix->size(), Numeric{0.0});
      std::fill_n(data.extinction_matrix.data_handle(), data.extinction_matrix.size(), Numeric{0.0});
      std::fill_n(data.absorption_vector.data_handle(), data.absorption_vector.size(), Numeric{0.0});
    } else {
      if (data.phase_matrix) *data.phase_matrix *= pnd[i];
      data.extinction_matrix *= pnd[i];
      data.absorption_vector *= pnd[i];
    }
    Bulk bulk{std::move(data.phase_matrix), std::move(data.extinction_matrix), std::move(data.absorption_vector)};
    if (result) {
      *result += bulk;
    } else {
      result = std::move(bulk);
    }
  }
  ARTS_USER_ERROR_IF(not result, "Cannot calculate bulk properties for an empty particle habit.")
  return std::move(*result);
}

BulkScatteringProperties<Format::TRO, Representation::Gridded>
ScatteringHabit::get_bulk_scattering_properties_tro_gridded_derivative(const AtmPoint&                  point,
                                                                       const Vector&                    f_grid,
                                                                       std::shared_ptr<ZenithAngleGrid> za_scat_grid,
                                                                       const AtmKeyVal&                 target) const {
  const auto sizes = particle_habit.get_sizes(std::visit([](const auto& p) { return p.get_size_parameter(); }, psd));
  const auto pnd   = std::visit(
      [&point, &sizes, this](const auto& p) {
        return p.evaluate_with_derivatives(point, sizes, mass_size_rel_a, mass_size_rel_b);
      },
      psd);
  ARTS_USER_ERROR_IF(pnd.values.size() != static_cast<Size>(particle_habit.size()),
                     "PSD and particle-habit sizes differ.")

  const auto*   property   = std::get_if<ScatteringSpeciesProperty>(&target);
  const Vector* derivative = nullptr;
  if (property) {
    const auto iter = pnd.derivatives.find(*property);
    if (iter != pnd.derivatives.end()) derivative = &iter->second;
  }

  auto grids = ScatteringDataGrids(
      std::make_shared<Vector>(Vector{point.temperature}), std::make_shared<Vector>(f_grid), za_scat_grid);
  using Bulk = BulkScatteringProperties<Format::TRO, Representation::Gridded>;
  std::optional<Bulk> result;
  for (Index i = 0; i < particle_habit.size(); ++i) {
    auto data = std::visit([&grids](const auto& ssd) { return ssd_to_tro_gridded(grids, ssd); }, particle_habit[i]);
    const Numeric weight = derivative ? (*derivative)[i] : 0.0;
    if (weight == 0.0) {
      if (data.phase_matrix) std::fill_n(data.phase_matrix->data_handle(), data.phase_matrix->size(), Numeric{0.0});
      std::fill_n(data.extinction_matrix.data_handle(), data.extinction_matrix.size(), Numeric{0.0});
      std::fill_n(data.absorption_vector.data_handle(), data.absorption_vector.size(), Numeric{0.0});
    } else {
      if (data.phase_matrix) *data.phase_matrix *= weight;
      data.extinction_matrix *= weight;
      data.absorption_vector *= weight;
    }
    Bulk bulk{std::move(data.phase_matrix), std::move(data.extinction_matrix), std::move(data.absorption_vector)};

    if (target == AtmKeyVal{AtmKey::t} and pnd.values[i] != 0.0) {
      Bulk temperature_derivative = std::visit(
          [&](const auto& ssd) -> Bulk {
            if constexpr (std::remove_cvref_t<decltype(ssd)>::get_format() == Format::ARO) {
              ARTS_USER_ERROR("TRO radar derivatives cannot use ARO particle data")
            } else {
              const auto gridded = [&]() {
                if constexpr (std::remove_cvref_t<decltype(ssd)>::get_representation() == Representation::Gridded)
                  return ssd;
                else
                  return ssd.to_gridded();
              }();
              ARTS_USER_ERROR_IF(not gridded.phase_matrix, "Temperature derivatives require particle phase matrices")
              const Vector& temperatures = *gridded.phase_matrix->get_t_grid();
              ARTS_USER_ERROR_IF(temperatures.size() < 2,
                                 "Temperature-dependent particle derivatives require at least two temperatures")

              auto value_at = [&](Numeric temperature) {
                ScatteringDataGrids local_grids(
                    std::make_shared<Vector>(Vector{temperature}), std::make_shared<Vector>(f_grid), za_scat_grid);
                auto value = ssd_to_tro_gridded(local_grids, ssd);
                return Bulk{std::move(value.phase_matrix),
                            std::move(value.extinction_matrix),
                            std::move(value.absorption_vector)};
              };
              auto       upper   = std::ranges::lower_bound(temperatures, point.temperature);
              const Size located = static_cast<Size>(upper - temperatures.begin());
              Size       right   = std::clamp<Size>(located, 1, temperatures.size() - 1);
              Size       left    = right - 1;

              if (upper != temperatures.end() and *upper == point.temperature and located == right and right > 0 and
                  right + 1 < temperatures.size()) {
                auto below   = value_at(temperatures[right - 1]);
                auto middle  = value_at(temperatures[right]);
                auto above   = value_at(temperatures[right + 1]);
                below       *= -0.5 / (temperatures[right] - temperatures[right - 1]);
                middle      *= 0.5 / (temperatures[right] - temperatures[right - 1]) -
                               0.5 / (temperatures[right + 1] - temperatures[right]);
                above       *= 0.5 / (temperatures[right + 1] - temperatures[right]);
                below       += middle;
                below       += above;
                return below;
              }

              auto          below            = value_at(temperatures[left]);
              auto          above            = value_at(temperatures[right]);
              const Numeric inverse_spacing  = 1.0 / (temperatures[right] - temperatures[left]);
              below                         *= -inverse_spacing;
              above                         *= inverse_spacing;
              below                         += above;
              return below;
            }
            std::unreachable();
          },
          particle_habit[i]);
      temperature_derivative *= pnd.values[i];
      bulk                   += temperature_derivative;
    }
    if (result)
      *result += bulk;
    else
      result = std::move(bulk);
  }
  ARTS_USER_ERROR_IF(not result, "Cannot calculate bulk properties for an empty particle habit.")
  return std::move(*result);
}

BulkScatteringProperties<Format::ARO, Representation::Gridded>
ScatteringHabit::get_bulk_scattering_properties_aro_gridded(const AtmPoint&                  point,
                                                            const Vector&                    f_grid,
                                                            const Vector&                    za_inc_grid,
                                                            const Vector&                    delta_aa_grid,
                                                            std::shared_ptr<ZenithAngleGrid> za_scat_grid) const {
  const auto sizes = particle_habit.get_sizes(std::visit([](const auto& p) { return p.get_size_parameter(); }, psd));
  const auto pnd   = std::visit(
      [&point, &sizes, this](const auto& p) { return p.evaluate(point, sizes, mass_size_rel_a, mass_size_rel_b); },
      psd);
  ARTS_USER_ERROR_IF(pnd.size() != static_cast<Size>(particle_habit.size()), "PSD and particle-habit sizes differ.")

  auto grids = ScatteringDataGrids(std::make_shared<Vector>(Vector{point.temperature}),
                                   std::make_shared<Vector>(f_grid),
                                   std::make_shared<Vector>(za_inc_grid),
                                   std::make_shared<Vector>(delta_aa_grid),
                                   std::move(za_scat_grid));
  using SSD  = SingleScatteringData<Numeric, Format::ARO, Representation::Gridded>;
  std::optional<BulkScatteringProperties<Format::ARO, Representation::Gridded>> result;
  for (Index i = 0; i < particle_habit.size(); ++i) {
    SSD data = std::visit([&grids](const auto& ssd) { return ssd_to_aro_gridded(grids, ssd); }, particle_habit[i]);
    if (pnd[i] == 0.0) {
      // A laboratory-frame conversion can contain an indeterminate 0*NaN at
      // angular-coordinate poles.  An absent particle population is exactly
      // zero regardless of its single-particle coordinate representation.
      if (data.phase_matrix) std::fill_n(data.phase_matrix->data_handle(), data.phase_matrix->size(), Numeric{0.0});
      std::fill_n(data.extinction_matrix.data_handle(), data.extinction_matrix.size(), Numeric{0.0});
      std::fill_n(data.absorption_vector.data_handle(), data.absorption_vector.size(), Numeric{0.0});
    } else {
      if (data.phase_matrix) *data.phase_matrix *= pnd[i];
      data.extinction_matrix *= pnd[i];
      data.absorption_vector *= pnd[i];
    }
    BulkScatteringProperties<Format::ARO, Representation::Gridded> bulk{
        std::move(data.phase_matrix), std::move(data.extinction_matrix), std::move(data.absorption_vector)};
    if (result) {
      *result += bulk;
    } else {
      result = std::move(bulk);
    }
  }
  ARTS_USER_ERROR_IF(not result, "Cannot calculate bulk properties for an empty particle habit.")
  return std::move(*result);
}

BulkScatteringProperties<Format::ARO, Representation::Gridded>
ScatteringHabit::get_bulk_scattering_properties_aro_gridded_derivative(const AtmPoint&                  point,
                                                                       const Vector&                    f_grid,
                                                                       const Vector&                    za_inc_grid,
                                                                       const Vector&                    delta_aa_grid,
                                                                       std::shared_ptr<ZenithAngleGrid> za_scat_grid,
                                                                       const AtmKeyVal&                 target) const {
  const auto sizes = particle_habit.get_sizes(std::visit([](const auto& p) { return p.get_size_parameter(); }, psd));
  const auto pnd   = std::visit(
      [&point, &sizes, this](const auto& p) {
        return p.evaluate_with_derivatives(point, sizes, mass_size_rel_a, mass_size_rel_b);
      },
      psd);
  ARTS_USER_ERROR_IF(pnd.values.size() != static_cast<Size>(particle_habit.size()),
                     "PSD and particle-habit sizes differ.")

  const auto*   property   = std::get_if<ScatteringSpeciesProperty>(&target);
  const Vector* derivative = nullptr;
  if (property) {
    const auto iter = pnd.derivatives.find(*property);
    if (iter != pnd.derivatives.end()) derivative = &iter->second;
  }

  auto grids = ScatteringDataGrids(std::make_shared<Vector>(Vector{point.temperature}),
                                   std::make_shared<Vector>(f_grid),
                                   std::make_shared<Vector>(za_inc_grid),
                                   std::make_shared<Vector>(delta_aa_grid),
                                   za_scat_grid);
  using Bulk = BulkScatteringProperties<Format::ARO, Representation::Gridded>;
  std::optional<Bulk> result;
  for (Index i = 0; i < particle_habit.size(); ++i) {
    auto data = std::visit([&grids](const auto& ssd) { return ssd_to_aro_gridded(grids, ssd); }, particle_habit[i]);
    const Numeric weight = derivative ? (*derivative)[i] : 0.0;
    if (weight == 0.0) {
      if (data.phase_matrix) std::fill_n(data.phase_matrix->data_handle(), data.phase_matrix->size(), Numeric{0.0});
      std::fill_n(data.extinction_matrix.data_handle(), data.extinction_matrix.size(), Numeric{0.0});
      std::fill_n(data.absorption_vector.data_handle(), data.absorption_vector.size(), Numeric{0.0});
    } else {
      if (data.phase_matrix) *data.phase_matrix *= weight;
      data.extinction_matrix *= weight;
      data.absorption_vector *= weight;
    }
    Bulk bulk{std::move(data.phase_matrix), std::move(data.extinction_matrix), std::move(data.absorption_vector)};

    if (target == AtmKeyVal{AtmKey::t} and pnd.values[i] != 0.0) {
      Bulk temperature_derivative = std::visit(
          [&](const auto& ssd) -> Bulk {
            const auto gridded = [&]() {
              if constexpr (std::remove_cvref_t<decltype(ssd)>::get_representation() == Representation::Gridded)
                return ssd;
              else
                return ssd.to_gridded();
            }();
            ARTS_USER_ERROR_IF(not gridded.phase_matrix, "Temperature derivatives require particle phase matrices")
            const Vector& temperatures = *gridded.phase_matrix->get_t_grid();
            ARTS_USER_ERROR_IF(temperatures.size() < 2,
                               "Temperature-dependent particle derivatives require at least two temperatures")

            auto value_at = [&](Numeric temperature) {
              ScatteringDataGrids local_grids(std::make_shared<Vector>(Vector{temperature}),
                                              std::make_shared<Vector>(f_grid),
                                              std::make_shared<Vector>(za_inc_grid),
                                              std::make_shared<Vector>(delta_aa_grid),
                                              za_scat_grid);
              auto                value = ssd_to_aro_gridded(local_grids, ssd);
              return Bulk{std::move(value.phase_matrix),
                          std::move(value.extinction_matrix),
                          std::move(value.absorption_vector)};
            };
            auto       upper   = std::ranges::lower_bound(temperatures, point.temperature);
            const Size located = static_cast<Size>(upper - temperatures.begin());
            const Size right   = std::clamp<Size>(located, 1, temperatures.size() - 1);
            const Size left    = right - 1;

            if (upper != temperatures.end() and *upper == point.temperature and located == right and right > 0 and
                right + 1 < temperatures.size()) {
              auto below   = value_at(temperatures[right - 1]);
              auto middle  = value_at(temperatures[right]);
              auto above   = value_at(temperatures[right + 1]);
              below       *= -0.5 / (temperatures[right] - temperatures[right - 1]);
              middle      *= 0.5 / (temperatures[right] - temperatures[right - 1]) -
                             0.5 / (temperatures[right + 1] - temperatures[right]);
              above       *= 0.5 / (temperatures[right + 1] - temperatures[right]);
              below       += middle;
              below       += above;
              return below;
            }

            auto          below            = value_at(temperatures[left]);
            auto          above            = value_at(temperatures[right]);
            const Numeric inverse_spacing  = 1.0 / (temperatures[right] - temperatures[left]);
            below                         *= -inverse_spacing;
            above                         *= inverse_spacing;
            below                         += above;
            return below;
          },
          particle_habit[i]);
      temperature_derivative *= pnd.values[i];
      bulk                   += temperature_derivative;
    }
    if (result)
      *result += bulk;
    else
      result = std::move(bulk);
  }
  ARTS_USER_ERROR_IF(not result, "Cannot calculate bulk properties for an empty particle habit.")
  return std::move(*result);
}

ScatteringTroSpectralVector ScatteringHabit::get_bulk_scattering_properties_tro_spectral(const AtmPoint& point,
                                                                                         const Vector&   f_grid,
                                                                                         const Index     degree
                                                                                         [[maybe_unused]]) const {
  auto   sizes = particle_habit.get_sizes(std::visit([](const auto& psd) { return psd.get_size_parameter(); }, psd));
  Index  n_particles = sizes.size();
  Vector bin_widths  = sizes;
  for (Index ind = 1; ind < n_particles - 1; ++ind) { bin_widths[ind] = 0.5 * (sizes[ind + 1] - sizes[ind - 1]); }
  bin_widths[0]               = sizes[1] - sizes[0];
  bin_widths[n_particles - 1] = sizes[n_particles - 1] - sizes[n_particles - 2];

  auto pnd = std::visit(
      [&point, &sizes, this](const auto& psd) { return psd.evaluate(point, sizes, mass_size_rel_a, mass_size_rel_b); },
      psd);
  for (Index ind = 0; ind < n_particles; ++ind) { pnd[ind] *= bin_widths[ind]; }

  if (!particle_habit.grids.has_value()) {
    ARTS_USER_ERROR("Particle habit must be brought on a shared grid before buld properties can be computed.")
  }

  auto    grids  = particle_habit.grids.value();
  GridPos interp = find_interp_weights(*grids.t_grid, point[AtmKey::t]);

  Index n_freqs = grids.f_grid->size();
  if (particle_habit.size() == 0) { ARTS_USER_ERROR("Encountered empty scattering habit without particles."); }
  auto  ssd      = std::get<SingleScatteringData<Numeric, Format::TRO, Representation::Spectral>>(particle_habit[0]);
  Index n_coeffs = ssd.phase_matrix.value().extent(2);

  SpecmatMatrix phase_matrix(n_freqs, n_coeffs, Specmat(0.0));
  PropmatVector extinction_matrix(n_freqs);
  StokvecVector absorption_vector(n_freqs);

  Numeric integral_fac = sqrt(4.0 * Constant::pi);

  for (Index part_ind = 0; part_ind < n_particles; ++part_ind) {
    try {
      auto ssd =
          std::get<SingleScatteringData<Numeric, Format::TRO, Representation::Spectral>>(particle_habit[part_ind]);
      for (Index f_ind = 0; f_ind < n_freqs; ++f_ind) {
        for (Index coeff_ind = 0; coeff_ind < n_coeffs; ++coeff_ind) {
          phase_matrix[f_ind, coeff_ind][0, 0] +=
              integral_fac * pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, coeff_ind, 0];
          phase_matrix[f_ind, coeff_ind][0, 0] += integral_fac * pnd[part_ind] * interp.fd[0] *
                                                  ssd.phase_matrix.value()[interp.idx + 1, f_ind, coeff_ind, 0];
          phase_matrix[f_ind, coeff_ind][0, 1] +=
              integral_fac * pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, coeff_ind, 1];
          phase_matrix[f_ind, coeff_ind][0, 1] += integral_fac * pnd[part_ind] * interp.fd[0] *
                                                  ssd.phase_matrix.value()[interp.idx + 1, f_ind, coeff_ind, 1];
          phase_matrix[f_ind, coeff_ind][1, 0] -=
              integral_fac * pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, coeff_ind, 1];
          phase_matrix[f_ind, coeff_ind][1, 0] -= integral_fac * pnd[part_ind] * interp.fd[0] *
                                                  ssd.phase_matrix.value()[interp.idx + 1, f_ind, coeff_ind, 1];
          phase_matrix[f_ind, coeff_ind][1, 1] +=
              integral_fac * pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, coeff_ind, 2];
          phase_matrix[f_ind, coeff_ind][1, 1] += integral_fac * pnd[part_ind] * interp.fd[0] *
                                                  ssd.phase_matrix.value()[interp.idx + 1, f_ind, coeff_ind, 2];
          phase_matrix[f_ind, coeff_ind][2, 2] +=
              integral_fac * pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, coeff_ind, 3];
          phase_matrix[f_ind, coeff_ind][2, 2] += integral_fac * pnd[part_ind] * interp.fd[0] *
                                                  ssd.phase_matrix.value()[interp.idx + 1, f_ind, coeff_ind, 3];
          phase_matrix[f_ind, coeff_ind][2, 3] +=
              integral_fac * pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, coeff_ind, 4];
          phase_matrix[f_ind, coeff_ind][2, 3] += integral_fac * pnd[part_ind] * interp.fd[0] *
                                                  ssd.phase_matrix.value()[interp.idx + 1, f_ind, coeff_ind, 4];
          phase_matrix[f_ind, coeff_ind][3, 2] -=
              integral_fac * pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, coeff_ind, 4];
          phase_matrix[f_ind, coeff_ind][3, 2] -= integral_fac * pnd[part_ind] * interp.fd[0] *
                                                  ssd.phase_matrix.value()[interp.idx + 1, f_ind, coeff_ind, 4];
          phase_matrix[f_ind, coeff_ind][3, 3] +=
              integral_fac * pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, coeff_ind, 5];
          phase_matrix[f_ind, coeff_ind][3, 3] += integral_fac * pnd[part_ind] * interp.fd[0] *
                                                  ssd.phase_matrix.value()[interp.idx + 1, f_ind, coeff_ind, 5];
        }

        extinction_matrix[f_ind].A() += pnd[part_ind] * interp.fd[1] * ssd.extinction_matrix[interp.idx, f_ind, 0];
        extinction_matrix[f_ind].A() += pnd[part_ind] * interp.fd[0] * ssd.extinction_matrix[interp.idx + 1, f_ind, 0];

        absorption_vector[f_ind][0] += pnd[part_ind] * interp.fd[1] * ssd.absorption_vector[interp.idx, f_ind, 0];
        absorption_vector[f_ind][0] += pnd[part_ind] * interp.fd[0] * ssd.absorption_vector[interp.idx + 1, f_ind, 0];
      }
    } catch (const std::bad_variant_access& e) {
      ARTS_USER_ERROR("Scattering habit must be in TRO gridded format to extract bulk scattering properties.");
    }
  }

  auto f_diff = detail::max_relative_difference(f_grid, *grids.f_grid);

  // If frequency grids are the same, return calculated bulk properties.
  if (f_diff < 1e-3) { return ScatteringTroSpectralVector(phase_matrix, extinction_matrix, absorption_vector); }

  Index n_freqs_new = f_grid.size();

  // Otherwise perform frequency interpolation.
  ArrayOfGridPos interp_weights{f_grid.size()};
  gridpos(interp_weights, *grids.f_grid, f_grid);
  SpecmatMatrix phase_matrix_new(n_freqs_new, n_coeffs, Specmat(0.0));
  PropmatVector extinction_matrix_new(n_freqs_new);
  StokvecVector absorption_vector_new(n_freqs_new);

  for (Size f_ind = 0; f_ind < f_grid.size(); ++f_ind) {
    auto weights = interp_weights[f_ind];
    for (Index coeff_ind = 0; coeff_ind < n_coeffs; ++coeff_ind) {
      phase_matrix_new[f_ind, coeff_ind] += weights.fd[1] * phase_matrix[weights.idx, coeff_ind];
      phase_matrix_new[f_ind, coeff_ind] += weights.fd[0] * phase_matrix[weights.idx + 1, coeff_ind];
    }
    extinction_matrix_new[f_ind] += weights.fd[1] * extinction_matrix[weights.idx];
    extinction_matrix_new[f_ind] += weights.fd[0] * extinction_matrix[weights.idx + 1];
    absorption_vector_new[f_ind] += weights.fd[1] * absorption_vector[weights.idx];
    absorption_vector_new[f_ind] += weights.fd[0] * absorption_vector[weights.idx + 1];
  }
  return ScatteringTroSpectralVector(phase_matrix_new, extinction_matrix_new, absorption_vector_new);
}

}  // namespace scattering

void xml_io_stream<scattering::ScatteringHabit>::write(std::ostream&,
                                                       const scattering::ScatteringHabit&,
                                                       bofstream*,
                                                       std::string_view) {
  throw std::runtime_error("private data not readable");
}

void xml_io_stream<scattering::ScatteringHabit>::read(std::istream&, scattering::ScatteringHabit&, bifstream*) {
  throw std::runtime_error("private data not writeable");
}

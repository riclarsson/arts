#include <arts_constants.h>
#include <physics_funcs.h>
#include <workspace.h>
#include <workspace_agenda_creator.h>

namespace {
Numeric water_mass_density(const AtmPoint& atm_point) {
  constexpr SpeciesEnum water = "H2O"_spec;

  ARTS_USER_ERROR_IF(not atm_point.has(water),
                     "Cannot derive water mass density: the atmospheric point has no H2O VMR. "
                     "Set H2O or provide water_mass_density explicitly.")

  return atm_point.number_density(water) * atm_point.mean_mass(water) * Constant::m_u;
}
}  // namespace

void single_dispersionAddWaterVisibleNIRHarvey98(Numeric&        single_dispersion,
                                                 const Numeric&  freq,
                                                 const AtmPoint& atm_point,
                                                 const Numeric&  water_mass_density_,
                                                 const Index&    check_validity) try {
  ARTS_TIME_REPORT

  const Numeric density = water_mass_density_ < 0.0 ? water_mass_density(atm_point) : water_mass_density_;
  single_dispersion +=
      refractive_index_water_visible_nir_harvey98(freq, atm_point.temperature, density, check_validity != 0) - 1.0;
}
ARTS_METHOD_ERROR_CATCH

void single_propmat_agendaSetWaterVisibleNIRHarvey98(Agenda&        single_propmat_agenda,
                                                     const Numeric& water_mass_density_,
                                                     const Index&   check_validity) {
  ARTS_TIME_REPORT

  AgendaCreator creator("single_propmat_agenda");
  creator.add("single_propmatInit");
  creator.add("single_dispersionAddWaterVisibleNIRHarvey98",
              SetWsv{"water_mass_density", water_mass_density_},
              SetWsv{"check_validity", check_validity});
  single_propmat_agenda = std::move(creator).finalize(true);
}

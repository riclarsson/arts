/* Copyright (C) 2020
 * Richard Larsson <ric.larsson@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 * USA. */

/**
 * @file m_linemixing.cc
 * @author Richard Larsson
 * @date 2020-06-23
 * 
 * @brief User interface for dealing with pure line mixing calculations
 * 
 * Note: If defined using parameterized form, the normal line-functions
 * approach is faster and more appropriate.  These functions should first
 * compute the relaxation, not simply use the relaxation
 */

#include "hitran_species.h"
#include "linemixing.h"
#include "linemixing_hitran.h"
#include "propagationmatrix.h"


void abs_hitran_relmat_dataReadHitranRelmatDataAndLines(HitranRelaxationMatrixData& abs_hitran_relmat_data,
                                                        ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                                        const ArrayOfArrayOfSpeciesTag& abs_species,
                                                        const String& basedir,
                                                        const Numeric& linemixinglimit,
                                                        const Numeric& fmin,
                                                        const Numeric& fmax,
                                                        const Numeric& stot,
                                                        const String& mode,
                                                        const String& hitran_type,
                                                        const Verbosity&)
{
  const lm_hitran_2017::ModeOfLineMixing intmode = lm_hitran_2017::toModeOfLineMixingOrThrow(mode);
  
  const SpeciesIsotopologueRatios isotopologue_ratios = Hitran::isotopologue_ratios(Hitran::toTypeOrThrow(hitran_type));
  
  ARTS_USER_ERROR_IF (abs_species.nelem() not_eq abs_lines_per_species.nelem(),
                      "Bad size of input species+lines");
  
  ArrayOfAbsorptionLines lines;
  lm_hitran_2017::read(abs_hitran_relmat_data, lines, isotopologue_ratios, basedir, linemixinglimit, fmin, fmax, stot, intmode);
  std::for_each(lines.begin(), lines.end(), [](auto& band){band.sort_by_frequency();});  // Sort so largest frequency is last
  ArrayOfIndex used(lines.nelem(), false);
  
  bool emptied=false;
  for (Index i=0; i<abs_species.nelem(); i++) {
    
    for (Index j=0; j<abs_species[i].nelem(); j++) {
      if (abs_species[i][j].Spec() not_eq Species::fromShortName("CO2")) 
        continue;
      
      if (not emptied) {
        abs_lines_per_species[i].resize(0);
        emptied = true;
      }
      
      for (Index k=0; k<lines.nelem(); k++) {
        if (used[k]) continue;
        
        const Numeric lf{abs_species[i][j].lower_freq > 0 ? abs_species[i][j].lower_freq : -std::numeric_limits<Numeric>::max()};
        const Numeric uf{abs_species[i][j].upper_freq > 0 ? abs_species[i][j].upper_freq : std::numeric_limits<Numeric>::max()};
        
        // Select lines with correct Isotopologue and one line center within the range
        if ((abs_species[i][j].Isotopologue() == lines[k].Isotopologue() or 
             abs_species[i][j].Isotopologue() == Species::select_joker("CO2")) and 
             (lines[k].AllLines().front().F0() <= uf and lines[k].AllLines().back().F0() >= lf)) {
          used[k] = true;  // The lines should not be copied into other places
          abs_lines_per_species[i].push_back(lines[k]);
        }
      }
    }
  }
}

void propmat_clearskyAddHitranLineMixingLines(PropagationMatrix& propmat_clearsky,
                                              const HitranRelaxationMatrixData& abs_hitran_relmat_data,
                                              const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                              const SpeciesIsotopologueRatios& isotopologue_ratios,
                                              const Vector& f_grid,
                                              const ArrayOfArrayOfSpeciesTag& abs_species,
                                              const ArrayOfRetrievalQuantity& jacobian_quantities,
                                              const Numeric& rtp_pressure,
                                              const Numeric& rtp_temperature,
                                              const Vector& rtp_vmr,
                                              const Verbosity&)
{
  ARTS_USER_ERROR_IF (jacobian_quantities.nelem(),
                      "Cannot support any Jacobian at this time");
  ARTS_USER_ERROR_IF (abs_species.nelem() not_eq abs_lines_per_species.nelem(),
                      "Bad size of input species+lines");
  ARTS_USER_ERROR_IF (abs_species.nelem() not_eq rtp_vmr.nelem(),
                      "Bad size of input species+vmrs");
  
  // vmrs should be [air, water, co2]
  Vector vmrs(3, 0);
  for (Index i=0; i<abs_species.nelem(); i++) {
    auto& specs = abs_species[i];
    for (auto& spec: specs) {
      if (Species::fromShortName("H2O") == spec.Spec()) {
        vmrs[1] = rtp_vmr[i];
      }
      else if (Species::fromShortName("CO2") == spec.Spec()) {
        vmrs[2] = rtp_vmr[i];
      }
    }
  }
  vmrs[0] = 1 - vmrs[1] - vmrs[2];
    
  for (Index i=0; i<abs_species.nelem(); i++) {
    if (abs_lines_per_species[i].nelem() and 
       (abs_lines_per_species[i].front().Population() == Absorption::PopulationType::ByHITRANFullRelmat or
        abs_lines_per_species[i].front().Population() == Absorption::PopulationType::ByHITRANRosenkranzRelmat))
      propmat_clearsky.Kjj() += lm_hitran_2017::compute(abs_hitran_relmat_data, abs_lines_per_species[i], isotopologue_ratios, rtp_pressure, rtp_temperature, vmrs, f_grid);
  }
}

void abs_lines_per_speciesAdaptHitranLineMixing(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                                const HitranRelaxationMatrixData& abs_hitran_relmat_data,
                                                const Vector& t_grid,
                                                const Numeric& pressure,
                                                const Index& order,
                                                const Verbosity&) 
{
  for (auto& abs_lines: abs_lines_per_species) {
    for (auto& band: abs_lines) {
      if (band.Population() == Absorption::PopulationType::ByHITRANFullRelmat or band.Population() == Absorption::PopulationType::ByHITRANRosenkranzRelmat) {
        lm_hitran_2017::hitran_lm_eigenvalue_adaptation(band, t_grid, abs_hitran_relmat_data, pressure, order);
      }
    }
  }
}

void abs_lines_per_speciesHitranLineMixingAdaptationData(ArrayOfTensor5& lm_data,
                                                         const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                                         const HitranRelaxationMatrixData& abs_hitran_relmat_data,
                                                         const Vector& t_grid,
                                                         const Vector& p_grid,
                                                         const Verbosity&) 
{
  lm_data.resize(0);
  for (auto& abs_lines: abs_lines_per_species) {
    for (auto& band: abs_lines) {
      if (band.Population() == Absorption::PopulationType::ByHITRANFullRelmat or band.Population() == Absorption::PopulationType::ByHITRANRosenkranzRelmat) {
        lm_data.emplace_back(lm_hitran_2017::hitran_lm_eigenvalue_adaptation_test(band, t_grid, abs_hitran_relmat_data, p_grid));
      }
    }
  }
}

void propmat_clearskyAddOnTheFlyLineMixing(PropagationMatrix& propmat_clearsky,
                                           ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
                                           const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                           const MapOfErrorCorrectedSuddenData& ecs_data,
                                           const SpeciesIsotopologueRatios& isotopologue_ratios,
                                           const Vector& f_grid,
                                           const ArrayOfArrayOfSpeciesTag& abs_species,
                                           const ArrayOfRetrievalQuantity& jacobian_quantities,
                                           const Numeric& rtp_pressure,
                                           const Numeric& rtp_temperature,
                                           const Vector& rtp_vmr,
                                           const Index& lbl_checked,
                                           const Verbosity&)
{
  ARTS_USER_ERROR_IF (abs_species.nelem() not_eq abs_lines_per_species.nelem(),
                      "Bad size of input species+lines");
  ARTS_USER_ERROR_IF (abs_species.nelem() not_eq rtp_vmr.nelem(),
                      "Bad size of input species+vmrs");
  ARTS_USER_ERROR_IF (not lbl_checked,
                      "Please set lbl_checked true to use this function");
  
  for (Index i=0; i<abs_species.nelem(); i++) {
    for (auto& band: abs_lines_per_species[i]) {
      if (band.OnTheFlyLineMixing() and band.DoLineMixing(rtp_pressure)) {
        // vmrs should be for the line
        const Vector line_shape_vmr = band.BroadeningSpeciesVMR(rtp_vmr, abs_species);
        const Numeric this_vmr = rtp_vmr[i] * isotopologue_ratios[band.Isotopologue()];
        const auto [abs, dabs] = Absorption::LineMixing::ecs_absorption(rtp_temperature,
                                                                        rtp_pressure,
                                                                        this_vmr,
                                                                        line_shape_vmr,
                                                                        ecs_data[band.QuantumIdentity()],
                                                                        f_grid,
                                                                        band,
                                                                        jacobian_quantities);
        propmat_clearsky.Kjj() += abs.real();
        
        // Sum up the resorted Jacobian
        for (Index j=0; j<jacobian_quantities.nelem(); j++) {
          const auto& deriv = jacobian_quantities[j];
          
          if (not deriv.propmattype()) continue;
          
          if (deriv == abs_species[i]) {
            dpropmat_clearsky_dx[j].Kjj() += abs.real();
          } else {
            dpropmat_clearsky_dx[j].Kjj() += dabs[j].real();
          }
        }
      }
    }
  }
}

void propmat_clearskyAddOnTheFlyLineMixingWithZeeman(PropagationMatrix& propmat_clearsky,
                                                     ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
                                                     const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                                     const MapOfErrorCorrectedSuddenData& ecs_data,
                                                     const SpeciesIsotopologueRatios& isotopologue_ratios,
                                                     const Vector& f_grid,
                                                     const ArrayOfArrayOfSpeciesTag& abs_species,
                                                     const ArrayOfRetrievalQuantity& jacobian_quantities,
                                                     const Numeric& rtp_pressure,
                                                     const Numeric& rtp_temperature,
                                                     const Vector& rtp_vmr,
                                                     const Vector& rtp_mag,
                                                     const Vector& rtp_los,
                                                     const Index& lbl_checked,
                                                     const Verbosity&)
{
  ARTS_USER_ERROR_IF (propmat_clearsky.StokesDimensions() not_eq 4,
                      "Only for stokes dim 4");
  ARTS_USER_ERROR_IF (abs_species.nelem() not_eq abs_lines_per_species.nelem(),
                      "Bad size of input species+lines");
  ARTS_USER_ERROR_IF (abs_species.nelem() not_eq rtp_vmr.nelem(),
                      "Bad size of input species+vmrs");
  ARTS_USER_ERROR_IF (not lbl_checked,
                      "Please set lbl_checked true to use this function");
  
  // Polarization
  const auto Z = Zeeman::FromGrids(rtp_mag[0], rtp_mag[1], rtp_mag[2], Conversion::deg2rad(rtp_los[0]), Conversion::deg2rad(rtp_los[1]));
  const auto polarization_scale_data = Zeeman::AllPolarization(Z.theta, Z.eta);
  const auto polarization_scale_dtheta_data = Zeeman::AllPolarization_dtheta(Z.theta, Z.eta);
  const auto polarization_scale_deta_data = Zeeman::AllPolarization_deta(Z.theta, Z.eta);
  
  for (Index i=0; i<abs_species.nelem(); i++) {
    for (auto& band: abs_lines_per_species[i]) {
      if (band.OnTheFlyLineMixing() and band.DoLineMixing(rtp_pressure)) {
        // vmrs should be for the line
        const Vector line_shape_vmr = band.BroadeningSpeciesVMR(rtp_vmr, abs_species);
        const Numeric this_vmr = rtp_vmr[i] * isotopologue_ratios[band.Isotopologue()];
        for (Zeeman::Polarization polarization : {Zeeman::Polarization::Pi, Zeeman::Polarization::SigmaMinus, Zeeman::Polarization::SigmaPlus}) {
          const auto [abs, dabs] = Absorption::LineMixing::ecs_absorption_zeeman(rtp_temperature,
                                                                                 Z.H,
                                                                                 rtp_pressure,
                                                                                 this_vmr,
                                                                                 line_shape_vmr,
                                                                                 ecs_data[band.QuantumIdentity()],
                                                                                 f_grid,
                                                                                 polarization,
                                                                                 band,
                                                                                 jacobian_quantities);
          
          // Sum up the propagation matrix
          Zeeman::sum(propmat_clearsky, abs, Zeeman::SelectPolarization(polarization_scale_data, polarization));
          
          // Sum up the resorted Jacobian
          for (Index j=0; j<jacobian_quantities.nelem(); j++) {
            const auto& deriv = jacobian_quantities[j];
            
            if (not deriv.propmattype()) continue;
            
            if (deriv == Jacobian::Atm::MagneticU) {
              Zeeman::dsum(dpropmat_clearsky_dx[j], abs, dabs[j],
                           Zeeman::SelectPolarization(polarization_scale_data, polarization),
                           Zeeman::SelectPolarization(polarization_scale_dtheta_data, polarization),
                           Zeeman::SelectPolarization(polarization_scale_deta_data, polarization),
                           Z.dH_du, Z.dtheta_du, Z.deta_du);
            } else if (deriv == Jacobian::Atm::MagneticV) {
              Zeeman::dsum(dpropmat_clearsky_dx[j], abs, dabs[j],
                           Zeeman::SelectPolarization(polarization_scale_data, polarization),
                           Zeeman::SelectPolarization(polarization_scale_dtheta_data, polarization),
                           Zeeman::SelectPolarization(polarization_scale_deta_data, polarization),
                           Z.dH_dv, Z.dtheta_dv, Z.deta_dv);
            } else if (deriv == Jacobian::Atm::MagneticW) {
              Zeeman::dsum(dpropmat_clearsky_dx[j], abs, dabs[j],
                           Zeeman::SelectPolarization(polarization_scale_data, polarization),
                           Zeeman::SelectPolarization(polarization_scale_dtheta_data, polarization),
                           Zeeman::SelectPolarization(polarization_scale_deta_data, polarization),
                           Z.dH_dw, Z.dtheta_dw, Z.deta_dw);
            } else if (deriv == abs_species[i]) {
              Zeeman::sum(dpropmat_clearsky_dx[j], abs, Zeeman::SelectPolarization(polarization_scale_data, polarization));
            } else {
              dpropmat_clearsky_dx[j].Kjj() += dabs[j].real();
            }
          }
        }
      }
    }
  }
}

void ecs_dataInit(MapOfErrorCorrectedSuddenData& ecs_data,
                  const Verbosity&)
{
  ecs_data.resize(0);
}

void ecs_dataSetSpeciesData(
  MapOfErrorCorrectedSuddenData& ecs_data,
  const SpeciesIsotopologueRatios& isotopologue_ratios,
  const Quantum::Identifier& qid,
  const String& species,
  const String& scaling_type,
  const Vector& scaling,
  const String& beta_type,
  const Vector& beta,
  const String& lambda_type,
  const Vector& lambda,
  const String& collisional_distance_type,
  const Vector& collisional_distance,
  const Verbosity&)
{
  const Species::Species spec = Species::fromShortName(species);
  ARTS_USER_ERROR_IF(not good_enum(spec), "Invalid species: ", species)
  auto& data = ecs_data[qid][spec];
  data.scaling = LineShapeModelParameters(LineShape::toTemperatureModelOrThrow(scaling_type), scaling);
  data.beta = LineShapeModelParameters(LineShape::toTemperatureModelOrThrow(beta_type), beta);
  data.lambda = LineShapeModelParameters(LineShape::toTemperatureModelOrThrow(lambda_type), lambda);
  data.collisional_distance = LineShapeModelParameters(LineShape::toTemperatureModelOrThrow(collisional_distance_type), collisional_distance);
  data.mass = Species::mean_mass(spec, isotopologue_ratios);
}

void ecs_dataSetMeanAir(
  MapOfErrorCorrectedSuddenData& ecs_data,
  const Vector& vmrs,
  const ArrayOfSpeciesTag& specs,
  const Verbosity&)
{
  ARTS_USER_ERROR_IF(specs.nelem() not_eq vmrs.nelem(), "Bad sizes of specs and vmrs\nspecs: [", specs, "]\nvmrs: [", vmrs, "]")
  
  static_assert(LineShapeModelParameters::N == 4);
  const auto scale = [](LineShapeModelParameters& mp, Numeric x){mp.X0 *= x; mp.X1 *= x; mp.X2 *= x; mp.X3 *= x;};
  const auto scale_and_add_or_throw = [](LineShapeModelParameters& mp1, const LineShapeModelParameters& mp2, Numeric x) {
    ARTS_USER_ERROR_IF(mp1.type not_eq mp2.type, "Can only scale and add same type\nmp1: ", mp1, "\nmp2: ", mp2)
    mp1.X0 += mp2.X0 * x;
    mp1.X1 += mp2.X1 * x;
    mp1.X2 += mp2.X2 * x;
    mp1.X3 += mp2.X3 * x;
  };
  
  for (auto& ecs: ecs_data) {
    SpeciesErrorCorrectedSuddenData airdata;
    
    bool found=false;
    Numeric sumvmr = 0;
    for (Index i=0; i<specs.nelem(); i++) {
      ARTS_USER_ERROR_IF(not specs[i].Isotopologue().joker(),
                         "Can only have joker species, finds: [", specs, ']')
      
      auto spec = specs[i].Spec();
      
      if (spec == Species::Species::Bath) continue;
      
      auto& data = ecs[spec];
      
      // This means we never had the data so we should skip
      if (std::isinf(data.mass)) {
        continue;
      }
      
      const Numeric vmr = vmrs[i];
      if (not found) {
        airdata.scaling = data.scaling;
        airdata.beta = data.beta;
        airdata.lambda = data.lambda;
        airdata.collisional_distance = data.collisional_distance;
        airdata.mass = data.mass * vmr;
        scale(airdata.scaling, vmr);
        scale(airdata.beta, vmr);
        scale(airdata.lambda, vmr);
        scale(airdata.collisional_distance, vmr);
        
        found = true;
      } else {
        scale_and_add_or_throw(airdata.scaling, data.scaling, vmr);
        scale_and_add_or_throw(airdata.beta, data.beta, vmr);
        scale_and_add_or_throw(airdata.lambda, data.lambda, vmr);
        scale_and_add_or_throw(airdata.collisional_distance, data.collisional_distance, vmr);
        airdata.mass += data.mass * vmr;
      }
      
      sumvmr += vmrs[i];
    }
    
    if (sumvmr not_eq 0 and sumvmr not_eq 1) {
      scale(airdata.scaling, 1.0 / sumvmr);
      scale(airdata.beta, 1.0 / sumvmr);
      scale(airdata.lambda, 1.0 / sumvmr);
      scale(airdata.collisional_distance, 1.0 / sumvmr);
      airdata.mass /= sumvmr;
    }
    
    ecs[Species::Species::Bath] = airdata;
  }
}

void ecs_dataAddMakarov2020(MapOfErrorCorrectedSuddenData& ecs_data,
                            const SpeciesIsotopologueRatios& isotopologue_ratios,
                            const Verbosity&)
{
  // The band is ignored
  auto& ecs = ecs_data[QuantumIdentifier("O2-66 ALL")];
  
  // All species have the same effect, so just copy the values but change the mass (allow new mass for Air)
  
  ecs[Species::Species::Oxygen].scaling = LineShapeModelParameters(LineShapeTemperatureModel::T0, 1.0, 0, 0, 0);
  ecs[Species::Species::Oxygen].collisional_distance = LineShapeModelParameters(LineShapeTemperatureModel::T0, Conversion::angstrom2meter(0.61), 0, 0, 0);
  ecs[Species::Species::Oxygen].lambda = LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.39, 0, 0, 0);
  ecs[Species::Species::Oxygen].beta = LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.567, 0, 0, 0);
  ecs[Species::Species::Oxygen].mass = Species::mean_mass(Species::Species::Oxygen, isotopologue_ratios);
  
  ecs[Species::Species::Nitrogen].scaling = LineShapeModelParameters(LineShapeTemperatureModel::T0, 1.0, 0, 0, 0);
  ecs[Species::Species::Nitrogen].collisional_distance = LineShapeModelParameters(LineShapeTemperatureModel::T0, Conversion::angstrom2meter(0.61), 0, 0, 0);
  ecs[Species::Species::Nitrogen].lambda = LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.39, 0, 0, 0);
  ecs[Species::Species::Nitrogen].beta = LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.567, 0, 0, 0);
  ecs[Species::Species::Nitrogen].mass = Species::mean_mass(Species::Species::Nitrogen, isotopologue_ratios);
}

void ecs_dataAddRodrigues1997(MapOfErrorCorrectedSuddenData& ecs_data,
                              const SpeciesIsotopologueRatios& isotopologue_ratios,
                              const Verbosity&)
{
  for (auto key: {"CO2-626 ALL", "CO2-628 ALL", "CO2-636 ALL"}) {
    auto& ecs = ecs_data[QuantumIdentifier(key)];
    
    ecs[Species::Species::Nitrogen].scaling = LineShapeModelParameters(LineShapeTemperatureModel::T1, Conversion::kaycm_per_atm2hz_per_pa(0.0180), 0.85, 0, 0);
    ecs[Species::Species::Nitrogen].lambda = LineShapeModelParameters(LineShapeTemperatureModel::T1, 0.81, 0.0152, 0, 0);
    ecs[Species::Species::Nitrogen].beta = LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.008, 0, 0, 0);
    ecs[Species::Species::Nitrogen].collisional_distance = LineShapeModelParameters(LineShapeTemperatureModel::T0, Conversion::angstrom2meter(2.2), 0, 0, 0);
    ecs[Species::Species::Nitrogen].mass = Species::mean_mass(Species::Species::Nitrogen, isotopologue_ratios);
    
    ecs[Species::Species::Oxygen].scaling = LineShapeModelParameters(LineShapeTemperatureModel::T1, Conversion::kaycm_per_atm2hz_per_pa(0.0168), 0.5, 0, 0);
    ecs[Species::Species::Oxygen].lambda = LineShapeModelParameters(LineShapeTemperatureModel::T1, 0.82, -0.091, 0, 0);
    ecs[Species::Species::Oxygen].beta = LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.007, 0, 0, 0);
    ecs[Species::Species::Oxygen].collisional_distance = LineShapeModelParameters(LineShapeTemperatureModel::T0, Conversion::angstrom2meter(2.4), 0, 0, 0);
    ecs[Species::Species::Oxygen].mass = Species::mean_mass(Species::Species::Oxygen, isotopologue_ratios);
  }
}

void ecs_dataAddTran2011(MapOfErrorCorrectedSuddenData& ecs_data,
                              const SpeciesIsotopologueRatios& isotopologue_ratios,
                              const Verbosity&)
{
  for (auto key: {"CO2-626 ALL", "CO2-628 ALL", "CO2-636 ALL"}) {
    auto& ecs = ecs_data[QuantumIdentifier(key)];
    
    ecs[Species::Species::CarbonDioxide].scaling = LineShapeModelParameters(LineShapeTemperatureModel::T0, Conversion::kaycm_per_atm2hz_per_pa(0.019), 0, 0, 0);
    ecs[Species::Species::CarbonDioxide].lambda = LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.61, 0, 0, 0);
    ecs[Species::Species::CarbonDioxide].beta = LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.052, 0, 0, 0);
    ecs[Species::Species::CarbonDioxide].collisional_distance = LineShapeModelParameters(LineShapeTemperatureModel::T0, Conversion::angstrom2meter(5.5), 0, 0, 0);
    ecs[Species::Species::CarbonDioxide].mass = Species::mean_mass(Species::Species::CarbonDioxide, isotopologue_ratios);
  }
}

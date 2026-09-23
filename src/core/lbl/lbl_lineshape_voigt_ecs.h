#pragma once

#include <array.h>
#include <rtepack.h>

#include "lbl_data.h"
#include "lbl_lineshape_linemixing.h"

//! FIXME: These functions should be elsewhere?
namespace Jacobian {
struct Targets;
}  // namespace Jacobian

namespace lbl::voigt::ecs {
//! Rotational quantum numbers of a line, prepared in matrix order.
//! N is only used by the Makarov kernel.
struct rotational_line {
  Rational Ju{}, Jl{}, Nu{}, Nl{};
};

//! Energies [J] prepared from one species model before evaluating collision partners.
struct energy_data {
  //! Resolved lower-state energies, in matrix line order.
  Vector e0{};
  //! Reference-rotor energies E(L) and E(L-2), indexed by angular momentum.
  //! These also supply the line-side adiabatic factors at J (CO2) or N (O2).
  Vector rotational{}, rotational_minus_two{};
};

//! Fill the reference ladder from the species energy model, independently of line order.
void prepare_rotational_ladder(energy_data& energies, int count, Numeric (*energy)(Rational));

//! Collision-partner basis rates, indexed by angular momentum; Q[0] is unused and zero.
struct basis_data {
  Vector Q{}, Omega{};
};

//! Evaluate and validate the common ECS basis from the prepared reference energies.
basis_data prepare_basis(int                             count,
                         const energy_data&              energies,
                         const linemixing::species_data& collision,
                         Numeric                         T0,
                         const SpeciesIsotope&           isot,
                         SpeciesEnum                     broadener,
                         const AtmPoint&                 atm);

//! Sequential truncated-band correction; vectors follow the matrix ordering.
void apply_sum_rule(MatrixView W, ConstVectorView dipr, ConstVectorView e0, Numeric T);

struct ComputeData {
  Numeric gd_fac{};  //! Gaussian 1/e half-width divided by line frequency

  //! Size of line shapes
  Vector       pop{};
  Vector       dip{};
  Vector       dipr{};
  ArrayOfIndex sort{};

  //! Collision energies share a species model; optical pop retains catalogue e0.
  energy_data                  energies{};
  std::vector<rotational_line> rotational_lines{};

  //! Size of line shapes x size of line shapes
  Matrix Wimag{};

  //! [1, or broadening species]
  Vector vmrs;

  //! [1, or broadening species] x size of line shapes
  ComplexMatrix eqv_strs{};
  ComplexMatrix eqv_vals{};

  //! [1, or broadening species] x size of line shapes x size of line shapes
  ComplexTensor3 Ws{};
  ComplexTensor3 Vs{};

  //! Reciprocal condition number of each equivalent-line eigenvector matrix.
  Vector eigenvector_rcond{};

  //! Per-broadener max |sum_j dipr[j] W[j,i]| / sum_j |dipr[j] W[j,i]|.
  //! A finite truncated band need not satisfy the optical sum rule exactly.
  Vector sum_rule_residual{};

  //! Size of frequency
  Vector        scl{};
  ComplexVector shape{};

  //! The orientation of the polarization
  Propmat npm{};

  //! Sizes scl, dscl, shape, dshape.  Sets scl, npm, dnpm_du, dnpm_dv, dnpm_dw
  ComputeData(const ConstVectorView&   f_grid,
              const AtmPoint&          atm,
              const Vector2&           los = {},
              const ZeemanPolarization pol = ZeemanPolarization::no);

  void update_zeeman(const Vector2& los, const Vector3& mag, const ZeemanPolarization pol);

  void core_calc_eqv();
  void core_calc(const ConstVectorView& f_grid);
  void adapt_single(const QuantumIdentifier&        bnd_qid,
                    const band_data&                bnd,
                    const LinemixingSpeciesEcsData& rovib_data,
                    const AtmPoint&                 atm,
                    const bool                      presorted = false);
  void adapt_multi(const QuantumIdentifier&        bnd_qid,
                   const band_data&                bnd,
                   const LinemixingSpeciesEcsData& rovib_data,
                   const AtmPoint&                 atm,
                   const bool                      presorted = false);

 private:
  void adapt(const QuantumIdentifier&        bnd_qid,
             const band_data&                bnd,
             const LinemixingSpeciesEcsData& rovib_data,
             const AtmPoint&                 atm,
             bool                            presorted,
             bool                            per_broadener);
};

void calculate(PropmatVectorView               pm,
               PropmatMatrixView               dpm,
               ComputeData&                    com_data,
               const ConstVectorView           f_grid,
               const Range&                    f_range,
               const Jacobian::Targets&        jac_targets,
               const QuantumIdentifier&        bnd_qid,
               const band_data&                bnd,
               const LinemixingSpeciesEcsData& rovib_data,
               const AtmPoint&                 atm,
               const ZeemanPolarization        pol,
               const bool                      no_negative_absorption);

void equivalent_values(ComplexTensor3View              eqv_str,
                       ComplexTensor3View              eqv_val,
                       ComputeData&                    com_data,
                       const QuantumIdentifier&        bnd_qid,
                       const band_data&                bnd,
                       const LinemixingSpeciesEcsData& rovib_data,
                       const AtmPoint&                 atm,
                       const Vector&                   T);
}  // namespace lbl::voigt::ecs

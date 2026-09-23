#pragma once

#include <matpack.h>

#include <span>

#include "lbl_data.h"
#include "lbl_lineshape_linemixing.h"

namespace lbl::voigt::ecs {
struct rotational_line;
struct energy_data;
}  // namespace lbl::voigt::ecs

namespace lbl::voigt::ecs::hartmann {
Numeric reduced_dipole(
    const Rational Jf, const Rational Ji, const Rational lf, const Rational li, const Rational k = Rational{1});

//! CO2-626 reference-rotor and resolved-state energies [J].
Numeric rotational_energy(Rational J);
Numeric level_energy(Rational J);

//! Prepare lower-state and reference-rotor energies before any angular-kernel swap.
void prepare_energies(energy_data& energies, const QuantumIdentifier& qid, std::span<const rotational_line> lines);

//! All per-line inputs follow the matrix ordering.
void relaxation_matrix_offdiagonal(MatrixView&                      W,
                                   const QuantumIdentifier&         bnd_qid,
                                   std::span<const rotational_line> lines,
                                   Numeric                          T0,
                                   const SpeciesEnum                broadening_species,
                                   const linemixing::species_data&  rovib_data,
                                   const Vector&                    dipr,
                                   const energy_data&               energies,
                                   const AtmPoint&                  atm);
}  // namespace lbl::voigt::ecs::hartmann

#pragma once

#include <atm.h>
#include <lbl.h>
#include <quantum_numbers.h>

namespace nonlte {
using data_t = std::unordered_map<QuantumIdentifier, AtmData>;

data_t level_density(const AtmField& atm,
                     const AbsorptionBands& bands);
}  // namespace nonlte

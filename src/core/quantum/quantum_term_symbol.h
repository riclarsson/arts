#include <string>

#include "quantum.h"

namespace Quantum::Helpers {
/** Returns the Molecular Term in LaTeX formatting of a QuantumIdentifier */
std::string molecular_term_symbol(const QuantumIdentifier& qid);
}  // namespace Quantum::Helpers

#pragma once

#include <matpack.h>

#include <variant>

namespace Atm::interp {
using namespace my_interp;

using altlag1 = Lagrange<1>;
using altlag0 = Lagrange<0>;

using latlag1 = Lagrange<1>;
using latlag0 = Lagrange<0>;

using lonlag1 = Lagrange<1, false, GridType::Cyclic, cycle_m180_p180>;
using lonlag0 = Lagrange<0, false, GridType::Cyclic, cycle_m180_p180>;

using altlags = std::variant<altlag0, altlag1>;
using latlags = std::variant<latlag0, latlag1>;
using lonlags = std::variant<lonlag0, lonlag1>;

altlags altlag(const AscendingGrid& xs, Numeric x);
latlags latlag(const AscendingGrid& xs, Numeric x);
lonlags lonlag(const AscendingGrid& xs, Numeric x);

Tensor3 interpweights(const altlags&, const latlags&, const lonlags&);
Matrix interpweights(const latlags&, const lonlags&);

Numeric get(const SortedGriddedField3& f,
            const Numeric,
            const Numeric,
            const Numeric);
Numeric get(const SortedGriddedField3& f,
            const Tensor3&,
            const altlags&,
            const latlags&,
            const lonlags&);
Numeric get(const SortedGriddedField3& f,
            Index ialt,
            const Matrix&,
            const latlags&,
            const lonlags&);

// Numeric interpolation
Numeric get(const Numeric, const Numeric, const Numeric, const Numeric);

// Functional interpolation
Numeric get(const std::function<Numeric(Numeric, Numeric, Numeric)>&,
            const Numeric,
            const Numeric,
            const Numeric);

std::vector<std::pair<Index, Numeric>> flat_weight(
    const SortedGriddedField3& gf3,
    const Numeric alt,
    const Numeric lat,
    const Numeric lon);
}  // namespace Atm::interp

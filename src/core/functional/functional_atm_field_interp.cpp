#include "functional_atm_field_interp.h"

namespace Atm::interp {
altlags altlag(const AscendingGrid& xs, Numeric x) {
  return xs.size() == 1 ? altlags{altlag0(0, x, xs)}
                        : altlags{altlag1(0, x, xs)};
}

latlags latlag(const AscendingGrid& xs, Numeric x) {
  return xs.size() == 1 ? latlags{latlag0(0, x, xs)}
                        : latlags{latlag1(0, x, xs)};
}

lonlags lonlag(const AscendingGrid& xs, Numeric x) {
  return xs.size() == 1 ? lonlags{lonlag0(0, x, xs)}
                        : lonlags{lonlag1(0, x, xs)};
}

Tensor3 interpweights(const altlags& a, const latlags& b, const lonlags& c) {
  return std::visit(
      [](auto& alt, auto& lat, auto& lon) -> Tensor3 {
        auto d = my_interp::interpweights(alt, lat, lon);
        Tensor3 out(alt.size(), lat.size(), lon.size());
        for (Size i0 = 0; i0 < alt.size(); i0++)
          for (Size i1 = 0; i1 < lat.size(); i1++)
            for (Size i2 = 0; i2 < lon.size(); i2++)
              out[i0, i1, i2] = d[i0, i1, i2];
        return out;
      },
      a,
      b,
      c);
}

Matrix interpweights(const latlags& b, const lonlags& c) {
  return std::visit(
      [](auto& lat, auto& lon) -> Matrix {
        auto d = my_interp::interpweights(lat, lon);
        Matrix out(lat.size(), lon.size());
        for (Size i1 = 0; i1 < lat.size(); i1++)
          for (Size i2 = 0; i2 < lon.size(); i2++) out[i1, i2] = d[i1, i2];
        return out;
      },
      b,
      c);
}

Numeric get(const SortedGriddedField3& f,
            const Numeric alt,
            const Numeric lat,
            const Numeric lon) {
  if (not f.ok()) throw std::runtime_error("bad field");
  return std::visit(
      [&data = f.data](auto&& al, auto&& la, auto&& lo) {
        return my_interp::interp(data, al, la, lo);
      },
      altlag(f.grid<0>(), alt),
      latlag(f.grid<1>(), lat),
      lonlag(f.grid<2>(), lon));
}

Numeric get(const SortedGriddedField3& f,
            const Tensor3& iw,
            const altlags& alt,
            const latlags& lat,
            const lonlags& lon) {
  if (not f.ok()) throw std::runtime_error("bad field");
  return std::visit(
      [&data = f.data, &iw](auto&& al, auto&& la, auto&& lo) {
        return my_interp::interp(data, iw, al, la, lo);
      },
      alt,
      lat,
      lon);
}

Numeric get(const SortedGriddedField3& f,
            Index ialt,
            const Matrix& iw,
            const latlags& lat,
            const lonlags& lon) {
  if (not f.ok()) throw std::runtime_error("bad field");
  return std::visit(
      [data = f.data[ialt], &iw](auto&& la, auto&& lo) {
        return my_interp::interp(data, iw, la, lo);
      },
      lat,
      lon);
}

Numeric get(const Numeric num, const Numeric, const Numeric, const Numeric) {
  return num;
}

Numeric get(const std::function<Numeric(Numeric, Numeric, Numeric)>& fd,
            const Numeric alt,
            const Numeric lat,
            const Numeric lon) {
  return fd(alt, lat, lon);
}

std::vector<std::pair<Index, Numeric>> flat_weight(
    const SortedGriddedField3 &gf3,
    const Numeric alt,
    const Numeric lat,
    const Numeric lon) {
  if (not gf3.ok()) throw std::runtime_error("bad field");

  const Index nalt = gf3.grid<0>().size();
  const Index nlat = gf3.grid<1>().size();
  const Index nlon = gf3.grid<2>().size();

  return std::visit(
      [NN = nlat * nlon, N = nlon]<typename ALT, typename LAT, typename LON>(
          const ALT &al, const LAT &la, const LON &lo) {
        const auto x = interpweights(al, la, lo);
        std::vector<std::pair<Index, Numeric>> out;
        out.reserve(al.size() * NN);

        for (Index i = 0; i < al.size(); i++) {
          for (Index j = 0; j < la.size(); j++) {
            for (Index k = 0; k < lo.size(); k++) {
              out.emplace_back(
                  (al.pos + i) * NN + (la.pos + j) * N + lo.pos + k,
                  x[i, j, k]);
            }
          }
        }
        return out;
      },
      nalt == 1 ? altlags{gf3.lag<0, altlag0>(alt)}
                : altlags{gf3.lag<0, altlag1>(alt)},
      nlat == 1 ? latlags{gf3.lag<1, latlag0>(lat)}
                : latlags{gf3.lag<1, latlag1>(lat)},
      nlon == 1 ? lonlags{gf3.lag<2, lonlag0>(lon)}
                : lonlags{gf3.lag<2, lonlag1>(lon)});
}
}  // namespace Atm::interp

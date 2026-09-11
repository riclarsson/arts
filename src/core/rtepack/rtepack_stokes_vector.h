#pragma once

#include <enumsPolarizationChoice.h>
#include <enumsSpectralRadianceUnitType.h>
#include <matpack.h>
#include <xml.h>

#include <utility>

#include "rtepack_common.h"

namespace rtepack {
struct stokvec final : Vector4 {
  explicit constexpr stokvec(Numeric i) : Vector4{i, 0., 0., 0.} {}

  constexpr stokvec() : stokvec{0.0} {}

  constexpr stokvec(Numeric i, Numeric q, Numeric u, Numeric v) : Vector4{i, q, u, v} {}

  template <stdr::forward_range T> requires std::is_convertible_v<stdr::range_value_t<T>, Numeric>
  explicit constexpr stokvec(const T &input) noexcept : stokvec{} {
    assert(stdr::size(input) == 4);
    stdr::copy(input, data.begin());
  }

  template <stdr::forward_range T> requires std::is_convertible_v<stdr::range_value_t<T>, Numeric>
  constexpr stokvec &operator=(const T &input) noexcept {
    return *this = stokvec{input};
  }

  [[nodiscard]] constexpr Numeric I() const { return data[0]; }
  [[nodiscard]] constexpr Numeric Q() const { return data[1]; }
  [[nodiscard]] constexpr Numeric U() const { return data[2]; }
  [[nodiscard]] constexpr Numeric V() const { return data[3]; }

  [[nodiscard]] constexpr Numeric &I() { return data[0]; }
  [[nodiscard]] constexpr Numeric &Q() { return data[1]; }
  [[nodiscard]] constexpr Numeric &U() { return data[2]; }
  [[nodiscard]] constexpr Numeric &V() { return data[3]; }

  [[nodiscard]] constexpr bool is_zero() const { return I() == 0.0 && Q() == 0.0 && U() == 0.0 && V() == 0.0; }

  [[nodiscard]] constexpr bool is_polarized() const { return Q() != 0.0 or U() != 0.0 or V() != 0.0; }

  [[nodiscard]] constexpr Size nonzero_components() const {
    return static_cast<Size>(I() != 0.0) + static_cast<Size>(Q() != 0.0) + static_cast<Size>(U() != 0.0) +
           static_cast<Size>(V() != 0.0);
  }

  constexpr stokvec operator-() const { return stokvec{-data[0], -data[1], -data[2], -data[3]}; }

  constexpr stokvec &operator+=(const stokvec &b) {
    for (Size i = 0; i < 4; ++i) { data[i] += b.data[i]; }
    return *this;
  }

  constexpr stokvec &operator-=(const stokvec &b) {
    for (Size i = 0; i < 4; ++i) { data[i] -= b.data[i]; }
    return *this;
  }

  constexpr stokvec &operator/=(const Numeric &b) {
    for (auto &x : data) { x /= b; }
    return *this;
  }

  constexpr stokvec &operator*=(const Numeric &b) {
    for (auto &x : data) { x *= b; }
    return *this;
  }

  constexpr stokvec &operator*=(const stokvec &b) {
    for (Size i = 0; i < 4; ++i) { data[i] *= b.data[i]; }
    return *this;
  }
};

constexpr stokvec to_stokvec(PolarizationChoice p) {
  using enum PolarizationChoice;
  switch (p) {
    case I:    return {1.0, 0.0, 0.0, 0.0};
    case Q:    return {0.0, 1.0, 0.0, 0.0};
    case U:    return {0.0, 0.0, 1.0, 0.0};
    case V:    return {0.0, 0.0, 0.0, 1.0};
    case Iv:   return {1.0, 1.0, 0.0, 0.0};
    case Ih:   return {1.0, -1.0, 0.0, 0.0};
    case Ip45: return {1.0, 0.0, 1.0, 0.0};
    case Im45: return {1.0, 0.0, -1.0, 0.0};
    case Ilhc: return {1.0, 0.0, 0.0, -1.0};
    case Irhc: return {1.0, 0.0, 0.0, 1.0};
  }
  std::unreachable();
}

constexpr stokvec operator+(stokvec a, const stokvec &b) { return a += b; }
constexpr stokvec operator-(stokvec a, const stokvec &b) { return a -= b; }
constexpr stokvec operator*(const Numeric &a, stokvec b) { return b *= a; }
constexpr stokvec operator*(stokvec a, const Numeric &b) { return a *= b; }
constexpr stokvec operator/(stokvec a, const Numeric &b) { return a /= b; }

constexpr stokvec fma(const Numeric &x, const stokvec &a, const stokvec &b) {
  return {std::fma(x, a.I(), b.I()), std::fma(x, a.Q(), b.Q()), std::fma(x, a.U(), b.U()), std::fma(x, a.V(), b.V())};
}

//! Take the average of two stokvec vectors
constexpr stokvec avg(const stokvec &a, const stokvec &b) {
  return {std::midpoint(a.I(), b.I()),
          std::midpoint(a.Q(), b.Q()),
          std::midpoint(a.U(), b.U()),
          std::midpoint(a.V(), b.V())};
}

/** Convertion methods for stokes vectors
 * 
  * @param type Spectral Radiance Unit Type
 * @return A function that takes a stokes vector and a frequency and returns a converted stokes vector
 */
std::function<stokvec(const stokvec, const Numeric)>                unit_converter(const SpectralRadianceUnitType type,
                                                                                   const Numeric                  n = 1.0);
std::function<stokvec(const stokvec, const stokvec, const Numeric)> dunit_converter(const SpectralRadianceUnitType type,
                                                                                    const Numeric n = 1.0);

using stokvec_vector            = matpack::data_t<stokvec, 1>;
using stokvec_vector_view       = matpack::view_t<stokvec, 1>;
using stokvec_vector_const_view = matpack::view_t<const stokvec, 1>;

using stokvec_matrix            = matpack::data_t<stokvec, 2>;
using stokvec_matrix_view       = matpack::view_t<stokvec, 2>;
using stokvec_matrix_const_view = matpack::view_t<const stokvec, 2>;

using stokvec_tensor3            = matpack::data_t<stokvec, 3>;
using stokvec_tensor3_view       = matpack::view_t<stokvec, 3>;
using stokvec_tensor3_const_view = matpack::view_t<const stokvec, 3>;

using stokvec_tensor4            = matpack::data_t<stokvec, 4>;
using stokvec_tensor4_view       = matpack::view_t<stokvec, 4>;
using stokvec_tensor4_const_view = matpack::view_t<const stokvec, 4>;

using stokvec_tensor5            = matpack::data_t<stokvec, 5>;
using stokvec_tensor5_view       = matpack::view_t<stokvec, 5>;
using stokvec_tensor5_const_view = matpack::view_t<const stokvec, 5>;

using stokvec_tensor6            = matpack::data_t<stokvec, 6>;
using stokvec_tensor6_view       = matpack::view_t<stokvec, 6>;
using stokvec_tensor6_const_view = matpack::view_t<const stokvec, 6>;
}  // namespace rtepack

template <> struct std::formatter<rtepack::stokvec> : std::formatter<Vector4> {};

template <> struct xml_io_stream<rtepack::stokvec> : xml_io_stream_inherit<Vector4, rtepack::stokvec> {};

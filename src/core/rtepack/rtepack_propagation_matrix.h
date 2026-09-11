#pragma once

#include <compare.h>
#include <matpack.h>
#include <xml.h>

#include <numeric>

#include "rtepack_common.h"

namespace rtepack {

struct propmat final : Vector7 {
  explicit constexpr propmat(Numeric a) : Vector7{a, 0., 0., 0., 0., 0., 0.} {}

  constexpr propmat() : propmat(0.0) {}

  constexpr propmat(Numeric a, Numeric b, Numeric c, Numeric d, Numeric u, Numeric v, Numeric w)
      : Vector7{a, b, c, d, u, v, w} {}

  template <matpack::exact_md<Numeric, 2> T> explicit constexpr propmat(const T &input) noexcept {
    assert(input.extent(0) == 4 and input.extent(1) == 4);
    assert((input[1, 0] == 0.0 and input[2, 0] == 0.0 and input[3, 0] == 0.0));
    assert((input[1, 0] == input[0, 1] and input[2, 0] == input[0, 2] and input[3, 0] == input[0, 3]));
    assert((input[2, 1] == -input[1, 2] and input[3, 1] == -input[1, 3] and input[3, 2] == -input[2, 3]));
    data[0] = input[0, 0];
    data[1] = input[0, 1];
    data[2] = input[0, 2];
    data[3] = input[0, 3];
    data[4] = input[1, 2];
    data[5] = input[1, 3];
    data[6] = input[2, 3];
  }

  template <stdr::forward_range T> requires std::is_convertible_v<stdr::range_value_t<T>, Numeric>
  explicit constexpr propmat(const T &input) noexcept : propmat{} {
    assert(stdr::size(input) == 7);
    stdr::copy(input, data.begin());
  }

  template <matpack::exact_md<Numeric, 2> T> constexpr propmat &operator=(const T &input) noexcept {
    return *this = propmat{input};
  }

  template <stdr::forward_range T> requires std::is_convertible_v<stdr::range_value_t<T>, Numeric>
  constexpr propmat &operator=(const T &input) noexcept {
    return *this = propmat{input};
  }

  [[nodiscard]] constexpr decltype(auto) A() const { return data[0]; }
  [[nodiscard]] constexpr decltype(auto) B() const { return data[1]; }
  [[nodiscard]] constexpr decltype(auto) C() const { return data[2]; }
  [[nodiscard]] constexpr decltype(auto) D() const { return data[3]; }
  [[nodiscard]] constexpr decltype(auto) U() const { return data[4]; }
  [[nodiscard]] constexpr decltype(auto) V() const { return data[5]; }
  [[nodiscard]] constexpr decltype(auto) W() const { return data[6]; }

  [[nodiscard]] constexpr decltype(auto) A() { return data[0]; }
  [[nodiscard]] constexpr decltype(auto) B() { return data[1]; }
  [[nodiscard]] constexpr decltype(auto) C() { return data[2]; }
  [[nodiscard]] constexpr decltype(auto) D() { return data[3]; }
  [[nodiscard]] constexpr decltype(auto) U() { return data[4]; }
  [[nodiscard]] constexpr decltype(auto) V() { return data[5]; }
  [[nodiscard]] constexpr decltype(auto) W() { return data[6]; }

  //! Check if the matrix is purely rotational
  [[nodiscard]] constexpr bool is_rotational() const { return A() == 0.0 and B() == 0.0 and C() == 0.0 and D() == 0.0; }

  //! Check if the matrix is polarized
  [[nodiscard]] constexpr bool is_polarized() const {
    return B() != 0 or C() != 0 or D() != 0 or U() != 0 or V() != 0 or W() != 0;
  }

  constexpr auto operator<=>(const propmat &pm) const { return A() <=> pm.A(); }

  [[nodiscard]] constexpr propmat operator-() const { return propmat{-A(), -B(), -C(), -D(), -U(), -V(), -W()}; }

  constexpr propmat &operator+=(const propmat &b) {
    for (Size i = 0; i < 7; ++i) { data[i] += b.data[i]; }
    return *this;
  }

  constexpr propmat &operator-=(const propmat &b) {
    for (Size i = 0; i < 7; ++i) { data[i] -= b.data[i]; }
    return *this;
  }

  constexpr propmat &operator/=(const Numeric &b) {
    for (auto &x : data) { x /= b; }
    return *this;
  }

  constexpr propmat &operator*=(const Numeric &b) {
    for (auto &x : data) { x *= b; }
    return *this;
  }
};

constexpr propmat operator+(propmat a, const propmat &b) { return a += b; }
constexpr propmat operator-(propmat a, const propmat &b) { return a -= b; }
constexpr propmat operator*(propmat a, const Numeric &b) { return a *= b; }
constexpr propmat operator*(const Numeric &a, propmat b) { return b *= a; }
constexpr propmat operator/(propmat a, const Numeric &b) { return a /= b; }

constexpr Numeric det(const propmat &k) {
  const auto &[a, b, c, d, u, v, w] = k.data;

  const Numeric a2 = a * a;
  const Numeric b2 = b * b;
  const Numeric c2 = c * c;
  const Numeric d2 = d * d;
  const Numeric u2 = u * u;
  const Numeric v2 = v * v;
  const Numeric w2 = w * w;

  const Numeric C = b * w - c * v + d * u;

  return a2 * (a2 - b2 - c2 - d2 + u2 + v2 + w2) - C * C;
}

//! Take the average of two propmat matrixes
constexpr propmat avg(const propmat &a, const propmat &b) {
  return {std::midpoint(a.A(), b.A()),
          std::midpoint(a.B(), b.B()),
          std::midpoint(a.C(), b.C()),
          std::midpoint(a.D(), b.D()),
          std::midpoint(a.U(), b.U()),
          std::midpoint(a.V(), b.V()),
          std::midpoint(a.W(), b.W())};
}

//! Element-wise dawson function of a propmat matrix (FIXME: implement as matrix?)
propmat dawson(const propmat &pm);

using propmat_vector            = matpack::data_t<propmat, 1>;
using propmat_vector_view       = matpack::view_t<propmat, 1>;
using propmat_vector_const_view = matpack::view_t<const propmat, 1>;

using propmat_matrix            = matpack::data_t<propmat, 2>;
using propmat_matrix_view       = matpack::view_t<propmat, 2>;
using propmat_matrix_const_view = matpack::view_t<const propmat, 2>;

using propmat_tensor3            = matpack::data_t<propmat, 3>;
using propmat_tensor3_view       = matpack::view_t<propmat, 3>;
using propmat_tensor3_const_view = matpack::view_t<const propmat, 3>;

propmat_vector operator*(Numeric x, const propmat_vector_const_view &y);
}  // namespace rtepack

template <> struct std::formatter<rtepack::propmat> : std::formatter<Vector7> {};

template <> struct xml_io_stream<rtepack::propmat> : xml_io_stream_inherit<Vector7, rtepack::propmat> {};

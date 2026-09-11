#include <matpack.h>

#include <concepts>

namespace {
struct Value {
  Numeric x{};

  constexpr Value() = default;
  explicit constexpr Value(Numeric value) : x(value) {}

  constexpr Value& operator+=(Numeric y) {
    x += y;
    return *this;
  }
  constexpr Value& operator-=(Numeric y) {
    x -= y;
    return *this;
  }
  constexpr Value& operator*=(Numeric y) {
    x *= y;
    return *this;
  }
  constexpr Value& operator/=(Numeric y) {
    x /= y;
    return *this;
  }
};

static_assert(not std::convertible_to<Numeric, Value>);
}  // namespace

int main() {
  matpack::data_t<Value, 2> values(2, 3);

  values  = 8.0;
  values += 2.0;
  values -= 1.0;
  values *= 2.0;
  values /= 3.0;
  for (const auto& value : values | by_elem)
    if (value.x != 6.0) return 1;

  auto strided  = values[joker, StridedRange(0, 2, 2)];
  strided      /= 2.0;
  for (Size row = 0; row < 2; ++row) {
    if (values[row, 0].x != 3.0) return 1;
    if (values[row, 1].x != 6.0) return 1;
    if (values[row, 2].x != 3.0) return 1;
  }
}

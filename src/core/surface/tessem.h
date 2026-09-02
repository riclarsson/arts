#pragma once

#include <matpack.h>
#include <xml.h>

struct TessemNN {
  Index  nb_inputs{};
  Index  nb_outputs{};
  Index  nb_cache{};
  Vector b1;
  Vector b2;
  Matrix w1;
  Matrix w2;
  Vector x_min;
  Vector x_max;
  Vector y_min;
  Vector y_max;
};

void   tessem_read_ascii(const String& filename, TessemNN& net);
void   tessem_read_ascii(std::istream& is, TessemNN& net);
Vector tessem_emissivity(const TessemNN& net, ConstVectorView input);

template <> struct xml_io_stream_name<TessemNN> {
  static constexpr std::string_view name = "TessemNN";
};

template <> struct xml_io_stream_aggregate<TessemNN> {
  static constexpr bool value = true;
};

template <> struct std::formatter<TessemNN> {
  format_tags                   tags;
  constexpr auto&               inner_fmt() { return *this; }
  constexpr const auto&         inner_fmt() const { return *this; }
  constexpr auto                parse(std::format_parse_context& ctx) { return parse_format_tags(tags, ctx); }
  template <class Context> auto format(const TessemNN& x, Context& ctx) const {
    return tags.format(
        ctx, "TessemNN(inputs="sv, x.nb_inputs, ", hidden="sv, x.nb_cache, ", outputs="sv, x.nb_outputs, ")"sv);
  }
};

#pragma once

#include <matpack.h>
#include <xml.h>

struct TelsemAtlas {
  Index                       ndat{};
  static constexpr Index      nchan{7};
  String                      name;
  Index                       month{};
  Numeric                     dlat{};
  IndexVector                 ncells;
  IndexVector                 firstcells;
  matpack::data_t<Vector7, 1> channel_emissivity;
  matpack::data_t<Vector7, 1> channel_emissivity_error;
  IndexVector                 classes1;
  IndexVector                 classes2;
  IndexVector                 cellnums;
  IndexVector                 correspondence;

  [[nodiscard]] bool                        contains(Index cellnumber) const;
  [[nodiscard]] Index                       calc_cellnum(Numeric lat, Numeric lon) const;
  [[nodiscard]] Index                       calc_cellnum_nearest_neighbor(Numeric lat, Numeric lon) const;
  [[nodiscard]] std::pair<Numeric, Numeric> get_coordinates(Index cellnum) const;
  [[nodiscard]] Vector2                     emissivity(
      Numeric lat, Numeric lon, Numeric incidence_angle, Numeric frequency, Numeric max_distance = -1) const;
  [[nodiscard]] Vector2 emis_interp(
      Numeric theta, Numeric freq_ghz, Index class1, Index class2, const Vector3& ev, const Vector3& eh) const;
  [[nodiscard]] Vector3 get_emis_v(Index cellnumber) const;
  [[nodiscard]] Vector3 get_emis_h(Index cellnumber) const;

  void read(std::istream& is);
  void equare();
  void rebuild_correspondence();

 private:
  [[nodiscard]] Numeric interp_freq2(
      Numeric emiss19, Numeric emiss37, Numeric emiss85, Numeric freq_ghz, Index class2) const;
};

void telsem_read_ascii(const String& filename, TelsemAtlas& atlas, Index month = 0);
void telsem_read_ascii(std::istream& is, TelsemAtlas& atlas, Index month = 0);

template <> struct xml_io_stream_name<TelsemAtlas> {
  static constexpr std::string_view name = "TelsemAtlas";
};

template <> struct xml_io_stream_aggregate<TelsemAtlas> {
  static constexpr bool value = true;
};

template <> struct std::formatter<TelsemAtlas> {
  format_tags                   tags;
  constexpr auto&               inner_fmt() { return *this; }
  constexpr const auto&         inner_fmt() const { return *this; }
  constexpr auto                parse(std::format_parse_context& ctx) { return parse_format_tags(tags, ctx); }
  template <class Context> auto format(const TelsemAtlas& x, Context& ctx) const {
    return tags.format(ctx, "TelsemAtlas(name="sv, x.name, ", month="sv, x.month, ", cells="sv, x.ndat, ")"sv);
  }
};

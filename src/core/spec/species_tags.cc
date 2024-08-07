#include "species_tags.h"

#include <debug.h>
#include <fast_float/fast_float.h>
#include <isotopologues.h>
#include <nonstd.h>
#include <partfun.h>

#include <algorithm>
#include <cfloat>
#include <iomanip>
#include <iterator>
#include <string_view>
#include <system_error>

#include "enums.h"
#include "species.h"

namespace Species {
namespace detail {
constexpr void trim(std::string_view& text) {
  while (text.size() and nonstd::isspace(text.front())) text.remove_prefix(1);
  while (text.size() and nonstd::isspace(text.back())) text.remove_suffix(1);
}

constexpr std::string_view next(std::string_view& text) {
  std::string_view next = text.substr(
      0, text.find_first_of('-', text.size() > 0 and text.front() == '-'));
  text.remove_prefix(std::min(text.size(), next.size() + 1));
  trim(next);
  trim(text);
  return next;
}

constexpr std::string_view next_tag(std::string_view& text) {
  std::string_view next = text.substr(0, text.find_first_of(','));
  text.remove_prefix(std::min(text.size(), next.size() + 1));
  trim(next);
  trim(text);
  return next;
}

constexpr SpeciesEnum spec(std::string_view part,
                           std::string_view orig [[maybe_unused]]) {
  return to<SpeciesEnum>(part);
}

constexpr Index isot(SpeciesEnum species,
                     std::string_view isot,
                     std::string_view orig [[maybe_unused]]) {
  Index spec_ind = find_species_index(species, isot);
  ARTS_USER_ERROR_IF(spec_ind < 0,
                     "Bad isotopologue: \"",
                     isot,
                     "\"\nValid options are:\n",
                     species,
                     "\nThe original tag reads \"",
                     orig,
                     '"')
  return spec_ind;
}

constexpr Numeric freq(std::string_view part,
                       std::string_view orig [[maybe_unused]]) {
  Numeric f;
  if (part.size() == 1 and part.front() == '*')
    f = -1;
  else {
    auto err =
        fast_float::from_chars(part.data(), part.data() + part.size(), f);
    ARTS_USER_ERROR_IF(err.ec == std::errc::invalid_argument,
                       "Invalid argument: \"",
                       part,
                       "\"\nThe original tag reads \"",
                       orig,
                       '"');
    ARTS_USER_ERROR_IF(err.ec == std::errc::result_out_of_range,
                       "Out-of-range: \"",
                       part,
                       "\"\nThe original tag reads \"",
                       orig,
                       '"');
  }
  return f;
}

constexpr void check(std::string_view text, std::string_view orig) {
  ARTS_USER_ERROR_IF(
      text.size(),
      "Parsing error.  The text \"",
      text,
      "\" remains to be parsed at end of a complete species tag\n"
      "The original tag reads \"",
      orig,
      '"')
}
}  // namespace detail

SpeciesTag parse_tag(std::string_view text) {
  using namespace detail;
  trim(text);

  const std::string_view orig = text;
  SpeciesTag tag;

  // The first part is a species --- we do not know what to do with it yet
  const SpeciesEnum species = spec(next(text), orig);

  // If there is no text remaining after the previous next(), then we are a
  // wild-tag species. Otherwise we have to process the tag a bit more
  if (text.size() == 0) {
    tag.spec_ind = isot(species, Joker, orig);
    tag.type     = SpeciesTagType::Plain;
    check(text, orig);
  } else {
    if (const std::string_view tag_key = next(text); tag_key == "CIA") {
      tag.spec_ind        = isot(species, Joker, orig);
      tag.cia_2nd_species = spec(next(text), orig);
      tag.type            = SpeciesTagType::Cia;
      check(text, orig);
    } else if (tag_key == "XFIT") {
      tag.spec_ind = isot(species, Joker, orig);
      tag.type     = SpeciesTagType::XsecFit;
      check(text, orig);
    } else {
      tag.spec_ind = isot(species, tag_key, orig);
      tag.type     = is_predefined_model(Isotopologues[tag.spec_ind])
                         ? SpeciesTagType::Predefined
                         : SpeciesTagType::Plain;
    }
  }

  if (text.size()) {
    check(text, orig);
  }

  return tag;
}

Array<Tag> parse_tags(std::string_view text) {
  using namespace detail;
  trim(text);

  Array<Tag> tags;
  while (text.size()) {
    tags.emplace_back(parse_tag(next_tag(text)));
  }

  return tags;
}

Numeric Tag::Q(Numeric T) const {
  return PartitionFunctions::Q(T, Isotopologue());
}

Numeric Tag::dQdT(Numeric T) const {
  return PartitionFunctions::dQdT(T, Isotopologue());
}

Tag::Tag(std::string_view text) : Tag(parse_tag(text)) {}

std::ostream& operator<<(std::ostream& os, const Tag& ot) {
  os << ot.Isotopologue().FullName();

  // Is this a CIA tag?
  if (ot.type == SpeciesTagType::Cia) {
    os << "-CIA-" << toString<1>(ot.cia_2nd_species);
  }

  else if (ot.type == SpeciesTagType::XsecFit) {
    os << "-XFIT";
  }

  return os;
}

/** Return the full name of the tag.
 * 
 * \return The tag name as a string.
 */
String Tag::Name() const { return var_string(*this); }
}  // namespace Species

ArrayOfSpeciesTag::ArrayOfSpeciesTag(std::string_view text)
    : ArrayOfSpeciesTag(Species::parse_tags(text)){ARTS_USER_ERROR_IF(
          size() and std::any_of(begin(),
                                 end(),
                                 [front_species =
                                      front().Spec()](const SpeciesTag& tag) {
                                   return tag.Spec() not_eq front_species;
                                 }),
          "All species in a group of tags must be the same\n"
          "Your list of tags have been parsed as: ",
          *this,
          "\nThe original tags-list read \"",
          text,
          '"')}

      Index find_next_species(const ArrayOfArrayOfSpeciesTag& specs,
                              SpeciesEnum spec,
                              Index i) noexcept {
  const Index n = static_cast<Index>(specs.size());
  for (; i < n; i++)
    if (specs[i].Species() == spec) return i;
  return -1;
}

Index find_first_species(const ArrayOfArrayOfSpeciesTag& specs,
                         SpeciesEnum spec) noexcept {
  return find_next_species(specs, spec, 0);
}

std::pair<Index, Index> find_first_species_tag(
    const ArrayOfArrayOfSpeciesTag& specs, const SpeciesTag& tag) noexcept {
  for (Size i = 0; i < specs.size(); i++) {
    if (auto ptr = std::find(specs[i].cbegin(), specs[i].cend(), tag);
        ptr not_eq specs[i].cend())
      return {i, std::distance(specs[i].cbegin(), ptr)};
  }
  return {-1, -1};
}

std::pair<Index, Index> find_first_isotologue(
    const ArrayOfArrayOfSpeciesTag& specs,
    const SpeciesIsotope& isot) noexcept {
  for (Size i = 0; i < specs.size(); i++) {
    if (auto ptr =
            std::find_if(specs[i].cbegin(),
                         specs[i].cend(),
                         [&](auto& tag) { return tag.Isotopologue() == isot; });
        ptr not_eq specs[i].cend())
      return {i, std::distance(specs[i].cbegin(), ptr)};
  }
  return {-1, -1};
}

String ArrayOfSpeciesTag::Name() const {
  String out = "";
  bool first = true;
  for (auto& x : *this) {
    if (not first) out += ", ";
    out   += x.Name();
    first  = false;
  }
  return out;
}

std::set<SpeciesEnum> lbl_species(
    const ArrayOfArrayOfSpeciesTag& abs_species) noexcept {
  std::set<SpeciesEnum> unique_species;
  for (auto& specs : abs_species) {
    if (specs.RequireLines()) unique_species.insert(specs.front().Spec());
  }
  return unique_species;
}

Numeric Species::first_vmr(const ArrayOfArrayOfSpeciesTag& abs_species,
                           const Vector& rtp_vmr,
                           const SpeciesEnum spec) ARTS_NOEXCEPT {
  ARTS_ASSERT(abs_species.size() == static_cast<Size>(rtp_vmr.size()))

  auto pos =
      std::find_if(abs_species.begin(),
                   abs_species.end(),
                   [spec](const ArrayOfSpeciesTag& tag_group) {
                     return tag_group.size() and tag_group.Species() == spec;
                   });
  return pos == abs_species.end()
             ? 0
             : rtp_vmr[std::distance(abs_species.begin(), pos)];
}

SpeciesTagTypeStatus::SpeciesTagTypeStatus(
    const ArrayOfArrayOfSpeciesTag& abs_species) {
  for (auto& species_list : abs_species) {
    for (auto& tag : species_list) {
      switch (tag.type) {
        case SpeciesTagType::Plain:
          Plain = true;
          break;
        case SpeciesTagType::Predefined:
          Predefined = true;
          break;
        case SpeciesTagType::Cia:
          Cia = true;
          break;
        case SpeciesTagType::XsecFit:
          XsecFit = true;
          break;
      }
    }
  }
}

std::ostream& operator<<(std::ostream& os, SpeciesTagTypeStatus val) {
  SpeciesTagType x{SpeciesTagType::Plain};
  os << "Species tag types:\n";
  switch (x) {
    case SpeciesTagType::Plain:
      os << "    Plain:            " << val.Plain << '\n';
      [[fallthrough]];
    case SpeciesTagType::Predefined:
      os << "    Predefined: " << val.Predefined << '\n';
      [[fallthrough]];
    case SpeciesTagType::Cia:
      os << "    Cia:              " << val.Cia << '\n';
      [[fallthrough]];
    case SpeciesTagType::XsecFit:
      os << "    XsecFit:       " << val.XsecFit << '\n';
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfSpeciesTag& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

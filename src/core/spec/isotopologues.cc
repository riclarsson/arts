#include "isotopologues.h"

#include <debug.h>
#include <xml_io_base.h>

#include <cassert>
#include <compare>
#include <sstream>

namespace Species {
ArrayOfSpeciesIsotope isotopologues(SpeciesEnum spec) {
#define deal_with_spec(SPEC)                                      \
  case SpeciesEnum::SPEC: {                                       \
    static constexpr auto v = isotopologues<SpeciesEnum::SPEC>(); \
    return {v.begin(), v.end()};                                  \
  } break

  switch (spec) {
    case SpeciesEnum::Bath:
      break;
      deal_with_spec(Water);
      deal_with_spec(CarbonDioxide);
      deal_with_spec(Ozone);
      deal_with_spec(NitrogenOxide);
      deal_with_spec(CarbonMonoxide);
      deal_with_spec(Methane);
      deal_with_spec(Oxygen);
      deal_with_spec(NitricOxide);
      deal_with_spec(SulfurDioxide);
      deal_with_spec(NitrogenDioxide);
      deal_with_spec(Ammonia);
      deal_with_spec(NitricAcid);
      deal_with_spec(Hydroxyl);
      deal_with_spec(HydrogenFluoride);
      deal_with_spec(HydrogenChloride);
      deal_with_spec(HydrogenBromide);
      deal_with_spec(HydrogenIodide);
      deal_with_spec(ChlorineMonoxide);
      deal_with_spec(CarbonylSulfide);
      deal_with_spec(Formaldehyde);
      deal_with_spec(HeavyFormaldehyde);
      deal_with_spec(VeryHeavyFormaldehyde);
      deal_with_spec(HypochlorousAcid);
      deal_with_spec(Nitrogen);
      deal_with_spec(HydrogenCyanide);
      deal_with_spec(Chloromethane);
      deal_with_spec(HydrogenPeroxide);
      deal_with_spec(Acetylene);
      deal_with_spec(Ethane);
      deal_with_spec(Phosphine);
      deal_with_spec(CarbonylFluoride);
      deal_with_spec(SulfurHexafluoride);
      deal_with_spec(HydrogenSulfide);
      deal_with_spec(FormicAcid);
      deal_with_spec(LeftHeavyFormicAcid);
      deal_with_spec(RightHeavyFormicAcid);
      deal_with_spec(Hydroperoxyl);
      deal_with_spec(OxygenAtom);
      deal_with_spec(ChlorineNitrate);
      deal_with_spec(NitricOxideCation);
      deal_with_spec(HypobromousAcid);
      deal_with_spec(Ethylene);
      deal_with_spec(Methanol);
      deal_with_spec(Bromomethane);
      deal_with_spec(Acetonitrile);
      deal_with_spec(HeavyAcetonitrile);
      deal_with_spec(CarbonTetrafluoride);
      deal_with_spec(Diacetylene);
      deal_with_spec(Cyanoacetylene);
      deal_with_spec(Hydrogen);
      deal_with_spec(CarbonMonosulfide);
      deal_with_spec(SulfurTrioxide);
      deal_with_spec(Cyanogen);
      deal_with_spec(Phosgene);
      deal_with_spec(SulfurMonoxide);
      deal_with_spec(CarbonDisulfide);
      deal_with_spec(Methyl);
      deal_with_spec(Cyclopropene);
      deal_with_spec(SulfuricAcid);
      deal_with_spec(HydrogenIsocyanide);
      deal_with_spec(BromineMonoxide);
      deal_with_spec(ChlorineDioxide);
      deal_with_spec(Propane);
      deal_with_spec(Helium);
      deal_with_spec(ChlorineMonoxideDimer);
      deal_with_spec(HydrogenAtom);
      deal_with_spec(Argon);
      deal_with_spec(Hexafluoroethane);
      deal_with_spec(Perfluoropropane);
      deal_with_spec(Perfluorobutane);
      deal_with_spec(Perfluoropentane);
      deal_with_spec(Perfluorohexane);
      deal_with_spec(Perfluorooctane);
      deal_with_spec(Perfluorocyclobutane);
      deal_with_spec(CarbonTetrachloride);
      deal_with_spec(CFC11);
      deal_with_spec(CFC113);
      deal_with_spec(CFC114);
      deal_with_spec(CFC115);
      deal_with_spec(CFC12);
      deal_with_spec(Dichloromethane);
      deal_with_spec(Trichloroethane);
      deal_with_spec(Trichloromethane);
      deal_with_spec(Bromochlorodifluoromethane);
      deal_with_spec(Bromotrifluoromethane);
      deal_with_spec(Dibromotetrafluoroethane);
      deal_with_spec(HCFC141b);
      deal_with_spec(HCFC142b);
      deal_with_spec(HCFC22);
      deal_with_spec(HFC125);
      deal_with_spec(HFC134a);
      deal_with_spec(HFC143a);
      deal_with_spec(HFC152a);
      deal_with_spec(HFC227ea);
      deal_with_spec(HFC23);
      deal_with_spec(HFC236fa);
      deal_with_spec(HFC245fa);
      deal_with_spec(HFC32);
      deal_with_spec(HFC365mfc);
      deal_with_spec(NitrogenTrifluoride);
      deal_with_spec(TrihydrogenCation);
      deal_with_spec(SulfurDimer);
      deal_with_spec(CarbonylChlorofluoride);
      deal_with_spec(NitrousAcid);
      deal_with_spec(NitrylChloride);
      deal_with_spec(SulfurylFluoride);
      deal_with_spec(HFC4310mee);
      deal_with_spec(Germane);
      deal_with_spec(Iodomethane);
      deal_with_spec(Fluoromethane);
      deal_with_spec(Arsine);
      deal_with_spec(Benzene);
      deal_with_spec(liquidcloud);
      deal_with_spec(icecloud);
      deal_with_spec(rain);
      deal_with_spec(free_electrons);
      deal_with_spec(particles);
      deal_with_spec(unused);
  }

#undef deal_with_spec

  ARTS_USER_ERROR("Cannot understand: {}", spec)
}

String isotopologues_names(SpeciesEnum spec) {
  auto x = isotopologues(spec);
  std::ostringstream os;
  for (auto& s : x) os << s.FullName() << '\n';
  return os.str();
}

String update_isot_name(const String& old_name) {
  if (old_name == "CH3CN-211224") return "CH2DCN-224";
  if (old_name == "CH3CN-211124") return "CH3CN-2124";
  if (old_name == "CH3CN-211125") return "CH3CN-2125";
  if (old_name == "CH3CN-211134") return "CH3CN-2134";
  if (old_name == "CH3CN-311124") return "CH3CN-3124";
  if (old_name == "CO2-728") return "CO2-827";
  if (old_name == "HCOOH-2261") return "DCOOH-266";
  if (old_name == "HCOOH-1262") return "HCOOD-266";
  if (old_name == "HCOOH-1261") return "HCOOH-126";
  if (old_name == "HCOOH-1361") return "HCOOH-136";
  if (old_name == "H2CO-1126") return "H2CO-126";
  if (old_name == "H2CO-1128") return "H2CO-128";
  if (old_name == "H2CO-1136") return "H2CO-136";
  if (old_name == "H2CO-1226") return "HDCO-26";
  if (old_name == "H2CO-2226") return "D2CO-26";
  return old_name;
}

const Isotope& Isotope::from_name(const std::string_view name) {
  return select(name);
}

String Isotope::FullName() const {
  return is_joker() ? String{toString<1>(spec)}
                    : std::format("{}-{}", toString<1>(spec), isotname);
}

std::ostream& operator<<(std::ostream& os, const Isotope& ir) {
  return os << ir.FullName();
}

std::ostream& operator<<(std::ostream& os, const std::vector<Isotope>& isots) {
  for (const auto& i : isots) os << i << ' ';
  return os;
}

Numeric IsotopologueRatios::operator[](const Index spec_ind) const {
  assert(spec_ind < maxsize and spec_ind >= 0);
  return data[spec_ind];
}

Numeric IsotopologueRatios::operator[](const Isotope& ir) const {
  const Index spec_ind = find_species_index(ir);
  ARTS_USER_ERROR_IF(
      spec_ind >= maxsize or spec_ind < 0, "Invalid species {}", ir.FullName())
  return data[spec_ind];
}

std::vector<std::string> IsotopologueRatios::valueless_isotopes() const {
  std::vector<std::string> names;

  for (Index i = 0; i < maxsize; i++) {
    if (not Isotopologues[i].is_predefined() and
        not Isotopologues[i].is_joker() and nonstd::isnan(data[i])) {
      names.push_back(Isotopologues[i].FullName());
    }
  }

  return names;
}

namespace {
consteval IsotopologueRatios from_builtin() {
  IsotopologueRatios isotopologue_ratios{};

  stdr::copy(Isotopologues | std::views::transform(&Isotope::builtin_ratio),
             isotopologue_ratios.data.begin());

  return isotopologue_ratios;
}

consteval bool all_values() {
  for (const auto& iso : Isotopologues) {
    if (not(iso.is_predefined() or iso.is_joker())) {
      if (nonstd::isnan(iso.builtin_ratio) or nonstd::isnan(iso.mass) or
          iso.gi == 0) {
        return false;
      }
    }
  }
  return true;
}

static_assert(all_values(),
              "Some isotopologues do not have all values defined.  "
              "Please check the source files.  "
              "You must define the missing values.");
}  // namespace

const IsotopologueRatios& isotopologue_ratiosInitFromBuiltin() {
  static constexpr IsotopologueRatios isotopologue_ratios = from_builtin();
  return isotopologue_ratios;
}

std::strong_ordering Isotope::operator<=>(const Isotope& other) const {
  if (std::strong_ordering test = spec <=> other.spec;
      test != std::strong_ordering::equal)
    return test;
  return isotname <=> other.isotname;
}
static_assert(std::strong_ordering::equal == std::strong_ordering::equivalent);

bool Isotope::operator==(const Isotope& other) const {
  return this->operator<=>(other) == std::strong_ordering::equal;
}

bool Isotope::operator!=(const Isotope& other) const {
  return not this->operator==(other);
}
}  // namespace Species

void xml_io_stream<SpeciesIsotope>::write(std::ostream& os,
                                          const SpeciesIsotope& x,
                                          bofstream*,
                                          std::string_view name) {
  XMLTag tag(type_name, "name", name, "isot", x.FullName());
  tag.write_to_stream(os);
  tag.write_to_end_stream(os);
}

void xml_io_stream<SpeciesIsotope>::read(std::istream& is,
                                         SpeciesIsotope& x,
                                         bifstream*) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  String v;
  tag.get_attribute_value("isot", v);
  x = SpeciesIsotope::from_name(v);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

#include "nonlte_atm.h"

#include <partfun.h>

namespace nonlte {
Numeric level_density(Numeric T,
                      Numeric g,
                      Numeric E,
                      const SpeciesIsotope& spec) {
  using Constant::k;
  return g * std::exp(-E / (k * T)) / PartitionFunctions::Q(T, spec);
}

data_t level_density(const AbsorptionBands&, const AtmFunctionalData&) {
  throw std::runtime_error("Not implemented for functional temperature field");

  return {};
}

data_t level_density(const AbsorptionBands&, const Numeric&) {
  throw std::runtime_error("Not implemented for Numeric temperature field");

  return {};
}

data_t level_density(const AbsorptionBands& bands, const GriddedField3& ts) {
  data_t out;

  for (auto& [qid, band] : bands) {
    ARTS_USER_ERROR_IF(
        band.size() != 1,
        "Only one line per band is supported:\nQID: {}\nBand: {}",
        qid,
        band);

    const auto& line = band.lines.front();

    if (const auto lvl = qid.LowerLevel(); not out.contains(lvl)) {
      const auto level = [g = line.gl, E = line.el(), &qid](Numeric T) {
        return level_density(T, g, E, qid.Isotopologue());
      };
      GriddedField3 data = ts;
      assert(false);
      // data.data.unary_transform(level);
      out[lvl] = std::move(data);
    }

    if (const auto lvl = qid.UpperLevel(); not out.contains(lvl)) {
      const auto level = [g = line.gu, E = line.eu(), &qid](Numeric T) {
        return level_density(T, g, E, qid.Isotopologue());
      };
      GriddedField3 data = ts;
      assert(false);
      // data.data.unary_transform(level);
      out[lvl] = std::move(data);
    }
  }

  return out;
}

data_t level_density(const AtmField& atm, const AbsorptionBands& bands) {
  return std::visit(
      [&bands](const auto& temperatures) {
        return level_density(bands, temperatures);
      },
      atm.other().at(AtmKey::t).data);
}
}  // namespace nonlte

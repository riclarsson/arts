#pragma once

#include <memory>

#include <predefined_absorption_models.h>
#include <species_tags.h>

namespace predef {
struct full {
  std::vector<Species::IsotopeRecord> tags;
  Numeric P;
  Numeric T;
  Absorption::PredefinedModel::VMRS vmrs;
  std::shared_ptr<PredefinedModelData> predefined_model_data;

  full() = default;

  full(Numeric p,
       Numeric t,
       const Vector& allvmrs,
       const ArrayOfArrayOfSpeciesTag& allspecs,
       const std::shared_ptr<PredefinedModelData>& data);

  [[nodiscard]] Complex at(Numeric f) const;
  void at(ComplexVector& abs, const Vector& fs) const;
  [[nodiscard]] ComplexVector at(const Vector& fs) const;
};
}  // namespace predef

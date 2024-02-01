#include "fwd_predef.h"

#include <algorithm>

#include "atm.h"
#include "jacobian.h"
#include "rtepack.h"

namespace fwd::predef {
void full::adapt() { vmrs = Absorption::PredefinedModel::VMRS(*atm); }

full::full(std::shared_ptr<AtmPoint> atm_,
           std::shared_ptr<PredefinedModelData> data_)
    : atm(std::move(atm_)), data(std::move(data_)) {
  adapt();
}

Complex full::operator()(const Numeric frequency) const {
  PropmatVector propmat_clearsky(1);
  PropmatMatrix dpropmat_clearsky_dx;
  JacobianTargets jacobian_targets;
  Vector f_grid{frequency};

  for (auto& [tag, mod] : *data) {
    Absorption::PredefinedModel::compute(propmat_clearsky,
                                         dpropmat_clearsky_dx,
                                         tag,
                                         f_grid,
                                         atm->pressure,
                                         atm->temperature,
                                         vmrs,
                                         jacobian_targets,
                                         mod);
  }

  return propmat_clearsky[0].A();
}

void full::set_model(std::shared_ptr<PredefinedModelData> data_) {
  data = std::move(data_);
  adapt();
}

void full::set_atm(std::shared_ptr<AtmPoint> atm_) {
  atm = std::move(atm_);
  adapt();
}
}  // namespace fwd::predef

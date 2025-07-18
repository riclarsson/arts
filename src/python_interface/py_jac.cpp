#include <python_interface.h>

#include "atm.h"
#include "hpy_arts.h"
#include "jacobian.h"
#include "lbl_data.h"
#include "nanobind/stl/bind_vector.h"
#include "nanobind/stl/function.h"
#include "nanobind/stl/variant.h"
#include "surf.h"

NB_MAKE_OPAQUE(std::vector<Jacobian::AtmTarget>)
NB_MAKE_OPAQUE(std::vector<Jacobian::SurfaceTarget>)
NB_MAKE_OPAQUE(std::vector<Jacobian::LineTarget>)
NB_MAKE_OPAQUE(std::vector<Jacobian::SensorTarget>)
NB_MAKE_OPAQUE(std::vector<Jacobian::ErrorTarget>)

namespace Python {
void py_jac(py::module_& m) try {
  py::class_<ErrorKey> errkey(m, "ErrorKey");
  errkey.doc() = "Error key";
  errkey.def_rw("y_start",
                &ErrorKey::y_start,
                "Start index of target in measurement vector");
  errkey.def_rw(
      "y_size", &ErrorKey::y_size, "Size of target in measurement vector");

  py::class_<Jacobian::AtmTarget> atm(m, "JacobianAtmTarget");
  atm.doc() = "Atmospheric target";
  atm.def_ro("type", &Jacobian::AtmTarget::type, "Type of target");
  atm.def_ro("d", &Jacobian::AtmTarget::d, "Perturbation magnitude");
  atm.def_ro("target_pos", &Jacobian::AtmTarget::target_pos, "Target position");
  atm.def_ro("x_start",
             &Jacobian::AtmTarget::x_start,
             "Start index of target in state vector");
  atm.def_ro(
      "x_size", &Jacobian::AtmTarget::x_size, "Size of target in state vector");
  str_interface(atm);

  py::class_<Jacobian::SurfaceTarget> surf(m, "JacobianSurfaceTarget");
  surf.doc() = "Surface target";
  surf.def_ro("type", &Jacobian::SurfaceTarget::type, "Type of target");
  surf.def_ro("d", &Jacobian::SurfaceTarget::d, "Perturbation magnitude");
  surf.def_ro(
      "target_pos", &Jacobian::SurfaceTarget::target_pos, "Target position");
  surf.def_ro("x_start",
              &Jacobian::SurfaceTarget::x_start,
              "Start index of target in state vector");
  surf.def_ro("x_size",
              &Jacobian::SurfaceTarget::x_size,
              "Size of target in state vector");
  str_interface(surf);

  py::class_<Jacobian::LineTarget> line(m, "JacobianLineTarget");
  line.doc() = "Line target";
  line.def_ro("type", &Jacobian::LineTarget::type, "Type of target");
  line.def_ro("d", &Jacobian::LineTarget::d, "Perturbation magnitude");
  line.def_ro(
      "target_pos", &Jacobian::LineTarget::target_pos, "Target position");
  line.def_ro("x_start",
              &Jacobian::LineTarget::x_start,
              "Start index of target in state vector");
  line.def_ro("x_size",
              &Jacobian::LineTarget::x_size,
              "Size of target in state vector");
  str_interface(line);

  py::class_<Jacobian::SensorTarget> sensor(m, "JacobianSensorTarget");
  sensor.doc() = "Sensor target";
  sensor.def_ro("type", &Jacobian::SensorTarget::type, "Type of target");
  sensor.def_ro("d", &Jacobian::SensorTarget::d, "Perturbation magnitude");
  sensor.def_ro(
      "target_pos", &Jacobian::SensorTarget::target_pos, "Target position");
  sensor.def_ro("x_start",
                &Jacobian::SensorTarget::x_start,
                "Start index of target in state vector");
  sensor.def_ro("x_size",
                &Jacobian::SensorTarget::x_size,
                "Size of target in state vector");
  str_interface(sensor);

  py::class_<Jacobian::ErrorTarget> error(m, "JacobianErrorTarget");
  error.doc() = "Error target";
  error.def_ro("type", &Jacobian::ErrorTarget::type, "Type of target");
  error.def_ro(
      "target_pos", &Jacobian::ErrorTarget::target_pos, "Target position");
  error.def_ro("x_start",
               &Jacobian::ErrorTarget::x_start,
               "Start index of target in state vector");
  error.def_ro("x_size",
               &Jacobian::ErrorTarget::x_size,
               "Size of target in state vector");
  str_interface(error);

  py::bind_vector<std::vector<Jacobian::AtmTarget>>(m, "ArrayOfAtmTargets")
      .doc() = "List of atmospheric targets";
  py::bind_vector<std::vector<Jacobian::SurfaceTarget>>(m,
                                                        "ArrayOfSurfaceTarget")
      .doc() = "List of surface targets";
  py::bind_vector<std::vector<Jacobian::LineTarget>>(m, "ArrayOfLineTarget")
      .doc() = "List of line targets";
  py::bind_vector<std::vector<Jacobian::SensorTarget>>(m, "ArrayOfSensorTarget")
      .doc() = "List of sensor targets";
  py::bind_vector<std::vector<Jacobian::ErrorTarget>>(m, "ArrayOfErrorTarget")
      .doc() = "List of error targets";

  py::class_<JacobianTargets> jacs(m, "JacobianTargets");
  generic_interface(jacs);
  jacs.def_prop_rw(
      "atm",
      [](JacobianTargets& j) { return j.atm(); },
      [](JacobianTargets& j, std::vector<Jacobian::AtmTarget> t) {
        j.atm() = std::move(t);
      },
      "List of atmospheric targets");
  jacs.def_prop_rw(
      "surf",
      [](JacobianTargets& j) { return j.surf(); },
      [](JacobianTargets& j, std::vector<Jacobian::SurfaceTarget> t) {
        j.surf() = std::move(t);
      },
      "List of surface targets");
  jacs.def_prop_rw(
      "line",
      [](JacobianTargets& j) { return j.line(); },
      [](JacobianTargets& j, std::vector<Jacobian::LineTarget> t) {
        j.line() = std::move(t);
      },
      "List of line targets");
  jacs.def_prop_rw(
      "sensor",
      [](JacobianTargets& j) { return j.sensor(); },
      [](JacobianTargets& j, std::vector<Jacobian::SensorTarget> t) {
        j.sensor() = std::move(t);
      },
      "List of sensor targets");
  jacs.def_prop_rw(
      "error",
      [](JacobianTargets& j) { return j.error(); },
      [](JacobianTargets& j, std::vector<Jacobian::ErrorTarget> t) {
        j.error() = std::move(t);
      },
      "List of error targets");
  jacs.def(
      "x_size",
      [](const JacobianTargets& j){return j.x_size();},
      "The number of elements required in the *model_state_vector* to represent the Jacobian.");
  jacs.def(
      "target_count",
      [](const JacobianTargets& j){return j.target_count();},
      "The number of targets added to the Jacobian.");
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize jac\n{}", e.what()));
}
}  // namespace Python
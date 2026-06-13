#include <lupnt/lupnt.h>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitMathUtils(py::module& m) {
  m.def("wrap_to_pi", py::overload_cast<Real>(&WrapToPi), py::arg("angle"));
  m.def("wrap_to_two_pi", py::overload_cast<Real>(&WrapToTwoPi), py::arg("angle"));

  m.def("wrap_to_pi", py::overload_cast<const VecX&>(&WrapToPi), py::arg("angle"));
  m.def("wrap_to_two_pi", py::overload_cast<const VecX&>(&WrapToTwoPi), py::arg("angle"));

  m.def("deg_min_sec_to_degrees", py::overload_cast<const Vec3&>(&DegMinSecToDeg), py::arg("hms"));
  m.def("degrees_to_deg_min_sec", py::overload_cast<Real>(&DegToDegMinSec), py::arg("deg"));

  // Decibel-Decimal conversion
  m.def("decibel_to_decimal", py::overload_cast<Real>(&DecibelToDecimal), py::arg("x"));
  m.def("decimal_to_decibel", py::overload_cast<Real>(&DecimalToDecibel), py::arg("x"));
  m.def("decimal_to_decibel", py::overload_cast<const ArrX&>(&DecimalToDecibel), py::arg("x"));
  m.def("decibel_to_decimal", py::overload_cast<const ArrX&>(&DecibelToDecimal), py::arg("x"));

  // Rotation
  m.def("rot_x", &RotX, py::arg("angle"));
  m.def("rot_y", &RotY, py::arg("angle"));
  m.def("rot_z", &RotZ, py::arg("angle"));
  m.def("skew", &Skew, py::arg("x"));

  m.def("rotation_angle", &RotationAngle, py::arg("R"));
}

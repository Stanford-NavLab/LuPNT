#include <lupnt/lupnt.h>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitMathUtils(py::module& m) {
  m.def("wrap_to_pi", py::overload_cast<Real>(&WrapToPi), py::arg("angle"),
        "Wrap an angle [rad] to (-pi, pi].");
  m.def("wrap_to_two_pi", py::overload_cast<Real>(&WrapToTwoPi), py::arg("angle"),
        "Wrap an angle [rad] to [0, 2*pi).");

  m.def("wrap_to_pi", py::overload_cast<const VecX&>(&WrapToPi), py::arg("angle"),
        "Element-wise wrap of angles [rad] to (-pi, pi].");
  m.def("wrap_to_two_pi", py::overload_cast<const VecX&>(&WrapToTwoPi), py::arg("angle"),
        "Element-wise wrap of angles [rad] to [0, 2*pi).");

  m.def("deg_min_sec_to_degrees", py::overload_cast<const Vec3&>(&DegMinSecToDeg), py::arg("hms"),
        "Convert (degrees, arcminutes, arcseconds) to decimal degrees.");
  m.def("degrees_to_deg_min_sec", py::overload_cast<Real>(&DegToDegMinSec), py::arg("deg"),
        "Convert decimal degrees to (degrees, arcminutes, arcseconds).");

  // Decibel-Decimal conversion
  m.def("decibel_to_decimal", py::overload_cast<Real>(&DecibelToDecimal), py::arg("x"),
        "Convert a value in decibels [dB] to a linear ratio: 10^(x/10).");
  m.def("decimal_to_decibel", py::overload_cast<Real>(&DecimalToDecibel), py::arg("x"),
        "Convert a linear ratio to decibels [dB]: 10*log10(x).");
  m.def("decimal_to_decibel", py::overload_cast<const ArrX&>(&DecimalToDecibel), py::arg("x"),
        "Element-wise conversion of linear ratios to decibels [dB].");
  m.def("decibel_to_decimal", py::overload_cast<const ArrX&>(&DecibelToDecimal), py::arg("x"),
        "Element-wise conversion of decibels [dB] to linear ratios.");

  // Rotation
  m.def("rot_x", &RotX, py::arg("angle"),
        "3x3 passive (frame) rotation matrix about the x-axis by angle [rad].");
  m.def("rot_y", &RotY, py::arg("angle"),
        "3x3 passive (frame) rotation matrix about the y-axis by angle [rad].");
  m.def("rot_z", &RotZ, py::arg("angle"),
        "3x3 passive (frame) rotation matrix about the z-axis by angle [rad].");
  m.def("skew", &Skew, py::arg("x"),
        "3x3 skew-symmetric cross-product matrix [x]_x of a 3-vector.");

  m.def("rotation_angle", &RotationAngle, py::arg("R"),
        "Rotation angle [rad] of a 3x3 rotation matrix.");
}

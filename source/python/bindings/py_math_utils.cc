#include <lupnt/numerics/math_utils.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lupnt;

void init_math_utils(py::module &m) {
  // m.def("Wrap2Pi", &Wrap2Pi, "Wrap angle to [-pi, pi]");
  m.def(
      "Wrap2Pi", [](double angle) -> double { return Wrap2Pi(angle).val(); },
      "Wrap angle to [-pi, pi]");
  m.def(
      "Wrap2TwoPi", [](double angle) -> double { return Wrap2TwoPi(angle).val(); },
      "Wrap angle to [0, 2pi]");
  m.def(
      "Wrap2Pi", [](const VecXd &angle) -> VecXd { return Wrap2Pi(angle).cast<double>(); },
      "Wrap angle to [-pi, pi]");
  m.def(
      "Wrap2TwoPi", [](const VecXd &angle) -> VecXd { return Wrap2TwoPi(angle).cast<double>(); },
      "Wrap angle to [0, 2pi]");
  m.def(
      "DegMinSec2Degrees", [](const Vec3d &hms) -> double { return DegMinSec2Degrees(hms).val(); },
      "Convert degrees, minutes, seconds to degrees", py::arg("hms"));
  m.def(
      "Degrees2DegMinSec",
      [](double deg) -> Vec3d { return Degrees2DegMinSec(deg).cast<double>(); },
      "Convert degrees to degrees, minutes, seconds", py::arg("deg"));
  m.def(
      "Decimal2Decibel", [](double x) -> double { return Decimal2Decibel(x).val(); },
      "Convert decimal to dB", py::arg("x"));
  m.def(
      "Decibel2Decimal", [](double x) -> double { return Decibel2Decimal(x).val(); },
      "Convert dB to decimal", py::arg("x"));
  m.def(
      "Decibel2Decimal", [](VecXd x) -> VecXd { return Decibel2Decimal(x).cast<double>(); },
      py::arg("x"));
  m.def(
      "Decimal2Decibel", [](VecXd x) -> VecXd { return Decimal2Decibel(x).cast<double>(); },
      py::arg("x"));
}

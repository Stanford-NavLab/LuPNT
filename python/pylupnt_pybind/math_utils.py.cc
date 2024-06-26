#include <lupnt/numerics/math_utils.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lupnt;

void init_math_utils(py::module &m) {
  // m.def("wrapToPi", &wrapToPi, "Wrap angle to [-pi, pi]");
  m.def(
      "wrapToPi", [](double angle) -> double { return wrapToPi(angle).val(); },
      "Wrap angle to [-pi, pi]");
  m.def(
      "wrapTo2Pi",
      [](double angle) -> double { return wrapTo2Pi(angle).val(); },
      "Wrap angle to [0, 2pi]");
  m.def(
      "wrapToPi",
      [](const VecXd &angle) -> VecXd {
        return wrapToPi(angle).cast<double>();
      },
      "Wrap angle to [-pi, pi]");
  m.def(
      "wrapTo2Pi",
      [](const VecXd &angle) -> VecXd {
        return wrapTo2Pi(angle).cast<double>();
      },
      "Wrap angle to [0, 2pi]");
  m.def(
      "dms2degrees",
      [](const Vec3d &hms) -> double { return dms2degrees(hms).val(); },
      "Convert degrees, minutes, seconds to degrees", py::arg("hms"));
  m.def(
      "degrees2dms",
      [](double deg) -> Vec3d { return degrees2dms(deg).cast<double>(); },
      "Convert degrees to degrees, minutes, seconds", py::arg("deg"));
  m.def(
      "decimal2dB", [](double x) -> double { return decimal2dB(x).val(); },
      "Convert decimal to dB", py::arg("x"));
  m.def(
      "dB2decimal", [](double x) -> double { return dB2decimal(x).val(); },
      "Convert dB to decimal", py::arg("x"));
  m.def(
      "dB2decimal",
      [](VecXd x) -> VecXd { return dB2decimal(x).cast<double>(); },
      py::arg("x"));
  m.def(
      "decimal2dB",
      [](VecXd x) -> VecXd { return decimal2dB(x).cast<double>(); },
      py::arg("x"));
}

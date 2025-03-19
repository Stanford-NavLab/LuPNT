#include <lupnt/core/constants.h>
#include <lupnt/numerics/math_utils.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lupnt;

void init_math_utils(py::module &m) {
  // m.def("Wrap2Pi", &Wrap2Pi, "Wrap angle to [-pi, pi]");
  m.def(
      "wrap2pi", [](double angle) -> double { return Wrap2Pi(angle).val(); },
      "Wrap angle to [-pi, pi]");
  m.def(
      "wrap2two_pi", [](double angle) -> double { return Wrap2TwoPi(angle).val(); },
      "Wrap angle to [0, 2pi]");
  m.def(
      "wrap2pi", [](const VecXd &angle) -> VecXd { return Wrap2Pi(angle).cast<double>(); },
      "Wrap angle to [-pi, pi]");
  m.def(
      "wrap2two_pi", [](const VecXd &angle) -> VecXd { return Wrap2TwoPi(angle).cast<double>(); },
      "Wrap angle to [0, 2pi]");

  m.def(
      "deg_min_sec2degrees",
      [](const Vec3d &hms) -> double { return DegMinSec2Degrees(hms).val(); },
      "Convert degrees, minutes, seconds to degrees", py::arg("hms"));
  m.def(
      "degrees2deg_min_sec",
      [](double deg) -> Vec3d { return Degrees2DegMinSec(deg).cast<double>(); },
      "Convert degrees to degrees, minutes, seconds", py::arg("deg"));

  // Decibel-Decimal conversion
  m.def(
      "decibel2decimal", [](double x) -> double { return Decibel2Decimal(x).val(); },
      "Convert dB to decimal", py::arg("x"));
  m.def(
      "decimal2decibel", [](double x) -> double { return Decimal2Decibel(x).val(); },
      "Convert decimal to dB", py::arg("x"));
  m.def(
      "decibel2decimal",
      [](VecXd x) -> VecXd { return Decibel2Decimal(x.cast<Real>().eval()).cast<double>(); },
      py::arg("x"));
  m.def(
      "decimal2decibel",
      [](VecXd x) -> VecXd { return Decimal2Decibel(x.cast<Real>().eval()).cast<double>(); },
      py::arg("x"));
  m.def(
      "decibel2decimal",
      [](MatXd x) -> MatXd { return Decibel2Decimal(x.cast<Real>().eval()).cast<double>(); },
      py::arg("x"));
  m.def(
      "decimal2decibel",
      [](MatXd x) -> MatXd { return Decimal2Decibel(x.cast<Real>().eval()).cast<double>(); },
      py::arg("x"));

  // Rotation
  m.def("rot_x", &RotX<double>, "Active rotation matrix about the x-axis", py::arg("angle"));
  m.def("rot_y", &RotY<double>, "Active rotation matrix about the y-axis", py::arg("angle"));
  m.def("rot_z", &RotZ<double>, "Active rotation matrix about the z-axis", py::arg("angle"));
  m.def("skew", &Skew<double>, "Skew symmetric matrix from a vector", py::arg("x"));
  m.def(
      "rot2quat",
      [](const Mat3d &rot) -> Vec4d {
        Quatd quat(rot);
        quat.normalize();
        if (quat.w() < 0) quat.coeffs() *= -1;
        return quat.coeffs();
      },
      "Convert rotation matrix to quaternion [x, y, z, w]", py::arg("R"));
  m.def(
      "quat2rot",
      [](const Vec4d &q) -> Mat3d {
        Quatd quat(q(3), q(0), q(1), q(2));
        return quat.toRotationMatrix();
      },
      "Convert quaternion to rotation matrix", py::arg("q"));
}

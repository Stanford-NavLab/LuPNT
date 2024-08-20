// lupnt
#include <lupnt/core/constants.h>
#include <lupnt/dynamics/forces.h>
#include <lupnt/numerics/math_utils.h>

// pybind11
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>

#include "py_vectorized_macros.cc"

namespace py = pybind11;
using namespace lupnt;

void init_forces(py::module &m) {
  m.def(
      "acceleration_gravity_field",
      [](const Vec3d &r, double GM, double R, MatXd CS, int n_max, int m_max) -> Vec3d {
        return AccelarationGravityField(r, GM, R, CS, n_max, m_max).cast<double>();
      },
      py::arg("r"), py::arg("GM"), py::arg("R"), py::arg("CS"), py::arg("n_max"), py::arg("m_max"));
  m.def(
      "acceleration_point_mass",
      [](const Vec3d &r, const Vec3d &s, double GM) -> Vec3d {
        return AccelerationPointMass(r, s, GM).cast<double>();
      },
      py::arg("r"), py::arg("s"), py::arg("GM"));
  m.def(
      "acceleration_solar_radiation",
      [](const Vec3d &r, const Vec3d &r_sun, double area, double mass, double CR, double P_SUN,
         double AU) -> Vec3d {
        return AccelerationSolarRadiation(r, r_sun, area, mass, CR, P_SUN, AU).cast<double>();
      },
      py::arg("r"), py::arg("r_sun"), py::arg("area"), py::arg("mass"), py::arg("CR"),
      py::arg("P_SUN"), py::arg("AU"));
}

// lupnt
#include <lupnt/core/constants.h>
#include <lupnt/environment/forces.h>
#include <lupnt/numerics/math_utils.h>

// pybind11
#include <string>

#include "py_pybind11.h"
#include "py_vectorized_macros.h"

namespace py = pybind11;
using namespace lupnt;

void init_forces(py::module& m) {
  m.def(
      "acceleration_gravity_field",
      [](const Vec3d& r, double GM, double R, MatXd CS, int n_max, int m_max) -> Vec3d {
        return AccelarationGravityField(r, GM, R, CS, n_max, m_max).cast<double>();
      },
      py::arg("r"), py::arg("GM"), py::arg("R"), py::arg("CS"), py::arg("n_max"), py::arg("m_max"));
  m.def(
      "acceleration_point_mass",
      [](const Vec3d& r, const Vec3d& s, double GM) -> Vec3d {
        return AccelerationPointMass(r, s, GM).cast<double>();
      },
      py::arg("r"), py::arg("s"), py::arg("GM"));
  m.def(
      "acceleration_solar_radiation",
      [](const Vec3d& r, const Vec3d& r_sun, double b_srp, double P_SUN, double AU) -> Vec3d {
        return AccelerationSolarRadiation(r, r_sun, b_srp, P_SUN, AU).cast<double>();
      },
      py::arg("r"), py::arg("r_sun"), py::arg("b_srp"), py::arg("P_SUN"), py::arg("AU"));
  m.def(
      "acceleration_relativistic_correction",
      [](const Vec3d& r, const Vec3d& v, double GM, double c_light) -> Vec3d {
        return AccelerationRelativisticCorrection(r, v, GM, c_light).cast<double>();
      },
      py::arg("r"), py::arg("v"), py::arg("GM"), py::arg("c_light") = C);
  m.def(
      "shadow_function",
      [](const Vec3d& r, const Vec3d& r_sun, double R_body, double R_sun) -> double {
        return ShadowFunction(r, r_sun, R_body, R_sun).val();
      },
      py::arg("r"), py::arg("r_sun"), py::arg("R_body"), py::arg("R_sun") = R_SUN);
  m.def(
      "illumination",
      [](const Vec3d& r, const Vec3d& r_sun, double R_body) -> double {
        return Illumination(r, r_sun, R_body).val();
      },
      py::arg("r"), py::arg("r_sun"), py::arg("R_body"));
}

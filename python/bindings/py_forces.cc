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
      py::arg("r"), py::arg("GM"), py::arg("R"), py::arg("CS"), py::arg("n_max"), py::arg("m_max"),
      "Spherical-harmonic gravity acceleration [m/s^2] at body-fixed position r [m], from GM "
      "[m^3/s^2], reference radius R [m], and unnormalized coefficients CS up to degree/order "
      "n_max/m_max.");
  m.def(
      "acceleration_point_mass",
      [](const Vec3d& r, const Vec3d& s, double GM) -> Vec3d {
        return AccelerationPointMass(r, s, GM).cast<double>();
      },
      py::arg("r"), py::arg("s"), py::arg("GM"),
      "Third-body point-mass acceleration [m/s^2] on a spacecraft at r [m] from a perturbing mass "
      "at s [m] with GM [m^3/s^2] (includes the indirect central-body term).");
  m.def(
      "acceleration_solar_radiation",
      [](const Vec3d& r, const Vec3d& r_sun, double b_srp, double P_SUN, double AU) -> Vec3d {
        return AccelerationSolarRadiation(r, r_sun, b_srp, P_SUN, AU).cast<double>();
      },
      py::arg("r"), py::arg("r_sun"), py::arg("b_srp"), py::arg("P_SUN"), py::arg("AU"),
      "Cannonball solar-radiation-pressure acceleration [m/s^2] at r [m] given Sun position r_sun "
      "[m], SRP ballistic coefficient b_srp [m^2/kg], pressure P_SUN [N/m^2] at 1 AU, and AU [m].");
  m.def(
      "acceleration_relativistic_correction",
      [](const Vec3d& r, const Vec3d& v, double GM, double c_light) -> Vec3d {
        return AccelerationRelativisticCorrection(r, v, GM, c_light).cast<double>();
      },
      py::arg("r"), py::arg("v"), py::arg("GM"), py::arg("c_light") = C,
      "First-order post-Newtonian (Schwarzschild) acceleration correction [m/s^2] at position r "
      "[m], velocity v [m/s], central-body GM [m^3/s^2], and speed of light c_light [m/s].");
  m.def(
      "shadow_function",
      [](const Vec3d& r, const Vec3d& r_sun, double R_body, double R_sun) -> double {
        return ShadowFunction(r, r_sun, R_body, R_sun).val();
      },
      py::arg("r"), py::arg("r_sun"), py::arg("R_body"), py::arg("R_sun") = R_SUN,
      "Solar shadow function nu (0 in umbra, 1 in sunlight, in between in penumbra) at r [m], Sun "
      "position r_sun [m], occulting-body radius R_body [m], and Sun radius R_sun [m].");
  m.def(
      "illumination",
      [](const Vec3d& r, const Vec3d& r_sun, double R_body) -> double {
        return Illumination(r, r_sun, R_body).val();
      },
      py::arg("r"), py::arg("r_sun"), py::arg("R_body"),
      "Solar-radiation-pressure illumination factor (0 to 1) at r [m], Sun position r_sun [m], and "
      "occulting-body radius R_body [m].");
}

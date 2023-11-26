// lupnt
#include <lupnt/core/constants.h>
#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/orbit_state.h>
#include <lupnt/physics/orbit_state_utils.h>

// pybind11
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lupnt;

void init_orbit_state_utils(py::module &m) {
  //   py::class_<TLE>(m, "TLE")
  //       .def("FromLines", &TLE::FromLines)
  //       .def("FromFile", &TLE::FromFile);

  //   // functions
  //   m.def("convert_state_representation", &ConvertOrbitStateRepresentation,
  //         py::arg("state"), py::arg("to"), py::arg("mu"));

  //   m.def("convert_state_coord_system", &ConvertOrbitStateCoordSystem,
  //         py::arg("state"), py::arg("epoch"), py::arg("to"));

  // coe <-> cart
  m.def(
      "coe_to_cart",
      [](const ClassicalOE &coe, double mu) -> CartesianOrbitState {
        return CoeToCart(coe, mu);
      },
      py::arg("coe"), py::arg("mu"));
  m.def(
      "coe_to_cart",
      [](const Vector6d &coeVec, double mu) -> Vector6d {
        return CoeToCart(coeVec, mu).cast<double>();
      },
      py::arg("coe"), py::arg("mu"));
  m.def(
      "cart_to_coe",
      [](const CartesianOrbitState &cart, double mu) -> ClassicalOE {
        return CartToCoe(cart, mu);
      },
      py::arg("cart"), py::arg("mu"));
  m.def(
      "cart_to_coe",
      [](const Vector6d &cartVec, double mu) -> Vector6d {
        return CartToCoe(cartVec, mu).cast<double>();
      },
      py::arg("cart"), py::arg("mu"));

  // coe <-> roe
  m.def(
      "roe_to_coe",
      [](const ClassicalOE &coe_chief, const QuasiNonsingularROE &roe) {
        return RoeToCoe(coe_chief, roe);
      },
      py::arg("coe_chief"), py::arg("roe"));
  m.def(
      "roe_to_coe",
      [](const Vector6d &coe_chief, const Vector6d &roe) -> Vector6d {
        return RoeToCoe(coe_chief, roe).cast<double>();
      },
      py::arg("coe_chief"), py::arg("roe"));

  // inertial <-> rtn
  m.def(
      "inertial_to_rtn",
      [](const CartesianOrbitState &rv_chief,
         const CartesianOrbitState &rv_deputy) {
        return InertialToRtn(rv_chief, rv_deputy);
      },
      py::arg("cart_orig"), py::arg("cart"));
  m.def(
      "inertial_to_rtn",
      [](const Vector6d &cart_orig, const Vector6d &cart) -> Vector6d {
        return InertialToRtn(cart_orig, cart).cast<double>();
      },
      py::arg("cart_orig"), py::arg("cart"));

  m.def(
      "coe_to_rtn",
      [](const ClassicalOE &coe_chief, const ClassicalOE &coe_deputy,
         double mu) { return CoeToRtn(coe_chief, coe_deputy, mu); },
      py::arg("coe_chief"), py::arg("coe_deputy"), py::arg("mu"));
  m.def(
      "coe_to_rtn",
      [](const Vector6d &coe_chief, const Vector6d &coe_deputy,
         double mu) -> Vector6d {
        return CoeToRtn(coe_chief, coe_deputy, mu).cast<double>();
      },
      py::arg("coe_chief"), py::arg("coe_deputy"), py::arg("mu"));

  // Anomaly Conversions
  m.def("eccentric_to_true", [](double E, double e) -> double {
    return EccentricAnomToTrueAnom(E, e).val();
  });
  m.def("eccentric_to_mean", [](double E, double e) -> double {
    return EccentricAnomToMeanAnom(E, e).val();
  });
  m.def("mean_to_eccentric", [](double M, double e) -> double {
    return MeanAnomToEccentricAnom(M, e).val();
  });
  m.def("mean_to_true", [](double M, double e) -> double {
    return MeanAnomToTrueAnom(M, e).val();
  });
  m.def("true_to_eccentric", [](double nu, double e) -> double {
    return TrueAnomToEccentricAnom(nu, e).val();
  });
  m.def("true_to_mean", [](double f, double e) -> double {
    return TrueAnomToMeanAnom(f, e).val();
  });
}
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
  pnt.coe_to_eqoe
      // coe <-> qnsoe
      m.def(
          "coe_to_qnsoe",
          [](const ClassicalOE &coe) -> QuasiNonsingularOE {
            return CoeToQnsoe(coe);
          },
          py::arg("coe"));
  m.def(
      "coe_to_qnsoe",
      [](const Vector6d &coeVec) -> Vector6d {
        return CoeToQnsoe(coeVec).cast<double>();
      },
      py::arg("coe"));
  m.def(
      "qnsoe_to_coe",
      [](const QuasiNonsingularOE &qnsoe) -> ClassicalOE {
        return QnsoeToCoe(qnsoe);
      },
      py::arg("qnsoe"));
  m.def(
      "qnsoe_to_coe",
      [](const Vector6d &qnsoeVec) -> Vector6d {
        return QnsoeToCoe(qnsoeVec).cast<double>();
      },
      py::arg("qnsoe"));

  // coe <-> roe
  m.def(
      "qnsroe_to_coe",
      [](const ClassicalOE &coe_chief, const QuasiNonsingularROE &roe) {
        return QnsroeToCoe(coe_chief, roe);
      },
      py::arg("coe_chief"), py::arg("roe"));
  m.def(
      "qnsroe_to_coe",
      [](const Vector6d &coe_chief, const Vector6d &roe) -> Vector6d {
        return QnsroeToCoe(coe_chief, roe).cast<double>();
      },
      py::arg("coe_chief"), py::arg("roe"));

  // coe <-> eqoe
  m.def(
      "coe_to_eqoe",
      [](const ClassicalOE &coe) -> EquinoctialOE { return CoeToEqoe(coe); },
      py::arg("coe"));
  m.def(
      "coe_to_eqoe",
      [](const Vector6d &coeVec) -> Vector6d {
        return CoeToEqoe(coeVec).cast<double>();
      },
      py::arg("coe"));
  m.def(
      "eqoe_to_coe",
      [](const EquinoctialOE &eqoe) -> ClassicalOE { return EqoeToCoe(eqoe); },
      py::arg("eqoe"));
  m.def(
      "eqoe_to_coe",
      [](const Vector6d &eqoeVec) -> Vector6d {
        return EqoeToCoe(eqoeVec).cast<double>();
      },
      py::arg("eqoe"));

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
    return EccentricToTrueAnomaly(E, e).val();
  });
  m.def("eccentric_to_mean", [](double E, double e) -> double {
    return EccentricToMeanAnomaly(E, e).val();
  });
  m.def("mean_to_eccentric", [](double M, double e) -> double {
    return MeanToEccentricAnomaly(M, e).val();
  });
  m.def("mean_to_true", [](double M, double e) -> double {
    return MeanToTrueAnomaly(M, e).val();
  });
  m.def("true_to_eccentric", [](double nu, double e) -> double {
    return TrueToEccentricAnomaly(nu, e).val();
  });
  m.def("true_to_mean", [](double f, double e) -> double {
    return TrueToMeanAnomaly(f, e).val();
  });
}
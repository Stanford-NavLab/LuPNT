// lupnt
#include <lupnt/core/constants.h>
#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/coordinates.h>
#include <lupnt/physics/orbit_state.h>

// pybind11
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include <string>

#include "py_vectorized_macros.cc"

namespace py = pybind11;
using namespace lupnt;

void init_orbit_state_utils(py::module &m) {
  // Orbit State Conversions
  m.def("convert_orbit_state",
        [](const Vec6d &state_in, OrbitStateRepres repres_in, OrbitStateRepres repres_out,
           double GM) -> Vec6d {
          return ConvertOrbitState(state_in, repres_in, repres_out, GM).cast<double>();
        });
  m.def("convert_orbit_state",
        [](const Vec6d &state_in_c, const Vec6 &state_in_d, OrbitStateRepres repres_in_c,
           OrbitStateRepres repres_in_d, OrbitStateRepres repres_out, double GM) -> Vec6d {
          return ConvertOrbitState(state_in_c, state_in_d, repres_in_c, repres_in_d, repres_out, GM)
              .cast<double>();
        });

  m.def(
      "cartesian_to_classical",
      [](const CartesianOrbitState &cart, double GM) -> ClassicalOE {
        return Cart2Classical(cart, GM);
      },
      py::arg("cart"), py::arg("GM"));
  m.def(
      "classical_to_cartesian",
      [](const ClassicalOE &coe, double GM) -> CartesianOrbitState {
        return Classical2Cart(coe, GM);
      },
      py::arg("coe"), py::arg("GM"));
  m.def(
      "classical_to_quasi_nonsingular",
      [](const ClassicalOE &coe, double GM) -> QuasiNonsingOE {
        return Classical2QuasiNonsing(coe, GM);
      },
      py::arg("coe"), py::arg("GM"));
  m.def("classical_to_equinoctial", [](const ClassicalOE &coe, double GM) -> EquinoctialOE {
    return Classical2Equinoctial(coe, GM);
  });
  m.def("classical_to_delaunay", [](const ClassicalOE &coe, double GM) -> DelaunayOE {
    return Classical2Delaunay(coe, GM);
  });
  m.def("quasi_nonsingular_to_classical",
        [](const QuasiNonsingOE &qnsoe, double GM) -> ClassicalOE {
          return QuasiNonsing2Classical(qnsoe, GM);
        });
  m.def(
      "equinoctial_to_classical",
      [](const EquinoctialOE &eqoe, double GM) -> ClassicalOE {
        return Equinoctial2Classical(eqoe, GM);
      },
      py::arg("eqoe"), py::arg("GM"));
  m.def(
      "delaunay_to_classical",
      [](const DelaunayOE &deloe, double GM) -> ClassicalOE {
        return Delaunay2Classical(deloe, GM);
      },
      py::arg("deloe"), py::arg("GM"));
  m.def(
      "relative_quasi_nonsingular_to_classical",
      [](const ClassicalOE &coe, const QuasiNonsingROE &rel_qnsoe) -> ClassicalOE {
        return RelQuasiNonsing2Classical(coe, rel_qnsoe);
      },
      py::arg("coe"), py::arg("rel_qnsoe"));

  // State Conversions
  VECTORIZED_BINDING_FROM_VECTOR_REAL("cartesian_to_classical", Cart2Classical, 6, "cart", "GM");
  VECTORIZED_BINDING_FROM_VECTOR_REAL("classical_to_cartesian", Classical2Cart, 6, "coe", "GM");
  VECTORIZED_BINDING_FROM_VECTOR_REAL("classical_to_quasi_nonsingular", Classical2QuasiNonsing, 6,
                                      "coe", "GM");
  VECTORIZED_BINDING_FROM_VECTOR_REAL("classical_to_equinoctial", Classical2Equinoctial, 6, "coe",
                                      "GM");
  VECTORIZED_BINDING_FROM_VECTOR_REAL("classical_to_delaunay", Classical2Delaunay, 6, "coe", "GM");
  VECTORIZED_BINDING_FROM_VECTOR_REAL("quasi_nonsingular_to_classical", QuasiNonsing2Classical, 6,
                                      "qnsoe", "GM");
  VECTORIZED_BINDING_FROM_VECTOR_REAL("equinoctial_to_classical", Equinoctial2Classical, 6, "eqoe",
                                      "GM");
  VECTORIZED_BINDING_FROM_VECTOR_REAL("delaunay_to_classical", Delaunay2Classical, 6, "deloe",
                                      "GM");
  VECTORIZED_BINDING_FROM_VECTOR_VECTOR("relative_quasi_nonsingular_to_classical",
                                        RelQuasiNonsing2Classical, 6, "coe", "rel_qnsoe");

  // Anomaly Conversions
  VECTORIZED_BINDING_FROM_REAL_REAL("eccentric_to_true_anomaly", Ecc2TrueAnomaly, "E", "e");
  VECTORIZED_BINDING_FROM_REAL_REAL("eccentric_to_mean_anomaly", Ecc2MeanAnomaly, "E", "e");
  VECTORIZED_BINDING_FROM_REAL_REAL("mean_to_eccentric_anomaly", Mean2EccAnomaly, "M", "e");
  VECTORIZED_BINDING_FROM_REAL_REAL("mean_to_true_anomaly", Mean2TrueAnomaly, "M", "e");
  VECTORIZED_BINDING_FROM_REAL_REAL("true_to_eccentric_anomaly", True2EccAnomaly, "nu", "e");
  VECTORIZED_BINDING_FROM_REAL_REAL("true_to_mean_anomaly", True2MeanAnomaly, "f", "e");

  // Coordinate System Conversions
  // VECTORIZED_BINDING_FROM_VECTOR("LatLonAlt2Cart", LatLonAlt2Cart, 3, "r_geo");
  m.def(
      "LatLonAlt2Cart",
      [](const Vecd<3> &x) -> Vecd<3> {
        return LatLonAlt2Cart(x.cast<Real>().eval()).cast<double>();
      },
      py::arg("r_geo"));
  m.def(
      "LatLonAlt2Cart",
      [](const Matd<-1, 3> &x) -> Matd<-1, 3> {
        return LatLonAlt2Cart(x.cast<Real>().eval()).cast<double>();
      },
      py::arg("r_geo"));

  // VECTORIZED_BINDING_FROM_VECTOR("Cart2LatLonAlt", Cart2LatLonAlt, 3, "r_cart");
  m.def(
      "Cart2LatLonAlt",
      [](const Vecd<3> &x) -> Vecd<3> {
        return Cart2LatLonAlt(x.cast<Real>().eval()).cast<double>();
      },
      py::arg("r_cart"));
  m.def(
      "Cart2LatLonAlt",
      [](const Matd<-1, 3> &x) -> Matd<-1, 3> {
        return Cart2LatLonAlt(x.cast<Real>().eval()).cast<double>();
      },
      py::arg("r_cart"));

  VECTORIZED_BINDING_FROM_VECTOR_VECTOR("EastNorthUp2Cart", EastNorthUp2Cart, 3, "enu", "xyz_ref");
  VECTORIZED_BINDING_FROM_VECTOR_VECTOR("Cart2EastNorthUp", Cart2EastNorthUp, 3, "xzy", "xyz_ref");
  VECTORIZED_BINDING_FROM_VECTOR_VECTOR("Cart2AzElRange", Cart2AzElRange, 3, "xzy", "xyz_ref");
  VECTORIZED_BINDING_FROM_VECTOR_VECTOR("AzElRange2Cart", AzElRange2Cart, 3, "aer", "xyz_ref");
}
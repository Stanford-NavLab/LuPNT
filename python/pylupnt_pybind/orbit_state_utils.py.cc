// lupnt
#include <lupnt/core/constants.h>
#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/orbit_state.h>
#include <lupnt/physics/orbit_state_utils.h>

// pybind11
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include <string>

#include "vectorized_macros.py.cc"

namespace py = pybind11;
using namespace lupnt;

void init_orbit_state_utils(py::module &m) {
  // Orbit State Conversions
  m.def("convert_orbit_state",
        [](const Vec6d &state_in, OrbitStateRepres repres_in,
           OrbitStateRepres repres_out, double mu) -> Vec6d {
          return ConvertOrbitState(state_in, repres_in, repres_out, mu)
              .cast<double>();
        });
  m.def("convert_orbit_state",
        [](const Vec6d &state_in_c, const Vec6 &state_in_d,
           OrbitStateRepres repres_in_c, OrbitStateRepres repres_in_d,
           OrbitStateRepres repres_out, double mu) -> Vec6d {
          return ConvertOrbitState(state_in_c, state_in_d, repres_in_c,
                                   repres_in_d, repres_out, mu)
              .cast<double>();
        });

  m.def(
      "cartesian_to_classical",
      [](const CartesianOrbitState &cart, double mu) -> ClassicalOE {
        return CartesianToClassical(cart, mu);
      },
      py::arg("cart"), py::arg("mu"));
  m.def(
      "classical_to_cartesian",
      [](const ClassicalOE &coe, double mu) -> CartesianOrbitState {
        return ClassicalToCartesian(coe, mu);
      },
      py::arg("coe"), py::arg("mu"));
  m.def(
      "classical_to_quasi_nonsingular",
      [](const ClassicalOE &coe, double mu) -> QuasiNonsingularOE {
        return ClassicalToQuasiNonsingular(coe, mu);
      },
      py::arg("coe"), py::arg("mu"));
  m.def("classical_to_equinoctial",
        [](const ClassicalOE &coe, double mu) -> EquinoctialOE {
          return ClassicalToEquinoctial(coe, mu);
        });
  m.def("classical_to_delaunay",
        [](const ClassicalOE &coe, double mu) -> DelaunayOE {
          return ClassicalToDelaunay(coe, mu);
        });
  m.def("quasi_nonsingular_to_classical",
        [](const QuasiNonsingularOE &qnsoe, double mu) -> ClassicalOE {
          return QuasiNonsingularToClassical(qnsoe, mu);
        });
  m.def(
      "equinoctial_to_classical",
      [](const EquinoctialOE &eqoe, double mu) -> ClassicalOE {
        return EquinoctialToClassical(eqoe, mu);
      },
      py::arg("eqoe"), py::arg("mu"));
  m.def(
      "delaunay_to_classical",
      [](const DelaunayOE &deloe, double mu) -> ClassicalOE {
        return DelaunayToClassical(deloe, mu);
      },
      py::arg("deloe"), py::arg("mu"));
  m.def(
      "relative_quasi_nonsingular_to_classical",
      [](const ClassicalOE &coe,
         const QuasiNonsingularROE &rel_qnsoe) -> ClassicalOE {
        return RelativeQuasiNonsingularToClassical(coe, rel_qnsoe);
      },
      py::arg("coe"), py::arg("rel_qnsoe"));

  // State Conversions
  VECTORIZED_BINDING_FROM_VECTOR_REAL("cartesian_to_classical",
                                      CartesianToClassical, 6, "cart", "mu");
  VECTORIZED_BINDING_FROM_VECTOR_REAL("classical_to_cartesian",
                                      ClassicalToCartesian, 6, "coe", "mu");
  VECTORIZED_BINDING_FROM_VECTOR_REAL("classical_to_quasi_nonsingular",
                                      ClassicalToQuasiNonsingular, 6, "coe",
                                      "mu");
  VECTORIZED_BINDING_FROM_VECTOR_REAL("classical_to_equinoctial",
                                      ClassicalToEquinoctial, 6, "coe", "mu");
  VECTORIZED_BINDING_FROM_VECTOR_REAL("classical_to_delaunay",
                                      ClassicalToDelaunay, 6, "coe", "mu");
  VECTORIZED_BINDING_FROM_VECTOR_REAL("quasi_nonsingular_to_classical",
                                      QuasiNonsingularToClassical, 6, "qnsoe",
                                      "mu");
  VECTORIZED_BINDING_FROM_VECTOR_REAL("equinoctial_to_classical",
                                      EquinoctialToClassical, 6, "eqoe", "mu");
  VECTORIZED_BINDING_FROM_VECTOR_REAL("delaunay_to_classical",
                                      DelaunayToClassical, 6, "deloe", "mu");
  VECTORIZED_BINDING_FROM_VECTOR_VECTOR(
      "relative_quasi_nonsingular_to_classical",
      RelativeQuasiNonsingularToClassical, 6, "coe", "rel_qnsoe");

  // Anomaly Conversions
  VECTORIZED_BINDING_FROM_REAL_REAL("eccentric_to_true_anomaly",
                                    EccentricToTrueAnomaly, "E", "e");
  VECTORIZED_BINDING_FROM_REAL_REAL("eccentric_to_mean_anomaly",
                                    EccentricToMeanAnomaly, "E", "e");
  VECTORIZED_BINDING_FROM_REAL_REAL("mean_to_eccentric_anomaly",
                                    MeanToEccentricAnomaly, "M", "e");
  VECTORIZED_BINDING_FROM_REAL_REAL("mean_to_true_anomaly", MeanToTrueAnomaly,
                                    "M", "e");
  VECTORIZED_BINDING_FROM_REAL_REAL("true_to_eccentric_anomaly",
                                    TrueToEccentricAnomaly, "nu", "e");
  VECTORIZED_BINDING_FROM_REAL_REAL("true_to_mean_anomaly", TrueToMeanAnomaly,
                                    "f", "e");

  // Coordinate System Conversions
  VECTORIZED_BINDING_FROM_VECTOR_REAL("geographical_to_cartesian",
                                      LatLonAltToEcef, 3, "r_geo", "radius");
  VECTORIZED_BINDING_FROM_VECTOR_REAL("cartesian_to_geographical",
                                      EcefToLatLonAlt, 3, "r_cart", "radius");
  VECTORIZED_BINDING_FROM_VECTOR("spherical_to_cartesian", SphericalToCartesian,
                                 3, "r_sph");
  VECTORIZED_BINDING_FROM_VECTOR("cartesian_to_spherical", CartesianToSpherical,
                                 3, "r_cart");
  VECTORIZED_BINDING_FROM_VECTOR_VECTOR("east_north_up_to_cartesian",
                                        EastNorthUpToCartesian, 3, "r_ref",
                                        "r_enu");
  VECTORIZED_BINDING_FROM_VECTOR_VECTOR("cartesian_to_east_north_up",
                                        CartesianToEastNorthUp, 3, "r_ref",
                                        "r_cart");
  VECTORIZED_BINDING_FROM_VECTOR_VECTOR("cartesian_to_azimuth_elevation_range",
                                        CartesianToAzimuthElevationRange, 3,
                                        "r_cart_ref", "r_cart");
  VECTORIZED_BINDING_FROM_VECTOR_VECTOR("azimuth_elevation_range_to_cartesian",
                                        AzimuthElevationRangeToCartesian, 3,
                                        "r_aer_ref", "r_aer");
}
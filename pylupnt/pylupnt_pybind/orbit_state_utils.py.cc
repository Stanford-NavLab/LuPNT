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
        [](const Vector6d &state_in, OrbitStateRepres repres_in,
           OrbitStateRepres repres_out, double mu) -> Vector6d {
          return ConvertOrbitState(state_in, repres_in, repres_out, mu)
              .cast<double>();
        });
  m.def("convert_orbit_state",
        [](const Vector6d &state_in_c, const Vector6 &state_in_d,
           OrbitStateRepres repres_in_c, OrbitStateRepres repres_in_d,
           OrbitStateRepres repres_out, double mu) -> Vector6d {
          return ConvertOrbitState(state_in_c, state_in_d, repres_in_c,
                                   repres_in_d, repres_out, mu)
              .cast<double>();
        });

  m.def(
      "cartesian_to_classical",
      [](const Vector6d &cart, double mu) -> Vector6d {
        return CartesianToClassical(cart, mu).cast<double>();
      },
      py::arg("cart"), py::arg("mu") = MU_MOON);
  m.def(
      "cartesian_to_classical",
      [](const CartesianOrbitState &cart, double mu) -> ClassicalOE {
        return CartesianToClassical(cart, mu);
      },
      py::arg("cart"), py::arg("mu") = MU_MOON);

  m.def(
      "classical_to_cartesian",
      [](const Vector6d &coe, double mu) -> Vector6d {
        return ClassicalToCartesian(coe, mu).cast<double>();
      },
      py::arg("coe"), py::arg("mu") = MU_MOON);
  m.def(
      "classical_to_cartesian",
      [](const ClassicalOE &coe, double mu) -> CartesianOrbitState {
        return ClassicalToCartesian(coe, mu);
      },
      py::arg("coe"), py::arg("mu") = MU_MOON);

  m.def(
      "classical_to_quasi_nonsingular",
      [](const Vector6d &coe, double mu) -> Vector6d {
        return ClassicalToQuasiNonsingular(coe, mu).cast<double>();
      },
      py::arg("coe"), py::arg("mu") = MU_MOON);
  m.def(
      "classical_to_quasi_nonsingular",
      [](const ClassicalOE &coe, double mu) -> QuasiNonsingularOE {
        return ClassicalToQuasiNonsingular(coe, mu);
      },
      py::arg("coe"), py::arg("mu") = MU_MOON);

  m.def(
      "classical_to_equinoctial",
      [](const Vector6d &coe, double mu) -> Vector6d {
        return ClassicalToEquinoctial(coe, mu).cast<double>();
      },
      py::arg("coe"), py::arg("mu") = MU_MOON);
  m.def("classical_to_equinoctial",
        [](const ClassicalOE &coe, double mu) -> EquinoctialOE {
          return ClassicalToEquinoctial(coe, mu);
        });

  m.def("classical_to_delaunay",
        [](const Vector6d &coe, double mu) -> Vector6d {
          return ClassicalToDelaunay(coe, mu).cast<double>();
        });
  m.def("classical_to_delaunay",
        [](const ClassicalOE &coe, double mu) -> DelaunayOE {
          return ClassicalToDelaunay(coe, mu);
        });

  m.def("quasi_nonsingular_to_classical",
        [](const Vector6d &qnsoe, double mu) -> Vector6d {
          return QuasiNonsingularToClassical(qnsoe, mu).cast<double>();
        });
  m.def("quasi_nonsingular_to_classical",
        [](const QuasiNonsingularOE &qnsoe, double mu) -> ClassicalOE {
          return QuasiNonsingularToClassical(qnsoe, mu);
        });

  m.def(
      "equinoctial_to_classical",
      [](const Vector6d &eqoe, double mu) -> Vector6d {
        return EquinoctialToClassical(eqoe, mu).cast<double>();
      },
      py::arg("eqoe"), py::arg("mu") = MU_MOON);
  m.def(
      "equinoctial_to_classical",
      [](const EquinoctialOE &eqoe, double mu) -> ClassicalOE {
        return EquinoctialToClassical(eqoe, mu);
      },
      py::arg("eqoe"), py::arg("mu") = MU_MOON);

  m.def("delaunay_to_classical",
        [](const Vector6d &deloe, double mu) -> Vector6d {
          return DelaunayToClassical(deloe, mu).cast<double>();
        });
  m.def(
      "delaunay_to_classical",
      [](const DelaunayOE &deloe, double mu) -> ClassicalOE {
        return DelaunayToClassical(deloe, mu);
      },
      py::arg("deloe"), py::arg("mu") = MU_MOON);

  m.def(
      "relative_quasi_nonsingular_to_classical",
      [](const Vector6d &coe, const Vector6d &rel_qnsoe) -> Vector6d {
        return RelativeQuasiNonsingularToClassical(coe, rel_qnsoe)
            .cast<double>();
      },
      py::arg("coe"), py::arg("rel_qnsoe"));
  m.def(
      "relative_quasi_nonsingular_to_classical",
      [](const ClassicalOE &coe,
         const QuasiNonsingularROE &rel_qnsoe) -> ClassicalOE {
        return RelativeQuasiNonsingularToClassical(coe, rel_qnsoe);
      },
      py::arg("coe"), py::arg("rel_qnsoe"));

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
                                      GeographicalToCartesian, 3, "r_geo",
                                      "radius");
  VECTORIZED_BINDING_FROM_VECTOR_REAL("cartesian_to_geographical",
                                      CartesianToGeographical, 3, "r_cart",
                                      "radius");
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
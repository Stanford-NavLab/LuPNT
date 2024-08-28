// lupnt
#include <lupnt/core/constants.h>
#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/coordinates.h>
#include <lupnt/physics/orbit_state.h>
#include <lupnt/physics/orbit_state/anomaly.h>

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
      "cart2classical",
      [](const CartesianOrbitState &cart, double GM) -> ClassicalOE {
        return Cart2Classical(cart, GM);
      },
      py::arg("cart"), py::arg("GM"));
  m.def(
      "classical2cart",
      [](const ClassicalOE &coe, double GM) -> CartesianOrbitState {
        return Classical2Cart(coe, GM);
      },
      py::arg("coe"), py::arg("GM"));
  m.def(
      "classical2quasi_nonsingular",
      [](const ClassicalOE &coe, double GM) -> QuasiNonsingOE {
        return Classical2QuasiNonsing(coe, GM);
      },
      py::arg("coe"), py::arg("GM"));
  m.def("classical2equinoctial", [](const ClassicalOE &coe, double GM) -> EquinoctialOE {
    return Classical2Equinoctial(coe, GM);
  });
  m.def("classical2delaunay", [](const ClassicalOE &coe, double GM) -> DelaunayOE {
    return Classical2Delaunay(coe, GM);
  });
  m.def("quasi_nonsingular2classical", [](const QuasiNonsingOE &qnsoe, double GM) -> ClassicalOE {
    return QuasiNonsing2Classical(qnsoe, GM);
  });
  m.def(
      "equinoctial2classical",
      [](const EquinoctialOE &eqoe, double GM) -> ClassicalOE {
        return Equinoctial2Classical(eqoe, GM);
      },
      py::arg("eqoe"), py::arg("GM"));
  m.def(
      "delaunay2classical",
      [](const DelaunayOE &deloe, double GM) -> ClassicalOE {
        return Delaunay2Classical(deloe, GM);
      },
      py::arg("deloe"), py::arg("GM"));
  m.def(
      "relative_quasi_nonsingular2classical",
      [](const ClassicalOE &coe, const QuasiNonsingROE &rel_qnsoe) -> ClassicalOE {
        return RelQuasiNonsing2Classical(coe, rel_qnsoe);
      },
      py::arg("coe"), py::arg("rel_qnsoe"));

  VEC_BIND_REAL_REAL("get_orbital_period", GetOrbitalPeriod, "a", "GM");

  // State Conversions
  VEC_BIND_VECTOR_REAL("cart2classical", Cart2Classical, 6, "cart", "GM");
  VEC_BIND_VECTOR_REAL("classical2cart", Classical2Cart, 6, "coe", "GM");
  VEC_BIND_VECTOR_REAL("classical2quasi_nonsingular", Classical2QuasiNonsing, 6, "coe", "GM");
  VEC_BIND_VECTOR_REAL("classical2equinoctial", Classical2Equinoctial, 6, "coe", "GM");
  VEC_BIND_VECTOR_REAL("classical2delaunay", Classical2Delaunay, 6, "coe", "GM");
  VEC_BIND_VECTOR_REAL("quasi_nonsingular2classical", QuasiNonsing2Classical, 6, "qnsoe", "GM");
  VEC_BIND_VECTOR_REAL("equinoctial2classical", Equinoctial2Classical, 6, "eqoe", "GM");
  VEC_BIND_VECTOR_REAL("delaunay2classical", Delaunay2Classical, 6, "deloe", "GM");
  VEC_BIND_VECTOR_VECTOR("relative_quasi_nonsingular2classical", RelQuasiNonsing2Classical, 6,
                         "coe", "rel_qnsoe");

  // Anomaly Conversions
  VEC_BIND_REAL_REAL("eccentric2true_anomaly", Ecc2TrueAnomaly, "E", "e");
  VEC_BIND_REAL_REAL("eccentric2mean_anomaly", Ecc2MeanAnomaly, "E", "e");
  VEC_BIND_REAL_REAL("mean2eccentric_anomaly", Mean2EccAnomaly, "M", "e");
  VEC_BIND_REAL_REAL("mean2true_anomaly", Mean2TrueAnomaly, "M", "e");
  VEC_BIND_REAL_REAL("true2eccentric_anomaly", True2EccAnomaly, "nu", "e");
  VEC_BIND_REAL_REAL("true2mean_anomaly", True2MeanAnomaly, "f", "e");

  // Coordinate System Conversions
  VEC_BIND_VECTOR("lat_lon_alt2cart", LatLonAlt2Cart, 3, "r_geo");
  VEC_BIND_VECTOR("cart2lat_lon_alt", Cart2LatLonAlt, 3, "r_cart");
  VEC_BIND_VECTOR_REAL("lat_lon_alt2cart", LatLonAlt2Cart, 3, "r_geo", "R_body");
  VEC_BIND_VECTOR_REAL("cart2lat_lon_alt", Cart2LatLonAlt, 3, "r_cart", "R_body");

  VEC_BIND_VECTOR_VECTOR("east_north_up2cart", EastNorthUp2Cart, 3, "enu", "xyz_ref");
  VEC_BIND_VECTOR_VECTOR("cart2east_north_up", Cart2EastNorthUp, 3, "xzy", "xyz_ref");
  VEC_BIND_VECTOR_VECTOR("cart2az_el_range", Cart2AzElRange, 3, "xzy", "xyz_ref");
  VEC_BIND_VECTOR_VECTOR("az_el_range2cart", AzElRange2Cart, 3, "aer", "xyz_ref");
  VEC_BIND_VECTOR_VECTOR("synodic2intertial", Synodic2Intertial, 6, "rv_c", "rv_d");
  VEC_BIND_VECTOR_VECTOR("inertial2synodic", Inertial2Synodic, 6, "rv_c", "rv_d");
}

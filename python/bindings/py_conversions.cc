#include <lupnt/lupnt.h>

#include <string>

#include "py_pybind11.h"
#include "py_vectorized_macros.h"

namespace py = pybind11;
using namespace lupnt;

void InitConversions(py::module& m) {
  // Anomaly conversions
  m.DEF_REAL_REAL("get_orbital_period", GetOrbitalPeriod, "a", "GM");
  m.DEF_REAL_REAL("eccentric_to_true_anomaly", EccToTrueAnomaly, "E", "e");
  m.DEF_REAL_REAL("eccentric_to_mean_anomaly", EccToMeanAnomaly, "E", "e");
  m.DEF_REAL_REAL("mean_to_eccentric_anomaly", MeanToEccAnomaly, "M", "e");
  m.DEF_REAL_REAL("mean_to_true_anomaly", MeanToTrueAnomaly, "M", "e");
  m.DEF_REAL_REAL("true_to_eccentric_anomaly", TrueToEccAnomaly, "nu", "e");
  m.DEF_REAL_REAL("true_to_mean_anomaly", TrueToMeanAnomaly, "f", "e");

  // Attitude conversions
  m.DEF_VECTOR("scalar_first_to_last", ScalarFirstToLast, 4, "qw");
  m.DEF_VECTOR("scalar_last_to_first", ScalarLastToFirst, 4, "qw");
  m.DEF_VECTOR("normalize_quat", NormalizeQuat, 4, "qw");
  m.def("quat_to_rot", &QuatToRot, py::arg("qw"));
  m.def("rot_to_quat", &RotToQuat, py::arg("rot"));
  m.def("roll_pitch_yaw_to_rot", &RollPitchYawToRot, py::arg("rpy"));
  m.def("rot_to_roll_pitch_yaw", &RotToRollPitchYaw, py::arg("rot"));

  // Coordinate conversions
  m.DEF_VECTOR_REAL("lat_lon_alt_to_cart", LatLonAltToCart, 3, "lla", "R_body");
  m.DEF_VECTOR_REAL("cart_to_lat_lon_alt", CartToLatLonAlt, 3, "xyz", "R_body");
  m.DEF_VECTOR("east_north_up_to_az_el_range", EastNorthUpToAzElRange, 3, "enu");
  m.DEF_VECTOR("az_el_range_to_east_north_up", AzElRangeToEastNorthUp, 3, "aer");
  m.DEF_VECTOR_VECTOR("east_north_up_to_cart", EastNorthUpToCart, 3, "enu", "xyz_ref");
  m.DEF_VECTOR_VECTOR("cart_to_east_north_up", CartToEastNorthUp, 3, "xyz", "xyz_ref");
  m.DEF_VECTOR_VECTOR("cart_to_az_el_range", CartToAzElRange, 3, "xyz", "xyz_ref");
  m.DEF_VECTOR_VECTOR("az_el_range_to_cart", AzElRangeToCart, 3, "aer", "xyz_ref");
  m.DEF_VECTOR_REAL("lat_lon_alt_to_stereographic", LatLonAltToStereographic, 3, "lla", "R_body");
  m.DEF_VECTOR_REAL("stereographic_to_lat_lon_alt", StereographicToLatLonAlt, 3, "xya", "R_body");
  m.DEF_VECTOR_REAL("stereographic_to_cart", StereographicToCart, 3, "xya", "R_body");
  m.DEF_VECTOR_REAL("cart_to_stereographic", CartToStereographic, 3, "xyz", "R_body");

  // State conversions
  m.DEF_VECTOR_REAL("classical_to_cart", ClassicalToCart, 6, "coe", "GM");
  m.DEF_VECTOR_REAL("cart_to_classical", CartToClassical, 6, "cart", "GM");
  m.DEF_VECTOR_REAL("classical_to_quasi_nonsingular", ClassicalToQuasiNonsing, 6, "coe", "GM");
  m.DEF_VECTOR_REAL("quasi_nonsingular_to_classical", QuasiNonsingToClassical, 6, "qnsoe", "GM");
  m.DEF_VECTOR_REAL("classical_to_equinoctial", ClassicalToEquinoctial, 6, "coe", "GM");
  m.DEF_VECTOR_REAL("equinoctial_to_classical", EquinoctialToClassical, 6, "eqoe", "GM");
  m.DEF_VECTOR_REAL("classical_to_delaunay", ClassicalToDelaunay, 6, "coe", "GM");
  m.DEF_VECTOR_REAL("delaunay_to_classical", DelaunayToClassical, 6, "deloe", "GM");
  m.DEF_VECTOR_REAL("classical_to_delaunay", ClassicalToDelaunay, 6, "coe", "GM");
  m.DEF_VECTOR_REAL("delaunay_to_classical", DelaunayToClassical, 6, "deloe", "GM");
  m.DEF_VECTOR_REAL("classical_to_delaunay", ClassicalToDelaunay, 6, "coe", "GM");
  m.DEF_VECTOR_REAL("delaunay_to_classical", DelaunayToClassical, 6, "deloe", "GM");
  m.DEF_VECTOR_REAL("classical_to_delaunay", ClassicalToDelaunay, 6, "coe", "GM");
  m.DEF_VECTOR_REAL("delaunay_to_classical", DelaunayToClassical, 6, "deloe", "GM");
  m.DEF_VECTOR_REAL("quasi_nonsingular_to_classical", QuasiNonsingToClassical, 6, "qnsoe", "GM");
  m.DEF_VECTOR_REAL("equinoctial_to_classical", EquinoctialToClassical, 6, "eqoe", "GM");
  m.DEF_VECTOR_REAL("delaunay_to_classical", DelaunayToClassical, 6, "deloe", "GM");
  m.DEF_VECTOR_REAL("classical_to_delaunay", ClassicalToDelaunay, 6, "coe", "GM");
  m.DEF_VECTOR_REAL("delaunay_to_classical", DelaunayToClassical, 6, "deloe", "GM");
  m.def("tle_to_classical", &TleToClassical, "tle", "GM");

  // State converter
  // m.def("convert_state",
  //       py::overload_cast<const State &, StateType, StateType, Real>(&ConvertState), "state",
  //       "repres_in", "repres_out", "GM");
  // m.def("convert_state",
  //       py::overload_cast<const MatX6 &, StateType, StateType, Real>(&ConvertState), "state",
  //       "repres_in", "repres_out", "GM");
}

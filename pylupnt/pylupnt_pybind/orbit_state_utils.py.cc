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
  m.def("eccentric_to_true_anomaly", [](double E, double e) -> double {
    return EccentricToTrueAnomaly(E, e).val();
  });
  m.def("eccentric_to_mean_anomaly", [](double E, double e) -> double {
    return EccentricToMeanAnomaly(E, e).val();
  });
  m.def("mean_to_eccentric_anomaly", [](double M, double e) -> double {
    return MeanToEccentricAnomaly(M, e).val();
  });
  m.def("mean_to_true_anomaly", [](double M, double e) -> double {
    return MeanToTrueAnomaly(M, e).val();
  });
  m.def("true_to_eccentric_anomaly", [](double nu, double e) -> double {
    return TrueToEccentricAnomaly(nu, e).val();
  });
  m.def("true_to_mean_anomaly", [](double f, double e) -> double {
    return TrueToMeanAnomaly(f, e).val();
  });

  m.def("geographical_to_cartesian",
        [](const Vector3d &r_geo, double radius) -> Vector3d {
          return GeographicalToCartesian(r_geo, radius).cast<double>();
        });
  m.def("cartesian_to_geographical",
        [](const Vector3d &r_cart, double radius) -> Vector3d {
          return CartesianToGeographical(r_cart, radius).cast<double>();
        });
  m.def("spherical_to_cartesian", [](const Vector3d &r_sph) -> Vector3d {
    return SphericalToCartesian(r_sph).cast<double>();
  });
  m.def("cartesian_to_spherical", [](const Vector3d &r_cart) -> Vector3d {
    return CartesianToSpherical(r_cart).cast<double>();
  });
  m.def("east_north_up_to_cartesian",
        [](const Vector3d &r_ref, const Vector3d &r_enu) -> Vector3d {
          return EastNortUpToCartesian(r_ref, r_enu).cast<double>();
        });
  m.def("cartesian_to_east_north_up",
        [](const Vector3d &r_ref, const Vector3d &r_cart) -> Vector3d {
          return CartesianToEastNortUp(r_ref, r_cart).cast<double>();
        });
  m.def("cartesian_to_azimuth_elevation_range",
        [](const Vector3d &r_cart_ref, const Vector3d &r_cart) -> Vector3d {
          return CartesianToAzimuthElevationRange(r_cart_ref, r_cart)
              .cast<double>();
        });
  m.def("azimuth_elevation_range_to_cartesian",
        [](const Vector3d &r_aer_ref, const Vector3d &r_aer) -> Vector3d {
          return AzimuthElevationRangeToCartesian(r_aer_ref, r_aer)
              .cast<double>();
        });
}
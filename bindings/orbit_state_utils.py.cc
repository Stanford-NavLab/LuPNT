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

  m.def("cartesian_to_classical",
        [](const Vector6d &cart, double mu) -> Vector6d {
          return CartesianToClassical(cart, mu).cast<double>();
        });
  m.def("cartesian_to_classical",
        [](const CartesianOrbitState &cart, double mu) -> ClassicalOE {
          return CartesianToClassical(cart, mu);
        });

  m.def("classical_to_cartesian",
        [](const Vector6d &coe, double mu) -> Vector6d {
          return ClassicalToCartesian(coe, mu).cast<double>();
        });
  m.def("classical_to_cartesian",
        [](const ClassicalOE &coe, double mu) -> CartesianOrbitState {
          return ClassicalToCartesian(coe, mu);
        });

  m.def("classical_to_quasi_nonsingular", [](const Vector6d &coe) -> Vector6d {
    return ClassicalToQuasiNonsingular(coe).cast<double>();
  });
  m.def("classical_to_quasi_nonsingular",
        [](const ClassicalOE &coe) -> QuasiNonsingularOE {
          return ClassicalToQuasiNonsingular(coe);
        });

  m.def("classical_to_equinoctial", [](const Vector6d &coe) -> Vector6d {
    return ClassicalToEquinoctial(coe).cast<double>();
  });
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

  m.def("equinoctial_to_classical", [](const Vector6d &eqoe) -> Vector6d {
    return EquinoctialToClassical(eqoe).cast<double>();
  });
  m.def("equinoctial_to_classical",
        [](const EquinoctialOE &eqoe, double mu) -> ClassicalOE {
          return EquinoctialToClassical(eqoe, mu);
        });

  m.def("delaunay_to_classical",
        [](const Vector6d &deloe, double mu) -> Vector6d {
          return DelaunayToClassical(deloe, mu).cast<double>();
        });
  m.def("delaunay_to_classical",
        [](const DelaunayOE &deloe, double mu) -> ClassicalOE {
          return DelaunayToClassical(deloe, mu);
        });

  m.def("relative_quasi_nonsingular_to_classical",
        [](const Vector6d &coe, const Vector6d &rel_qnsoe) -> Vector6d {
          return RelativeQuasiNonsingularToClassical(coe, rel_qnsoe)
              .cast<double>();
        });
  m.def("relative_quasi_nonsingular_to_classical",
        [](const ClassicalOE &coe,
           const QuasiNonsingularROE &rel_qnsoe) -> ClassicalOE {
          return RelativeQuasiNonsingularToClassical(coe, rel_qnsoe);
        });

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
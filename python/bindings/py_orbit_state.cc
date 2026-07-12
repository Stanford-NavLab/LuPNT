// lupnt
#include <lupnt/core/constants.h>
#include <lupnt/numerics/math_utils.h>

#include "lupnt/states/state.h"

// pybind11

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

#define DEFINE_GETSET(class, attribute) &class ::Get##attribute, &class## ::Set##attribute

#define DEFINE_GETSET_REAL(class, attribute)                    \
  [](const class& s) -> double { return s.attribute().val(); }, \
      [](class& s, double val) { s.Set_##attribute(val); }

#define DEFINE_GETSET_REALVEC(class, attribute, type)                       \
  [](const class& s) -> type { return s.Get##attribute().cast<double>(); }, \
      [](class& s, type val) { s.Set##attribute(val.cast<real>()); }

#define DEFINE_REPR(class)                                                                       \
  [](const class& s) -> std::string {                                                            \
    std::stringstream ss;                                                                        \
    ss << "<pylupnt." << #class << "(" << s.transpose().format(FMT_COMPACT) << ", " << s.frame() \
       << ">";                                                                                   \
    return ss.str();                                                                             \
  }

void init_orbit_state(py::module& m) {
  py::class_<State>(m, "State", "Generic labeled state vector base class.");
  py::class_<ClassicalOE, State>(m, "ClassicalOE",
                                 "Classical Keplerian orbital elements [a, e, i, Omega, w, M].");
  py::class_<Cart6, State>(m, "Cart6", "Cartesian position+velocity state [r (m); v (m/s)].");
  py::class_<QuasiNonsingularOE, State>(m, "QuasiNonsingularOE",
                                        "Quasi-nonsingular orbital elements [a, u, ex, ey, i, "
                                        "Omega] (singularity-free at e=0, i=0).");
  py::class_<EquinoctialOE, State>(
      m, "EquinoctialOE",
      "Equinoctial orbital elements [a, h, k, p, q, lon] (singularity-free at e=0, i=0).");
  py::class_<SingularROE, State>(
      m, "SingularROE", "Singular relative orbital elements, scaled by semi-major axis [m].");
  py::class_<QuasiNonsingROE, State>(
      m, "QuasiNonsingROE",
      "Quasi-nonsingular relative orbital elements, scaled by semi-major axis [m].");
}

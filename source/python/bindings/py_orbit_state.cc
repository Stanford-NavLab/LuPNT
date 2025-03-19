// lupnt
#include <lupnt/core/constants.h>
#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/orbit_state.h>

// pybind11
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lupnt;

#define DEFINE_GETSET(class, attribute) &class ::Get##attribute, &class## ::Set##attribute

#define DEFINE_GETSET_REAL(class, attribute)                    \
  [](const class &s) -> double { return s.attribute().val(); }, \
      [](class &s, double val) { s.Set_##attribute(val); }

#define DEFINE_GETSET_REALVEC(class, attribute, type)                       \
  [](const class &s) -> type { return s.Get##attribute().cast<double>(); }, \
      [](class &s, type val) { s.Set##attribute(val.cast<real>()); }

#define DEFINE_REPR(class)                                                                   \
  [](const class &s) -> std::string {                                                        \
    std::stringstream ss;                                                                    \
    ss << "<pylupnt." << #class << "(" << s.GetVec().transpose().format(FMT_COMPACT) << ", " \
       << s.GetFrame() << ">";                                                               \
    return ss.str();                                                                         \
  }

void init_orbit_state(py::module &m) {
  // OrbitStateRepres
  py::enum_<OrbitStateRepres>(m, "OrbitStateRepres")
      .value("CARTESIAN", OrbitStateRepres::CARTESIAN)
      .value("CLASSICAL_OE", OrbitStateRepres::CLASSICAL_OE)
      .value("QUASI_NONSINGULAR_OE", OrbitStateRepres::QUASI_NONSINGULAR_OE)
      .value("EQUINOTICAL_OE", OrbitStateRepres::EQUINOCTIAL_OE)
      .value("SINGULAR_ROE", OrbitStateRepres::SINGULAR_ROE)
      .value("QUASINONSINGULAR_ROE", OrbitStateRepres::QUASINONSINGULAR_ROE)
      .value("DELAUNAY_OE", OrbitStateRepres::DELAUNAY_OE)
      .export_values();

  // OrbitState
  py::class_<OrbitState>(m, "OrbitState")
      .def(py::init<const Vec6d &, Frame, OrbitStateRepres, const std::array<const char *, 6> &,
                    const std::array<const char *, 6> &>(),
           py::arg("vector"), py::arg("frame"), py::arg("state_repres"), py::arg("names"),
           py::arg("units"))
      .def_property(
          "vector", [](const OrbitState &s) -> Vec6d { return s.GetVec().cast<double>(); },
          [](OrbitState &s, const Vec6d &vec) { s.SetVec(vec.cast<Real>()); })
      .def_property("frame", &OrbitState::GetFrame, &OrbitState::SetCoordSystem)
      .def_property("state_repres", &OrbitState::GetOrbitStateRepres,
                    &OrbitState::SetOrbitStateRepres)
      .def_property_readonly("size", &OrbitState::GetSize)
      .def_property_readonly("names", &OrbitState::GetNames)
      .def_property_readonly("units", &OrbitState::GetUnits)
      .def("__repr__", [](const OrbitState &s) {
        std::stringstream ss;
        ss << "<pylupnt.OrbitState(" << s.GetVec().transpose() << ", " << s.GetFrame() << ", "
           << s.GetOrbitStateRepres() << ")>";
        return ss.str();
      });

  // ClassicalOE
  py::class_<ClassicalOE, OrbitState>(m, "ClassicalOE")
      .def(py::init<const Vec6d &, const Frame>(), py::arg("[a, e, i, Omega, w, M]"),
           py::arg("frame") = Frame::MOON_CI)
      .def_property("a", DEFINE_GETSET_REAL(ClassicalOE, a))
      .def_property("e", DEFINE_GETSET_REAL(ClassicalOE, e))
      .def_property("i", DEFINE_GETSET_REAL(ClassicalOE, i))
      .def_property("Omega", DEFINE_GETSET_REAL(ClassicalOE, Omega))
      .def_property("w", DEFINE_GETSET_REAL(ClassicalOE, w))
      .def_property("M", DEFINE_GETSET_REAL(ClassicalOE, M))
      .def("__repr__", DEFINE_REPR(ClassicalOE));

  // CartesianOrbitState
  py::class_<CartesianOrbitState, OrbitState>(m, "CartesianOrbitState")
      .def(py::init<const Vec6d &, const Frame>(), py::arg("rv"), py::arg("frame") = Frame::GCRF)
      .def_property(
          "r", [](const CartesianOrbitState &s) -> Vec3d { return s.r().cast<double>(); },
          [](CartesianOrbitState &s, const Vec3d &vec) { s.Set_r(vec.cast<Real>()); })
      .def_property(
          "v", [](const CartesianOrbitState &s) -> Vec3d { return s.v().cast<double>(); },
          [](CartesianOrbitState &s, const Vec3d &vec) { s.Set_v(vec.cast<Real>()); })
      .def("__repr__", DEFINE_REPR(CartesianOrbitState));

  // QuasiNonsingOE
  py::class_<QuasiNonsingOE, OrbitState>(m, "QuasiNonsingOE")
      .def(py::init<const Vec6d &, const Frame>(), py::arg("rv"), py::arg("frame") = Frame::GCRF)
      .def_property("a", DEFINE_GETSET_REAL(QuasiNonsingOE, a))
      .def_property("u", DEFINE_GETSET_REAL(QuasiNonsingOE, u))
      .def_property("ex", DEFINE_GETSET_REAL(QuasiNonsingOE, ex))
      .def_property("ey", DEFINE_GETSET_REAL(QuasiNonsingOE, ey))
      .def_property("i", DEFINE_GETSET_REAL(QuasiNonsingOE, i))
      .def_property("Omega", DEFINE_GETSET_REAL(QuasiNonsingOE, Omega))
      .def("__repr__", DEFINE_REPR(QuasiNonsingOE));

  // EquinoctialOE
  py::class_<EquinoctialOE, OrbitState>(m, "EquinoctialOE")
      .def(py::init<const Vec6d &, const Frame>(), py::arg("rv"), py::arg("frame") = Frame::GCRF)
      .def_property("a", DEFINE_GETSET_REAL(EquinoctialOE, a))
      .def_property("h", DEFINE_GETSET_REAL(EquinoctialOE, h))
      .def_property("k", DEFINE_GETSET_REAL(EquinoctialOE, k))
      .def_property("p", DEFINE_GETSET_REAL(EquinoctialOE, p))
      .def_property("q", DEFINE_GETSET_REAL(EquinoctialOE, q))
      .def_property("lon", DEFINE_GETSET_REAL(EquinoctialOE, lon))
      .def("__repr__", DEFINE_REPR(EquinoctialOE));

  // SingularROE
  py::class_<SingularROE, OrbitState>(m, "SingularROE")
      .def(py::init<const Vec6d &, const Frame>())
      .def_property("ada", DEFINE_GETSET_REAL(SingularROE, ada))
      .def_property("adM", DEFINE_GETSET_REAL(SingularROE, adM))
      .def_property("ade", DEFINE_GETSET_REAL(SingularROE, ade))
      .def_property("adw", DEFINE_GETSET_REAL(SingularROE, adw))
      .def_property("adi", DEFINE_GETSET_REAL(SingularROE, adi))
      .def_property("adOmega", DEFINE_GETSET_REAL(SingularROE, adOmega))
      .def("__repr__", DEFINE_REPR(SingularROE));

  // QuasiNonsingROE
  py::class_<QuasiNonsingROE, OrbitState>(m, "QuasiNonsingROE")
      .def(py::init<const Vec6d &, const Frame>())
      .def_property("ada", DEFINE_GETSET_REAL(QuasiNonsingROE, ada))
      .def_property("adl", DEFINE_GETSET_REAL(QuasiNonsingROE, adl))
      .def_property("adex", DEFINE_GETSET_REAL(QuasiNonsingROE, adex))
      .def_property("adey", DEFINE_GETSET_REAL(QuasiNonsingROE, adey))
      .def_property("adix", DEFINE_GETSET_REAL(QuasiNonsingROE, adix))
      .def_property("adiy", DEFINE_GETSET_REAL(QuasiNonsingROE, adiy))
      .def("__repr__", DEFINE_REPR(QuasiNonsingROE));
}

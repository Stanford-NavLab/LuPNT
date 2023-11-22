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

#define DEFINE_GETSET(class, name) &class ::Get##name, &class ::Set##name

#define DEFINE_GETSET_REAL(class, name) \
  [](const class &s) -> double { return s.name().val(); }, &class ::Set_##name

#define DEFINE_GETSET_REALVEC(class, name, type)                       \
  [](const class &s) -> type { return s.Get##name().cast<double>(); }, \
      &class ::Set##name

#define DEFINE_REPR(class)                                                    \
  [](const class &s) {                                                        \
    std::stringstream ss;                                                     \
    ss << "<pylupnt." << #class << " [" << s.GetVector().transpose() << "]>"; \
    return ss.str();                                                          \
  }

class PyOrbitState : public OrbitState {
 public:
  /* inherit the constructors*/
  PyOrbitState(const Vector6d &vec, const CoordSystem cs,
               const OrbitStateRepres sr)
      : OrbitState(vec, cs, sr) {}

  /* Trampoline (need one for each virtual functions)*/
  void Print(const bool deg) const override {
    PYBIND11_OVERLOAD_PURE(
        void,       /* Return type */
        OrbitState, /* Parent class */
        Print,      /* Name of function in C++ (must match Python name) */
        deg         /* Argument(s) */
    );
  }

  std::shared_ptr<OrbitState> Clone() const override {
    PYBIND11_OVERLOAD_PURE(
        std::shared_ptr<OrbitState>, /* Return type */
        OrbitState,                  /* Parent class */
        Clone, /* Name of function in C++ (must match Python name) */
    );
  }
};

void init_orbit_state(py::module &m) {
  // OrbitStateRepres
  py::enum_<OrbitStateRepres>(m, "OrbitStateRepres")
      .value("CARTESIAN", OrbitStateRepres::CARTESIAN)
      .value("CLASSICAL_OE", OrbitStateRepres::CLASSICAL_OE)
      .value("QUASI_NONSINGULAR_OE", OrbitStateRepres::QUASI_NONSINGULAR_OE)
      .value("NONSINGULAR_OE", OrbitStateRepres::NONSINGULAR_OE)
      .value("EQUINOTICAL_OE", OrbitStateRepres::EQUINOCTIAL_OE)
      .value("SINGULAR_ROE", OrbitStateRepres::SINGULAR_ROE)
      .value("QUASINONSINGULAR_ROE", OrbitStateRepres::QUASINONSINGULAR_ROE)
      .value("DELAUNAY_OE", OrbitStateRepres::DELAUNAY_OE)
      .export_values();

  // OrbitState
  py::class_<OrbitState, PyOrbitState>(m, "OrbitState")
      .def(py::init<const Vector6d &, const CoordSystem,
                    const OrbitStateRepres>())
      .def_property(
          "vector",
          [](const OrbitState &s) -> Vector6d {
            return s.GetVector().cast<double>();
          },
          &OrbitState::SetVector)
      .def_property("coord_sys", &OrbitState::GetCoordSystem,
                    &OrbitState::SetCoordSystem)
      .def_property("state_repres", &OrbitState::GetOrbitStateRepres,
                    &OrbitState::SetOrbitStateRepres)
      .def_property_readonly("size", &OrbitState::GetSize)
      .def("__copy__", &OrbitState::Clone)
      .def("__repr__", [](const OrbitState &s) {
        std::stringstream ss;
        ss << "<pylupnt.OrbitState [" << s.GetVector().transpose() << "]>";
        return ss.str();
      });

  // ClassicalOE
  py::class_<ClassicalOE, OrbitState>(m, "ClassicalOE")
      .def(py::init<const Vector6d &, const CoordSystem>(),
           py::arg("[a, e, i, Omega, w, M]"),
           py::arg("coord_sys") = CoordSystem::NONE)
      .def("print", &ClassicalOE::Print, py::arg("deg") = true)
      .def("clone", &ClassicalOE::Clone)
      .def_property("a", DEFINE_GETSET_REAL(ClassicalOE, a))
      .def_property("e", DEFINE_GETSET_REAL(ClassicalOE, e))
      .def_property("i", DEFINE_GETSET_REAL(ClassicalOE, i))
      .def_property("Omega", DEFINE_GETSET_REAL(ClassicalOE, Omega))
      .def_property("W", DEFINE_GETSET_REAL(ClassicalOE, w))
      .def("__repr__", DEFINE_REPR(ClassicalOE));

  // CartesianOrbitState
  py::class_<CartesianOrbitState, OrbitState>(m, "CartesianOrbitState")
      .def(py::init<const Vector6d &, const CoordSystem>(), py::arg("rv"),
           py::arg("coord_sys") = CoordSystem::NONE)
      .def("print", &CartesianOrbitState::Print, py::arg("deg") = true)
      .def("clone", &CartesianOrbitState::Clone)
      .def_property(
          "r",
          [](const CartesianOrbitState &s) -> Vector3d {
            return s.r().cast<double>();
          },
          &CartesianOrbitState::Set_r)
      .def_property(
          "v",
          [](const CartesianOrbitState &s) -> Vector3d {
            return s.v().cast<double>();
          },
          &CartesianOrbitState::Set_v)
      .def("__repr__", DEFINE_REPR(CartesianOrbitState));

  // QuasiNonsingularOE
  py::class_<QuasiNonsingularOE, OrbitState>(m, "QuasiNonsingularOE")
      .def(py::init<const Vector6d &, const CoordSystem>())
      .def("print", &QuasiNonsingularOE::Print, py::arg("deg") = true)
      .def("clone", &QuasiNonsingularOE::Clone)
      .def_property("a", DEFINE_GETSET_REAL(QuasiNonsingularOE, a))
      .def_property("u", DEFINE_GETSET_REAL(QuasiNonsingularOE, u))
      .def_property("ex", DEFINE_GETSET_REAL(QuasiNonsingularOE, ex))
      .def_property("ey", DEFINE_GETSET_REAL(QuasiNonsingularOE, ey))
      .def_property("i", DEFINE_GETSET_REAL(QuasiNonsingularOE, i))
      .def_property("Omega", DEFINE_GETSET_REAL(QuasiNonsingularOE, Omega))
      .def("__repr__", DEFINE_REPR(QuasiNonsingularOE));

  // NonsingularOE
  py::class_<NonsingularOE, OrbitState>(m, "NonsingularOE")
      .def(py::init<const Vector6d &, const CoordSystem>())
      .def("print", &NonsingularOE::Print, py::arg("deg") = true)
      .def("clone", &NonsingularOE::Clone)
      .def_property("a", DEFINE_GETSET_REAL(NonsingularOE, a))
      .def_property("e1", DEFINE_GETSET_REAL(NonsingularOE, e1))
      .def_property("e2", DEFINE_GETSET_REAL(NonsingularOE, e2))
      .def_property("e3", DEFINE_GETSET_REAL(NonsingularOE, e3))
      .def_property("e4", DEFINE_GETSET_REAL(NonsingularOE, e4))
      .def_property("e5", DEFINE_GETSET_REAL(NonsingularOE, e5))
      .def("__repr__", DEFINE_REPR(NonsingularOE));

  // EquinoctialOE
  py::class_<EquinoctialOE, OrbitState>(m, "EquinoctialOE")
      .def(py::init<const Vector6d &, const CoordSystem>())
      .def("print", &EquinoctialOE::Print, py::arg("deg") = true)
      .def("clone", &EquinoctialOE::Clone)
      .def_property("a", DEFINE_GETSET_REAL(EquinoctialOE, a))
      .def_property("h", DEFINE_GETSET_REAL(EquinoctialOE, h))
      .def_property("k", DEFINE_GETSET_REAL(EquinoctialOE, k))
      .def_property("p", DEFINE_GETSET_REAL(EquinoctialOE, p))
      .def_property("q", DEFINE_GETSET_REAL(EquinoctialOE, q))
      .def_property("lon", DEFINE_GETSET_REAL(EquinoctialOE, lon))
      .def("__repr__", DEFINE_REPR(EquinoctialOE));

  // SingularROE
  py::class_<SingularROE, OrbitState>(m, "SingularROE")
      .def(py::init<const Vector6d &, const CoordSystem>())
      .def("print", &SingularROE::Print, py::arg("deg") = true)
      .def("clone", &SingularROE::Clone)
      .def_property("da", DEFINE_GETSET_REAL(SingularROE, da))
      .def_property("dM", DEFINE_GETSET_REAL(SingularROE, dM))
      .def_property("de", DEFINE_GETSET_REAL(SingularROE, de))
      .def_property("dw", DEFINE_GETSET_REAL(SingularROE, dw))
      .def_property("di", DEFINE_GETSET_REAL(SingularROE, di))
      .def_property("dOmega", DEFINE_GETSET_REAL(SingularROE, dOmega))
      .def("__repr__", DEFINE_REPR(SingularROE));

  // QuasiNonsingularROE
  py::class_<QuasiNonsingularROE, OrbitState>(m, "QuasiNonsingularROE")
      .def(py::init<const Vector6d &, const CoordSystem>())
      .def("print", &QuasiNonsingularROE::Print, py::arg("deg") = true)
      .def("clone", &QuasiNonsingularROE::Clone)
      .def_property("da", DEFINE_GETSET_REAL(QuasiNonsingularROE, da))
      .def_property("dl", DEFINE_GETSET_REAL(QuasiNonsingularROE, dl))
      .def_property("dex", DEFINE_GETSET_REAL(QuasiNonsingularROE, dex))
      .def_property("dey", DEFINE_GETSET_REAL(QuasiNonsingularROE, dey))
      .def_property("dix", DEFINE_GETSET_REAL(QuasiNonsingularROE, dix))
      .def_property("diy", DEFINE_GETSET_REAL(QuasiNonsingularROE, diy))
      .def("__repr__", DEFINE_REPR(QuasiNonsingularROE));

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

  // coe <-> roe
  m.def(
      "roe_to_coe",
      [](const ClassicalOE &coe_chief, const QuasiNonsingularROE &roe) {
        return RoeToCoe(coe_chief, roe);
      },
      py::arg("coe_chief"), py::arg("roe"));
  m.def(
      "roe_to_coe",
      [](const Vector6d &coe_chief, const Vector6d &roe) -> Vector6d {
        return RoeToCoe(coe_chief, roe).cast<double>();
      },
      py::arg("coe_chief"), py::arg("roe"));

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
    return EccentricAnomToTrueAnom(E, e).val();
  });
  m.def("eccentric_to_mean", [](double E, double e) -> double {
    return EccentricAnomToMeanAnom(E, e).val();
  });
  m.def("mean_to_eccentric", [](double M, double e) -> double {
    return MeanAnomToEccentricAnom(M, e).val();
  });
  m.def("mean_to_true", [](double M, double e) -> double {
    return MeanAnomToTrueAnom(M, e).val();
  });
  m.def("true_to_eccentric", [](double nu, double e) -> double {
    return TrueAnomToEccentricAnom(nu, e).val();
  });
  m.def("true_to_mean", [](double f, double e) -> double {
    return TrueAnomToMeanAnom(f, e).val();
  });
}
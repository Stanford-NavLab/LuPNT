#include <lupnt/core/Constants.h>
#include <lupnt/physics/OrbitState.h>
#include <lupnt/physics/OrbitStateUtils.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include <Eigen/Dense>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

#include "eigen.hxx"

namespace py = pybind11;
using namespace LPT;

class PyOrbitState : public OrbitState {
 public:
  /* inherit the constructors*/
  PyOrbitState(const ad::Vector6real &vec, const CoordSystem cs,
          const OrbitStateRepres sr)
      : OrbitState(vec, cs, sr) {}
  PyOrbitState(const ad::real x, const ad::real y, const ad::real z,
          const ad::real vx, const ad::real vy, const ad::real vz,
          const CoordSystem cs, const OrbitStateRepres sr)
      : OrbitState(x, y, z, vx, vy, vz, cs, sr) {}

  /* Trampoline (need one for each virtual functions)*/
  void Print(const bool deg) const override {
    PYBIND11_OVERLOAD_PURE(
        void,  /* Return type */
        OrbitState, /* Parent class */
        Print, /* Name of function in C++ (must match Python name) */
        deg    /* Argument(s) */
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

void init_state(py::module &m) {
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

  py::class_<OrbitState, PyOrbitState>(m, "OrbitState")
      .def(py::init<const ad::Vector6real &, const CoordSystem,
                    const OrbitStateRepres>())
      .def(py::init<const ad::real, const ad::real, const ad::real,
                    const ad::real, const ad::real, const ad::real,
                    const CoordSystem, const OrbitStateRepres>())

      .def("get_vector", &OrbitState::GetVector)
      .def("set_vector", &OrbitState::SetVector)
      .def("get_coord_system", &OrbitState::GetCoordSystem)
      .def("set_coord_system", &OrbitState::SetCoordSystem)
      .def("get_state_repres", &OrbitState::GetOrbitStateRepres)
      .def("set_state_repres", &OrbitState::SetOrbitStateRepres)
      .def("set_value", &OrbitState::SetValue)

      .def("print", &OrbitState::Print, py::arg("deg") = true)
      .def("clone", &OrbitState::Clone)
      .def("__repr__", [](const OrbitState &s) {
        return "<pylupnt.OrbitState " + std::to_string(s.GetVector().size()) + ">";
      });

  py::class_<CartesianOrbitState, OrbitState>(m, "CartesianOrbitState")
      .def(py::init<const ad::Vector6real &, const CoordSystem>(),
           py::arg("rv"), py::arg("cs") = CoordSystem::NONE)
      .def("print", &CartesianOrbitState::Print, py::arg("deg") = true)
      .def("clone", &CartesianOrbitState::Clone)
      .def("r", &CartesianOrbitState::r)
      .def("v", &CartesianOrbitState::v)
      .def("set_r", &CartesianOrbitState::Set_r)
      .def("set_v", &CartesianOrbitState::Set_v)
      .def("__repr__",
           [](const CartesianOrbitState &s) { return "<pylupnt.CartesianOrbitState>"; });

  py::class_<ClassicalOE, OrbitState>(m, "ClassicalOE")
      .def(py::init<const ad::Vector6real &, const CoordSystem>(),
           py::arg("oe"), py::arg("cs") = CoordSystem::NONE)
      .def(py::init<const ad::real, const ad::real, const ad::real,
                    const ad::real, const ad::real, const ad::real,
                    const CoordSystem>(),
           py::arg("a"), py::arg("e"), py::arg("i"), py::arg("Omega"),
           py::arg("omega"), py::arg("M"), py::arg("cs") = CoordSystem::NONE)

      .def("print", &ClassicalOE::Print, py::arg("deg") = true)
      .def("clone", &ClassicalOE::Clone)
      .def("a", &ClassicalOE::a)
      .def("e", &ClassicalOE::e)
      .def("i", &ClassicalOE::i)
      .def("Omega", &ClassicalOE::Omega)
      .def("w", &ClassicalOE::w)
      .def("M", &ClassicalOE::M)
      .def("__repr__",
           [](const ClassicalOE &s) { return "<pylupnt.ClassicalOE>"; });

  py::class_<QuasiNonsingularOE, OrbitState>(m, "QuasiNonsingularOE")
      .def(py::init<const ad::Vector6real &, const CoordSystem>())
      .def(py::init<const ad::real, const ad::real, const ad::real,
                    const ad::real, const ad::real, const ad::real,
                    const CoordSystem>(),
           py::arg("a"), py::arg("u"), py::arg("ex"), py::arg("ey"),
           py::arg("i"), py::arg("Omega"), py::arg("cs") = CoordSystem::NONE)
      .def("print", &QuasiNonsingularOE::Print, py::arg("deg") = true)
      .def("clone", &QuasiNonsingularOE::Clone)
      .def("a", &QuasiNonsingularOE::a)
      .def("u", &QuasiNonsingularOE::u)
      .def("ex", &QuasiNonsingularOE::ex)
      .def("ey", &QuasiNonsingularOE::ey)
      .def("i", &QuasiNonsingularOE::i)
      .def("Omega", &QuasiNonsingularOE::Omega);

  py::class_<NonsingularOE, OrbitState>(m, "NonsingularOE")
      .def(py::init<const ad::Vector6real &, const CoordSystem>())
      .def(py::init<const ad::real, const ad::real, const ad::real,
                    const ad::real, const ad::real, const ad::real,
                    const CoordSystem>(),
           py::arg("a"), py::arg("e1"), py::arg("e2"), py::arg("e3"),
           py::arg("e4"), py::arg("e5"), py::arg("cs") = CoordSystem::NONE)
      .def("print", &NonsingularOE::Print, py::arg("deg") = true)
      .def("clone", &NonsingularOE::Clone)
      .def("a", &NonsingularOE::a)
      .def("e1", &NonsingularOE::e1)
      .def("e2", &NonsingularOE::e2)
      .def("e3", &NonsingularOE::e3)
      .def("e4", &NonsingularOE::e4)
      .def("e5", &NonsingularOE::e5);

  py::class_<EquinoctialOE, OrbitState>(m, "EquinoctialOE")
      .def(py::init<const ad::Vector6real &, const CoordSystem>())
      .def(py::init<const ad::real, const ad::real, const ad::real,
                    const ad::real, const ad::real, const ad::real,
                    const CoordSystem>(),
           py::arg("a"), py::arg("h"), py::arg("k"), py::arg("p"), py::arg("q"),
           py::arg("lon"), py::arg("cs") = CoordSystem::NONE)
      .def("print", &EquinoctialOE::Print, py::arg("deg") = true)
      .def("clone", &EquinoctialOE::Clone)
      .def("a", &EquinoctialOE::a)
      .def("h", &EquinoctialOE::h)
      .def("k", &EquinoctialOE::k)
      .def("p", &EquinoctialOE::p)
      .def("q", &EquinoctialOE::q)
      .def("lon", &EquinoctialOE::lon);

  // py::class_<SingularROE, OrbitState>(m, "SingularROE")
  //     .def(py::init<const ad::Vector6real &, const CoordSystem>())
  //     .def(py::init<const ad::real, const ad::real, const ad::real, const
  //     ad::real,
  //          const ad::real, const ad::real, const CoordSystem>(),
  //          py::arg("da"), py::arg("dM"), py::arg("de"), py::arg("dw"),
  //          py::arg("di"), py::arg("dOmega"), py::arg("cs") =
  //          CoordSystem::NONE)
  //     .def("print", &SingularROE::Print, py::arg("deg") = true)
  //     .def("clone", &SingularROE::Clone)
  //     .def("da", &SingularROE::da)
  //     .def("dM", &SingularROE::dM)
  //     .def("de", &SingularROE::de)
  //     .def("dw", &SingularROE::dw)
  //     .def("di", &SingularROE::di)
  //     .def("dOmega", &SingularROE::dOmega);

  py::class_<QuasiNonsingularROE, OrbitState>(m, "QuasiNonsingularROE")
      .def(py::init<const ad::Vector6real &, const CoordSystem>())
      .def(py::init<const ad::real, const ad::real, const ad::real,
                    const ad::real, const ad::real, const ad::real,
                    const CoordSystem>(),
           py::arg("da"), py::arg("dl"), py::arg("dex"), py::arg("dey"),
           py::arg("dix"), py::arg("diy"), py::arg("cs") = CoordSystem::NONE)
      .def("print", &QuasiNonsingularROE::Print, py::arg("deg") = true)
      .def("clone", &QuasiNonsingularROE::Clone)
      .def("da", &QuasiNonsingularROE::da)
      .def("dl", &QuasiNonsingularROE::dl)
      .def("dex", &QuasiNonsingularROE::dex)
      .def("dey", &QuasiNonsingularROE::dey)
      .def("dix", &QuasiNonsingularROE::dix)
      .def("diy", &QuasiNonsingularROE::diy);

  py::class_<TLE>(m, "TLE")
      .def("FromLines", &TLE::FromLines)
      .def("FromFile", &TLE::FromFile);

  // functions
  m.def("convert_state_representation", &ConvertOrbitStateRepresentation,
        py::arg("state"), py::arg("to"), py::arg("mu"));

  m.def("convert_state_coord_system", &ConvertOrbitStateCoordSystem,
        py::arg("state"), py::arg("epoch"), py::arg("to"));

  // OrbitState Conversion Functions
  m.def(
      "coe_to_cart",
      [](const ClassicalOE &coe, double mu) -> CartesianOrbitState {
        return CoeToCart(coe, mu);
      },
      py::arg("coe"), py::arg("mu"));

  m.def(
      "coe_to_cart",
      [](const ad::Vector6real &coeVec, double mu) -> ad::Vector6real {
        return CoeToCart(coeVec, mu);
      },
      py::arg("coe"), py::arg("mu"));
  m.def(
      "coe_to_cart",
      [](const Eigen::Vector6d &coeVec, double mu) -> Eigen::Vector6d {
        return CoeToCart(coeVec, mu);
      },
      py::arg("coe"), py::arg("mu"));

  //   m.def("cart_to_coe",
  //         py::overload_cast<const CartesianOrbitState, const double>(&CartToCoe),
  //         py::arg("cart"), py::arg("mu"));
  //   //   m.def("cart_to_coe",
  //         py::overload_cast<const ad::Vector6real &, const
  //         double>(&CartToCoe), py::arg("cart"), py::arg("mu"));

  m.def("roe_to_coe",
        py::overload_cast<const ClassicalOE, const QuasiNonsingularROE>(
            &RoeToCoe),
        py::arg("coe_chief"), py::arg("roe"));
  m.def("roe_to_coe",
        py::overload_cast<const ad::Vector6real &, const ad::Vector6real &>(
            &RoeToCoe),
        py::arg("coe_chief"), py::arg("roe"));

  m.def("inertial_to_rtn",
        py::overload_cast<const CartesianOrbitState &, const CartesianOrbitState &>(
            &InertialToRtn),
        py::arg("cart_orig"), py::arg("cart"));
  m.def("inertial_to_rtn",
        py::overload_cast<const ad::Vector6real &, const ad::Vector6real &>(
            &InertialToRtn),
        py::arg("cart_orig"), py::arg("cart"));

  m.def("coe_to_rtn",
        py::overload_cast<const ClassicalOE &, const ClassicalOE &, double>(
            &CoeToRtn),
        py::arg("coe_chief"), py::arg("coe_deputy"), py::arg("mu"));
  m.def("coe_to_rtn",
        py::overload_cast<const ad::Vector6real &, const ad::Vector6real &,
                          double>(&CoeToRtn),
        py::arg("coe_chief"), py::arg("coe_deputy"), py::arg("mu"));

  // Todo: add overloads for other conversions

  // Anomaly Conversions
  m.def("eccentric_to_true", &EccentricAnomToTrueAnom, py::arg("E"),
        py::arg("e"));
  m.def("eccentric_to_mean", &EccentricAnomToMeanAnom, py::arg("E"),
        py::arg("e"));
  m.def("mean_to_eccentric", &MeanAnomToEccentricAnom, py::arg("M"),
        py::arg("e"));
  m.def("mean_to_true", &MeanAnomToTrueAnom, py::arg("M"), py::arg("e"));
  m.def("true_to_eccentric", &TrueAnomToEccentricAnom, py::arg("f"),
        py::arg("e"));
  m.def("true_to_mean", &TrueAnomToMeanAnom, py::arg("f"), py::arg("e"));
}
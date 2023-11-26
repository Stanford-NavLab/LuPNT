#include <lupnt/core/constants.h>
#include <lupnt/dynamics/dynamics.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lupnt;

void init_dynamics(py::module &m) {
  // KeplerianDynamics
  py::class_<KeplerianDynamics>(m, "KeplerianDynamics")
      .def(py::init<const double>())
      .def(
          "propagate",
          py::overload_cast<ClassicalOE &, real>(&KeplerianDynamics::Propagate),
          py::arg("state"), py::arg("dt"))
      .def("propagate",
           py::overload_cast<QuasiNonsingularOE &, real>(
               &KeplerianDynamics::Propagate),
           py::arg("state"), py::arg("dt"))
      .def("propagate",
           py::overload_cast<EquinoctialOE &, real>(
               &KeplerianDynamics::Propagate),
           py::arg("state"), py::arg("dt"))
      .def(
          "propagate_with_stm",
          [](KeplerianDynamics &dyn, ClassicalOE &state, real dt) {
            Matrix6d stm;
            dyn.PropagateWithStm(state, dt, stm);
            return stm;
          },
          py::arg("state"), py::arg("dt"), py::return_value_policy::move)
      .def(
          "propagate_with_stm",
          [](KeplerianDynamics &dyn, QuasiNonsingularOE &state, real dt) {
            Matrix6d stm;
            dyn.PropagateWithStm(state, dt, stm);
            return stm;
          },
          py::arg("state"), py::arg("dt"), py::return_value_policy::move)
      .def(
          "propagate_with_stm",
          [](KeplerianDynamics &dyn, EquinoctialOE &state, real dt) {
            Matrix6d stm;
            dyn.PropagateWithStm(state, dt, stm);
            return stm;
          },
          py::arg("state"), py::arg("dt"), py::return_value_policy::move)
      .def("__repr__", [](const KeplerianDynamics &dyn) {
        return "<pylupnt.KeplerianDynamics>";
      });

  // NumericalDynamics
  py::class_<NumericalDynamics>(m, "NumericalDynamics")
      .def("propagate",
           py::overload_cast<OrbitState &, real, real, real>(
               &NumericalDynamics::Propagate),
           py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("dt"))
      .def("propagate",
           py::overload_cast<Vector6 &, real, real, real>(
               &NumericalDynamics::Propagate),
           py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("dt"))
      .def(
          "propagate_with_stm",
          [](NumericalDynamics &dyn, CartesianOrbitState &state, real t0,
             real tf, real dt) {
            Matrix6d stm;
            dyn.PropagateWithStm(state, t0, tf, dt, stm);
            return stm;
          },
          py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("dt"),
          py::return_value_policy::move)
      .def(
          "propagate_with_stm",
          [](NumericalDynamics &dyn, Vector6 &state, real t0, real tf,
             real dt) {
            Matrix6d stm;
            dyn.PropagateWithStm(state, t0, tf, dt, stm);
            return stm;
          },
          py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("dt"),
          py::return_value_policy::move);

  // CartesianTwoBodyDynamics
  py::class_<CartesianTwoBodyDynamics, NumericalDynamics>(
      m, "CartesianTwoBodyDynamics")
      .def(py::init<double, std::string>(), py::arg("mu"),
           py::arg("integrator") = "RK4")
      .def("__repr__", [](const CartesianTwoBodyDynamics &dyn) {
        return "<pylupnt.CartesianTwoBodyDynamics>";
      });
}
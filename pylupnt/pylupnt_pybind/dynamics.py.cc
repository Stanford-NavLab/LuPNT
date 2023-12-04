#include <lupnt/core/constants.h>
#include <lupnt/dynamics/dynamics.h>
#include <lupnt/numerics/math_utils.h>
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
          [](KeplerianDynamics &dyn, ClassicalOE &state, double dt) -> void {
            dyn.Propagate(state, dt);
          },
          py::arg("state"), py::arg("dt"))
      .def(
          "propagate",
          [](KeplerianDynamics &dyn, QuasiNonsingularOE &state,
             double dt) -> void { dyn.Propagate(state, dt); },
          py::arg("state"), py::arg("dt"))
      .def(
          "propagate",
          [](KeplerianDynamics &dyn, EquinoctialOE &state, double dt) -> void {
            dyn.Propagate(state, dt);
          },
          py::arg("state"), py::arg("dt"))
      .def(
          "propagate_with_stm",
          [](KeplerianDynamics &dyn, ClassicalOE &state, double dt) {
            Matrix6d stm;
            dyn.PropagateWithStm(state, dt, stm);
            return stm;
          },
          py::arg("state"), py::arg("dt"), py::return_value_policy::move)
      .def(
          "propagate_with_stm",
          [](KeplerianDynamics &dyn, QuasiNonsingularOE &state, double dt) {
            Matrix6d stm;
            dyn.PropagateWithStm(state, dt, stm);
            return stm;
          },
          py::arg("state"), py::arg("dt"), py::return_value_policy::move)
      .def(
          "propagate_with_stm",
          [](KeplerianDynamics &dyn, EquinoctialOE &state, double dt) {
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
      .def(
          "propagate",
          [](NumericalDynamics &dyn, OrbitState &state, double t0, double tf,
             double dt) -> void { dyn.Propagate(state, t0, tf, dt); },
          py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("dt"))
      .def(
          "propagate",
          [](NumericalDynamics &dyn, Vector6d &state, double t0, double tf,
             double dt) -> Vector6d {
            Vector6 x = state.cast<real>();
            dyn.Propagate(x, t0, tf, dt);
            return x.cast<double>();
          },
          py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("dt"))
      .def(
          "propagate_with_stm",
          [](NumericalDynamics &dyn, CartesianOrbitState &state, double t0,
             double tf, double dt) {
            Matrix6d stm;
            dyn.PropagateWithStm(state, t0, tf, dt, stm);
            return stm;
          },
          py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("dt"),
          py::return_value_policy::move)
      .def(
          "propagate_with_stm",
          [](NumericalDynamics &dyn, Vector6 &state, double t0, double tf,
             double dt) -> std::tuple<Vector6d, Matrix6d> {
            Matrix6d stm;
            dyn.PropagateWithStm(state, t0, tf, dt, stm);
            return std::make_tuple(state.cast<double>(), stm);
          },
          py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("dt"));

  // CartesianTwoBodyDynamics
  py::class_<CartesianTwoBodyDynamics, NumericalDynamics>(
      m, "CartesianTwoBodyDynamics")
      .def(py::init<double, std::string>(), py::arg("mu"),
           py::arg("integrator") = "RK4")
      .def("__repr__", [](const CartesianTwoBodyDynamics &dyn) {
        return "<pylupnt.CartesianTwoBodyDynamics>";
      });
}
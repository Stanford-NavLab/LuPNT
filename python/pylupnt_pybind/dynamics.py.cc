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

  // NumericalOrbitDynamics
  py::class_<NumericalOrbitDynamics>(m, "NumericalOrbitDynamics")
      .def(
          "propagate",
          [](NumericalOrbitDynamics &dyn, OrbitState &state, double t0,
             double tf,
             double dt) -> void { dyn.Propagate(state, t0, tf, dt); },
          py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("dt"))
      .def(
          "propagate",
          [](NumericalOrbitDynamics &dyn, Vector6d &x, double t0, double tf,
             double dt) -> Vector6d {
            Vector6 x_real = x.cast<real>();
            dyn.Propagate(x_real, t0, tf, dt);
            return x.cast<double>();
          },
          py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("dt"))
      .def(
          "propagate",
          [](NumericalOrbitDynamics &dyn, Vector6d &x, double t0, VectorXd &tfs,
             bool progress) -> MatrixXd {
            Vector6 x_real = x.cast<real>();
            VectorX tfs_real = tfs.cast<real>();
            return dyn.Propagate(x_real, t0, tfs_real, progress).cast<double>();
          },
          py::arg("state"), py::arg("t0"), py::arg("tfs"),
          py::arg("progress") = false)
      .def(
          "propagate_with_stm",
          [](NumericalOrbitDynamics &dyn, CartesianOrbitState &state, double t0,
             double tf, double dt) {
            Matrix6d stm;
            dyn.PropagateWithStm(state, t0, tf, dt, stm);
            return stm;
          },
          py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("dt"),
          py::return_value_policy::move)
      .def(
          "propagate_with_stm",
          [](NumericalOrbitDynamics &dyn, Vector6d &state, double t0, double tf,
             double dt) -> std::tuple<Vector6d, Matrix6d> {
            Matrix6d stm;
            Vector6 state_real = state.cast<real>();
            dyn.PropagateWithStm(state_real, t0, tf, dt, stm);
            return std::make_tuple(state.cast<double>(), stm);
          },
          py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("dt"));

  // CartesianTwoBodyDynamics
  py::class_<CartesianTwoBodyDynamics, NumericalOrbitDynamics>(
      m, "CartesianTwoBodyDynamics")
      .def(py::init<double, std::string>(), py::arg("mu"),
           py::arg("integrator") = "RK4")
      .def("__repr__", [](const CartesianTwoBodyDynamics &dyn) {
        return "<pylupnt.CartesianTwoBodyDynamics>";
      });

  // NBodyDynamics
  py::class_<NBodyDynamics, NumericalOrbitDynamics>(m, "NBodyDynamics")
      .def(py::init<>())
      .def("set_primary_body", &NBodyDynamics::SetPrimaryBody, py::arg("body"))
      .def("add_body", &NBodyDynamics::AddBody, py::arg("body"))
      .def(
          "set_time_step",
          [](NBodyDynamics &dyn, double dt) { dyn.SetTimeStep(dt); },
          py::arg("dt"))
      .def("__repr__",
           [](const NBodyDynamics &dyn) { return "<pylupnt.NBodyDynamics>"; });

  // Body
  py::class_<Body>(m, "Body")
      .def(py::init<>())
      .def_static("Moon", &Body::Moon, py::arg("n_max") = 0,
                  py::arg("m_max") = 0)
      .def_static("Earth", &Body::Earth, py::arg("n_max") = 0,
                  py::arg("m_max") = 0)
      .def("__repr__", [](const Body &body) { return "<pylupnt.Body>"; });
}
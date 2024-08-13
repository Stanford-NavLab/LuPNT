#include <lupnt/core/constants.h>
#include <lupnt/dynamics/dynamics.h>
#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/body.h>
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
          [](KeplerianDynamics &dyn, QuasiNonsingOE &state, double dt) -> void {
            dyn.Propagate(state, dt);
          },
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
            Mat6d stm;
            dyn.PropagateWithStm(state, dt, stm);
            return stm;
          },
          py::arg("state"), py::arg("dt"), py::return_value_policy::move)
      .def(
          "propagate_with_stm",
          [](KeplerianDynamics &dyn, QuasiNonsingOE &state, double dt) {
            Mat6d stm;
            dyn.PropagateWithStm(state, dt, stm);
            return stm;
          },
          py::arg("state"), py::arg("dt"), py::return_value_policy::move)
      .def(
          "propagate_with_stm",
          [](KeplerianDynamics &dyn, EquinoctialOE &state, double dt) {
            Mat6d stm;
            dyn.PropagateWithStm(state, dt, stm);
            return stm;
          },
          py::arg("state"), py::arg("dt"), py::return_value_policy::move)
      .def("__repr__", [](const KeplerianDynamics &dyn) { return "<pylupnt.KeplerianDynamics>"; });

  // NumericalOrbitDynamics
  py::class_<NumericalOrbitDynamics>(m, "NumericalOrbitDynamics")
      .def(
          "propagate",
          [](NumericalOrbitDynamics &dyn, OrbitState &state, double t0, double tf,
             double dt) -> void { dyn.Propagate(state, t0, tf, dt); },
          py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("dt") = 0.0)
      .def(
          "propagate",
          [](NumericalOrbitDynamics &dyn, Vec6d &x, double t0, double tf, double dt) -> Vec6d {
            Vec6 x_real = x.cast<Real>();
            dyn.Propagate(x_real, t0, tf, dt);
            return x_real.cast<double>();
          },
          py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("dt") = 0.0)
      .def(
          "propagate",
          [](NumericalOrbitDynamics &dyn, Vec6d &x, double t0, VecXd &tfs, double dt,
             bool progress) -> MatXd {
            Vec6 x_real = x.cast<Real>();
            VecX tfs_real = tfs.cast<Real>();
            return dyn.Propagate(x_real, t0, tfs_real, progress).cast<double>();
          },
          py::arg("state"), py::arg("t0"), py::arg("tfs"), py::arg("dt") = 0.0,
          py::arg("progress") = false)
      .def(
          "propagate_with_stm",
          [](NumericalOrbitDynamics &dyn, CartesianOrbitState &state, double t0, double tf,
             double dt) {
            Mat6d stm;
            dyn.PropagateWithStm(state, t0, tf, dt, stm);
            return stm;
          },
          py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("dt"),
          py::return_value_policy::move)
      .def(
          "propagate_with_stm",
          [](NumericalOrbitDynamics &dyn, Vec6d &state, double t0, double tf,
             double dt) -> std::tuple<Vec6d, Mat6d> {
            Mat6d stm;
            Vec6 state_real = state.cast<Real>();
            dyn.PropagateWithStm(state_real, t0, tf, dt, stm);
            return std::make_tuple(state_real.cast<double>(), stm);
          },
          py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("dt"))
      .def("set_time_step", [](NumericalOrbitDynamics &dyn, double dt) { dyn.SetTimeStep(dt); });
  // .def(
  //     "propagate_with_stm",
  //     [](NumericalOrbitDynamics &dyn, RowVec6d &state, double t0,
  //        double tf, double dt) -> std::tuple<Vec6d, Mat6d> {
  //       Mat6d stm;
  //       Vec6 state_real = state.cast<real>();
  //       dyn.PropagateWithStm(state_real, t0, tf, dt, stm);
  //       return std::make_tuple(state_real.cast<double>(), stm);
  //     },
  //     py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("dt"));

  // CartesianTwoBodyDynamics
  py::class_<CartesianTwoBodyDynamics, NumericalOrbitDynamics>(m, "CartesianTwoBodyDynamics")
      .def(py::init<double, std::string>(), py::arg("GM"), py::arg("integrator") = "RK4")
      .def("__repr__", [](const CartesianTwoBodyDynamics &dyn) {
        return "<pylupnt.CartesianTwoBodyDynamics>";
      });

  // NBodyDynamics
  py::class_<NBodyDynamics, NumericalOrbitDynamics>(m, "NBodyDynamics")
      .def(py::init<>())
      .def("set_primary_body", &NBodyDynamics::SetPrimaryBody, py::arg("body"))
      .def("add_body", &NBodyDynamics::AddBody, py::arg("body"))
      .def("__repr__", [](const NBodyDynamics &dyn) { return "<pylupnt.NBodyDynamics>"; });

  // Body
  py::class_<Body>(m, "Body")
      .def(py::init<>())
      .def_static("Moon", &Body::Moon, py::arg("n_max") = 0, py::arg("m_max") = 0,
                  py::arg("gravity_file") = "EGM96.cof")
      .def_static("Earth", &Body::Earth, py::arg("n_max") = 0, py::arg("m_max") = 0,
                  py::arg("gravity_file") = "grgm900c.cof")
      .def_static("Mars", &Body::Mars, py::arg("n_max") = 0, py::arg("m_max") = 0,
                  py::arg("gravity_file") = "GMM1.cof")
      .def_static("Venus", &Body::Venus, py::arg("n_max") = 0, py::arg("m_max") = 0,
                  py::arg("gravity_file") = "MGN75HSAAP.cof")
      .def_static("Sun", &Body::Sun)
      .def("__repr__", [](const Body &body) { return "<pylupnt.Body>"; });
}

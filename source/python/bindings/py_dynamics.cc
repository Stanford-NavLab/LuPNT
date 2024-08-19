#include <lupnt/core/constants.h>
#include <lupnt/dynamics/dynamics.h>
#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/body.h>

// pybind11
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lupnt;

class PyIDynamics : public IDynamics {
public:
  using IDynamics::IDynamics;
  // Ptr<IState>
  Ptr<IState> PropagateState(const Ptr<IState> &state, Real t0, Real tf,
                             MatXd *stm = nullptr) override {
    PYBIND11_OVERLOAD_PURE(Ptr<IState>, IDynamics, PropagateState, state, t0, tf, stm);
  }
  // VecX
  VecX Propagate(const VecX &x0, Real t0, Real tf, MatXd *stm = nullptr) override {
    PYBIND11_OVERLOAD_PURE(VecX, IDynamics, Propagate, x0, t0, tf, stm);
  }
};

class PyIOrbitDynamics : public IOrbitDynamics {
public:
  using IOrbitDynamics::IOrbitDynamics;
  // OrbitState
  OrbitState PropagateState(const OrbitState &state, Real t0, Real tf,
                            Mat6d *stm = nullptr) override {
    PYBIND11_OVERLOAD_PURE(OrbitState, IOrbitDynamics, PropagateState, state, t0, tf, stm);
  }
  // Vec6 - Real
  Vec6 Propagate(const Vec6 &x0, Real t0, Real tf, Mat6d *stm = nullptr) override {
    PYBIND11_OVERLOAD_PURE(Vec6, IOrbitDynamics, Propagate, x0, t0, tf, stm);
  }
  // Vec6 - VecX
  MatX6 Propagate(const Vec6 &x0, Real t0, const VecX &tf, bool progress = false) override {
    PYBIND11_OVERLOAD_PURE(MatX6, IOrbitDynamics, Propagate, x0, t0, tf, progress);
  }
};

class PyIAnalyticalOrbitDynamics : public IAnalyticalOrbitDynamics {
public:
  using IAnalyticalOrbitDynamics::IAnalyticalOrbitDynamics;
  // OrbitState
  OrbitState PropagateState(const OrbitState &state, Real t0, Real tf,
                            Mat6d *stm = nullptr) override {
    PYBIND11_OVERLOAD_PURE(OrbitState, IAnalyticalOrbitDynamics, PropagateState, state, t0, tf,
                           stm);
  }
  // Vec6 - Real
  Vec6 Propagate(const Vec6 &x0, Real t0, Real tf, Mat6d *stm = nullptr) override {
    PYBIND11_OVERLOAD_PURE(Vec6, IAnalyticalOrbitDynamics, Propagate, x0, t0, tf, stm);
  }
};

void init_dynamics(py::module &m) {
  // IDynamics
  py::class_<IDynamics, PyIDynamics>(m, "IDynamics")
      .def(py::init<>())
      // Ptr<IState>
      .def(
          "propagate_state",
          [](IDynamics &dyn, const Ptr<IState> &state, double t0, double tf, MatXd *stm) {
            return dyn.PropagateState(state, t0, tf, stm);
          },
          py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("stm") = nullptr)
      // VecX
      .def(
          "propagate",
          [](IDynamics &dyn, const VecXd &x0, double t0, double tf, MatXd *stm) -> VecXd {
            return dyn.Propagate(x0, t0, tf, stm).cast<double>();
          },
          py::arg("x0"), py::arg("t0"), py::arg("tf"), py::arg("stm") = nullptr);

  // IOrbitDynamics
  py::class_<IOrbitDynamics, PyIOrbitDynamics, IDynamics>(m, "IOrbitDynamics")
      .def(py::init<>())
      // OrbitState
      .def(
          "propagate_state",
          [](IOrbitDynamics &dyn, const OrbitState &state, double t0, double tf, Mat6d *stm) {
            return dyn.PropagateState(state, t0, tf, stm);
          },
          py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("stm") = nullptr)
      // Vec6 - Real
      .def(
          "propagate",
          [](IOrbitDynamics &dyn, const Vec6d &x0, double t0, double tf, Mat6d *stm) -> Vec6d {
            return dyn.Propagate(x0, t0, tf, stm).cast<double>();
          },
          py::arg("x0"), py::arg("t0"), py::arg("tf"), py::arg("stm") = nullptr)
      .def(
          "propagate",
          [](IOrbitDynamics &dyn, const RowVec6d &x0, double t0, double tf, Mat6d *stm) -> Vec6d {
            return dyn.Propagate(x0, t0, tf, stm).cast<double>();
          },
          py::arg("x0"), py::arg("t0"), py::arg("tf"), py::arg("stm") = nullptr)
      // MatX6 - Real
      .def(
          "propagate",
          [](IOrbitDynamics &dyn, const Vec6d &x0, double t0, const VecXd &tf, bool progress)
              -> MatX6d { return dyn.Propagate(x0, t0, tf.cast<Real>(), progress).cast<double>(); },
          py::arg("x0"), py::arg("t0"), py::arg("tf"), py::arg("progress") = false)
      // Vec6 - VecX
      .def(
          "propagate",
          [](IOrbitDynamics &dyn, const Vec6d &x0, double t0, VecXd tf, Mat6d *stm) -> MatX6d {
            return dyn.Propagate(x0, t0, tf, stm).cast<double>();
          },
          py::arg("x0"), py::arg("t0"), py::arg("tf"), py::arg("stm") = nullptr);

  // IAnalyticalOrbitDynamics
  py::class_<IAnalyticalOrbitDynamics, PyIAnalyticalOrbitDynamics, IOrbitDynamics>(
      m, "IAnalyticalOrbitDynamics")
      .def(py::init<>());

  // KeplerianDynamics
  py::class_<KeplerianDynamics, IAnalyticalOrbitDynamics>(m, "KeplerianDynamics")
      .def(py::init<double>(), py::arg("GM"));

  // // NumericalOrbitDynamics
  // py::class_<NumericalOrbitDynamics>(m, "NumericalOrbitDynamics")
  //     .def(py::init<ODE, IntegratorType>(), py::arg("odefunc"),
  //     py::arg("integrator")) .def("set_time_step", [](NumericalOrbitDynamics &dyn,
  //     double dt) { dyn.SetTimeStep(dt); }) .def("get_time_step",
  //          [](const NumericalOrbitDynamics &dyn) { return dyn.GetTimeStep().val();
  //          })
  //     .def(
  //         "propagate",
  //         [](NumericalOrbitDynamics &dyn, const Vec6 &x0, double t0, double tf,
  //            MatXd *stm) -> Vec6 { return dyn.Propagate(x0, t0, tf,
  //            stm).cast<double>(); },
  //         py::arg("x0"), py::arg("t0"), py::arg("tf"), py::arg("stm") = nullptr)
  //     .def(
  //         "propagate",
  //         [](NumericalOrbitDynamics &dyn, const Vec6 &x0, double t0, const VecXd
  //         &tf, bool progress)
  //             -> MatX6d { return dyn.Propagate(x0, t0, tf.cast<Real>(),
  //             progress).cast<double>();
  //             },
  //         py::arg("x0"), py::arg("t0"), py::arg("tf"), py::arg("progress") = false)
  //     .I_ORBIT_DYNAMICS_METHODS;

  // // CartesianTwoBodyDynamics
  // py::class_<CartesianTwoBodyDynamics, NumericalOrbitDynamics>(m,
  // "CartesianTwoBodyDynamics")
  //     .def(py::init<double, IntegratorType>(), py::arg("GM"),
  //          py::arg("integ") = default_integrator);

  // // NBodyDynamics
  // py::class_<NBodyDynamics, NumericalOrbitDynamics>(m, "NBodyDynamics")
  //     .def(py::init<>())
  //     .def("set_primary_body", &NBodyDynamics::SetPrimaryBody,
  //     py::arg("body")) .def("add_body", &NBodyDynamics::AddBody,
  //     py::arg("body"));

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

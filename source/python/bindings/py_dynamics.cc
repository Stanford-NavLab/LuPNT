// lupnt
#include <lupnt/core/constants.h>
#include <lupnt/core/definitions.h>
#include <lupnt/dynamics/dynamics.h>
#include <lupnt/numerics/integrator.h>
#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/body.h>

// pybind11
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <atomic>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>

namespace py = pybind11;
using namespace lupnt;

template <class T = IDynamics> class PyIDyn : public T {
public:
  using T::T;
  // Interface
  Ptr<IState> PropagateState(const Ptr<IState> &state, Real t0, Real tf,
                             MatXd *stm = nullptr) override {
    PYBIND11_OVERRIDE_PURE(Ptr<IState>, T, PropagateState, state, t0, tf, stm);
  }
  VecX Propagate(const VecX &x0, Real t0, Real tf, MatXd *stm = nullptr) override {
    PYBIND11_OVERRIDE_PURE(VecX, T, Propagate, x0, t0, tf, stm);
  }
};

template <class T = IOrbitDynamics> class PyIOrbDyn : public PyIDyn<T> {
public:
  using PyIDyn<T>::PyIDyn;
  // Overrides
  // Ptr<IState> PropagateState(const Ptr<IState> &state, Real t0, Real tf,
  //                            MatXd *stm = nullptr) override {
  //   PYBIND11_OVERRIDE(Ptr<IState>, T, PropagateState, state, t0, tf, stm);
  // }
  VecX Propagate(const VecX &x0, Real t0, Real tf, MatXd *stm = nullptr) override {
    PYBIND11_OVERRIDE(VecX, T, Propagate, x0, t0, tf, stm);
  }
  // Interface
  OrbitState PropagateState(const OrbitState &state, Real t0, Real tf,
                            Mat6d *stm = nullptr) override {
    PYBIND11_OVERRIDE_PURE(OrbitState, T, PropagateState, state, t0, tf, stm);
  }
  Vec6 Propagate(const Vec6 &x0, Real t0, Real tf, Mat6d *stm = nullptr) override {
    PYBIND11_OVERRIDE_PURE(Vec6, T, Propagate, x0, t0, tf, stm);
  }
  MatX6 Propagate(const Vec6 &x0, Real t0, const VecX &tf, bool progress = false) override {
    PYBIND11_OVERRIDE_PURE(MatX6, T, Propagate, x0, t0, tf, progress);
  }
  // Implementations
  MatX6 Propagate(const MatX6 &x0, Real t0, Real tf) { return T::Propagate(x0, t0, tf); }
};

template <class T = IAnalyticalOrbitDynamics> class PyIAnOrbDyn : public PyIOrbDyn<T> {
public:
  using PyIOrbDyn<T>::PyIOrbDyn;
  // Overrides
  MatX6 Propagate(const Vec6 &x0, Real t0, const VecX &tf, bool progress = false) override {
    PYBIND11_OVERRIDE(MatX6, T, Propagate, x0, t0, tf, progress);
  }
  // Interface
  OrbitState PropagateState(const OrbitState &state, Real t0, Real tf,
                            Mat6d *stm = nullptr) override {
    PYBIND11_OVERRIDE_PURE(OrbitState, T, PropagateState, state, t0, tf, stm);
  }
  Vec6 Propagate(const Vec6 &x0, Real t0, Real tf, Mat6d *stm = nullptr) override {
    PYBIND11_OVERRIDE_PURE(Vec6, T, Propagate, x0, t0, tf, stm);
  }
};

template <class T> class PyAnOrbDyn : public PyIAnOrbDyn<T> {
public:
  using PyIAnOrbDyn<T>::PyIAnOrbDyn;
  // Overrides
  OrbitState PropagateState(const OrbitState &state, Real t0, Real tf,
                            Mat6d *stm = nullptr) override {
    PYBIND11_OVERRIDE(OrbitState, T, PropagateState, state, t0, tf, stm);
  }
  Vec6 Propagate(const Vec6 &x0, Real t0, Real tf, Mat6d *stm = nullptr) override {
    PYBIND11_OVERRIDE(Vec6, T, Propagate, x0, t0, tf, stm);
  }
};

template <class T = NumericalOrbitDynamics> class PyINumOrbDyn : public PyIOrbDyn<T> {
public:
  using PyIOrbDyn<T>::PyIOrbDyn;
  // Overrides
  Vec6 Propagate(const Vec6 &x0, Real t0, Real tf, Mat6d *stm = nullptr) override {
    PYBIND11_OVERRIDE(Vec6, T, Propagate, x0, t0, tf, stm);
  }
  MatX6 Propagate(const Vec6 &x0, Real t0, const VecX &tf, bool progress = false) override {
    PYBIND11_OVERRIDE(MatX6, T, Propagate, x0, t0, tf, progress);
  }
  // Interface
  OrbitState PropagateState(const OrbitState &state, Real t0, Real tf,
                            Mat6d *stm = nullptr) override {
    PYBIND11_OVERRIDE_PURE(OrbitState, T, PropagateState, state, t0, tf, stm);
  }
  Vec6 ComputeRates(Real t, const Vec6 &x) const override {
    PYBIND11_OVERRIDE_PURE(Vec6, T, ComputeRates, t, x);
  }
};

template <class T> class PyNumOrbDyn : public PyINumOrbDyn<T> {
public:
  using PyINumOrbDyn<T>::PyINumOrbDyn;
  // Overrides
  OrbitState PropagateState(const OrbitState &state, Real t0, Real tf,
                            Mat6d *stm = nullptr) override {
    PYBIND11_OVERRIDE(OrbitState, T, PropagateState, state, t0, tf, stm);
  }
  Vec6 ComputeRates(Real t, const Vec6 &x) const override {
    PYBIND11_OVERRIDE(Vec6, T, ComputeRates, t, x);
  }
};

#define I_ORBIT_DYNAMICS_METHODS(class)                                                           \
  def(                                                                                            \
      "propagate",                                                                                \
      [](class &dyn, const Vec6d &x0, double t0, double tf, bool stm) -> py::object {             \
        Vec6 x0_ = x0.cast<Real>();                                                               \
        Real t0_ = Real(t0), tf_ = Real(tf);                                                      \
        if (stm) {                                                                                \
          Mat6d stm_out;                                                                          \
          Vec6d xf = dyn.Propagate(x0_, t0_, tf_, &stm_out).cast<double>();                       \
          return py::make_tuple(xf.cast<double>(), stm_out.cast<double>());                       \
        } else {                                                                                  \
          Vec6d xf = dyn.Propagate(x0_, t0_, tf_, nullptr).cast<double>();                        \
          return py::cast(xf);                                                                    \
        }                                                                                         \
      },                                                                                          \
      py::arg("x0"), py::arg("t0"), py::arg("tf"), py::arg("stm") = false)                        \
      .def(                                                                                       \
          "propagate",                                                                            \
          [](class &dyn, const RowVec6d &x0, double t0, double tf, bool stm) -> py::object {      \
            Vec6 x0_ = x0.transpose().cast<Real>();                                               \
            Real t0_ = Real(t0), tf_ = Real(tf);                                                  \
            if (stm) {                                                                            \
              Mat6d stm_out;                                                                      \
              RowVec6d xf = dyn.Propagate(x0_, t0_, tf_, &stm_out).transpose().cast<double>();    \
              return py::make_tuple(xf.cast<double>(), stm_out.cast<double>());                   \
            } else {                                                                              \
              RowVec6d xf = dyn.Propagate(x0_, t0_, tf_, nullptr).transpose().cast<double>();     \
              return py::cast(xf);                                                                \
            }                                                                                     \
          },                                                                                      \
          py::arg("x0"), py::arg("t0"), py::arg("tf"), py::arg("stm") = false)                    \
      .def(                                                                                       \
          "propagate",                                                                            \
          [](class &dyn, const MatX6d &x0, double t0, double tf) -> MatX6d {                      \
            MatX6 x0_ = x0.cast<Real>();                                                          \
            Real t0_ = Real(t0), tf_ = Real(tf);                                                  \
            return dyn.Propagate(x0_, t0_, tf_).cast<double>();                                   \
          },                                                                                      \
          py::arg("x0"), py::arg("t0"), py::arg("tf"))                                            \
      .def(                                                                                       \
          "propagate",                                                                            \
          [](class &dyn, const Vec6d &x0, double t0, VecXd tf, bool progress) -> MatX6d {         \
            Vec6 x0_ = x0.cast<Real>();                                                           \
            Real t0_ = Real(t0);                                                                  \
            return dyn.Propagate(x0_, t0_, tf.cast<Real>()).cast<double>();                       \
          },                                                                                      \
          py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>(),             \
          py::arg("x0"), py::arg("t0"), py::arg("tf"), py::arg("progress") = false)               \
      .def(                                                                                       \
          "propagate",                                                                            \
          [](class &dyn, const RowVec6d &x0, double t0, VecXd tf, bool progress) -> MatX6d {      \
            Vec6 x0_ = x0.transpose().cast<Real>();                                               \
            Real t0_ = Real(t0);                                                                  \
            return dyn.Propagate(x0_, t0_, tf.cast<Real>()).cast<double>();                       \
          },                                                                                      \
          py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>(),             \
          py::arg("x0"), py::arg("t0"), py::arg("tf"), py::arg("progress") = false)               \
      .def(                                                                                       \
          "propagate_state",                                                                      \
          [](class &dyn, const OrbitState &state, double t0, double tf, bool stm) -> py::object { \
            OrbitState state_ = state;                                                            \
            Real t0_ = Real(t0), tf_ = Real(tf);                                                  \
            if (stm) {                                                                            \
              Mat6d stm_out;                                                                      \
              OrbitState state_out = dyn.PropagateState(state_, t0_, tf_, &stm_out);              \
              return py::make_tuple(state_out, stm_out);                                          \
            } else {                                                                              \
              OrbitState state_out = dyn.PropagateState(state_, t0_, tf_, nullptr);               \
              return py::cast(state_out);                                                         \
            }                                                                                     \
          },                                                                                      \
          py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("stm") = false);

template <class T> Vec6d def_compute_rates(T &dyn, double t, const Vec6d &x) {
  Real t_ = Real(t);
  Vec6 x_ = x.cast<Real>().array();
  return dyn.ComputeRates(t_, x_).template cast<double>();
}

using ODEWrapper = std::function<VecXd(double, const VecXd &)>;

void init_dynamics(py::module &m) {
  // IntegratorType
  py::enum_<IntegratorType>(m, "IntegratorType")
      .value("RK4", IntegratorType::RK4)
      .value("RK8", IntegratorType::RK8)
      .value("RKF45", IntegratorType::RKF45)
      .export_values();

  // IDynamics
  py::class_<IDynamics, PyIDyn<>>(m, "IDynamics");

  // IOrbitDynamics
  py::class_<IOrbitDynamics, IDynamics, PyIOrbDyn<>>(m, "IOrbitDynamics");

  // IAnalyticalOrbitDynamics
  py::class_<IAnalyticalOrbitDynamics, IOrbitDynamics, PyIAnOrbDyn<>>(m,
                                                                      "IAnalyticalOrbitDynamics");

  // KeplerianDynamics
  py::class_<KeplerianDynamics, IAnalyticalOrbitDynamics, PyIAnOrbDyn<KeplerianDynamics>>(
      m, "KeplerianDynamics")
      .def(py::init<double>(), py::arg("GM"))
      .I_ORBIT_DYNAMICS_METHODS(KeplerianDynamics);

  // NumericalOrbitDynamics
  py::class_<NumericalOrbitDynamics, IOrbitDynamics, PyINumOrbDyn<NumericalOrbitDynamics>>(
      m, "NumericalOrbitDynamics")
      .def("set_ode_function",
           [](NumericalOrbitDynamics &dyn, ODEWrapper odefunc) {
             dyn.SetODEFunction([odefunc](Real t, const VecX &x) -> VecX {
               return odefunc(t.val(), x.cast<double>()).cast<Real>();
             });
           })
      .def("get_time_step", [](NumericalOrbitDynamics &dyn) { return dyn.GetTimeStep().val(); })
      .def("set_time_step", [](NumericalOrbitDynamics &dyn, double dt) { dyn.SetTimeStep(dt); });

  // CartesianTwoBodyDynamics
  py::class_<CartesianTwoBodyDynamics, NumericalOrbitDynamics,
             PyNumOrbDyn<CartesianTwoBodyDynamics>>(m, "CartesianTwoBodyDynamics")
      .def(py::init<double, IntegratorType>(), py::arg("GM"),
           py::arg("integ_type") = default_integrator)
      .def("compute_rates", &def_compute_rates<CartesianTwoBodyDynamics>, py::arg("t"),
           py::arg("x"))
      .I_ORBIT_DYNAMICS_METHODS(CartesianTwoBodyDynamics);

  // J2CartTwoBodyDynamics
  py::class_<J2CartTwoBodyDynamics, NumericalOrbitDynamics, PyNumOrbDyn<J2CartTwoBodyDynamics>>(
      m, "J2CartTwoBodyDynamics")
      .def(py::init<double, double, double, IntegratorType>(), py::arg("GM"), py::arg("J2"),
           py::arg("R_body"), py::arg("integ_type") = default_integrator)
      .def("compute_rates", &def_compute_rates<J2CartTwoBodyDynamics>, py::arg("t"), py::arg("x"))
      .I_ORBIT_DYNAMICS_METHODS(J2CartTwoBodyDynamics);

  // J2KeplerianDynamics
  py::class_<J2KeplerianDynamics, NumericalOrbitDynamics, PyNumOrbDyn<J2KeplerianDynamics>>(
      m, "J2KeplerianDynamics")
      .def(py::init<double, double, double, IntegratorType>(), py::arg("GM"), py::arg("J2"),
           py::arg("R_body"), py::arg("integ_type") = default_integrator)
      .def("compute_rates", &def_compute_rates<J2KeplerianDynamics>, py::arg("t"), py::arg("x"))
      .I_ORBIT_DYNAMICS_METHODS(J2KeplerianDynamics);

  // NBodyDynamics
  py::class_<NBodyDynamics<double>, NumericalOrbitDynamics, PyNumOrbDyn<NBodyDynamics<double>>>(
      m, "NBodyDynamics")
      .def(py::init<IntegratorType>(), py::arg("integ_type") = default_integrator)
      .def("add_body", &NBodyDynamics<double>::AddBody, py::arg("body"))
      .def("get_bodies", &NBodyDynamics<double>::GetBodies)
      .def("set_frame", &NBodyDynamics<double>::SetFrame, py::arg("frame"))
      .I_ORBIT_DYNAMICS_METHODS(NBodyDynamics<double>);

  // Body
  py::class_<Body>(m, "Body")
      .def(py::init<>())
      .def_static("Moon", &Body::Moon, py::arg("n") = 0, py::arg("m") = 0,
                  py::arg("gravity_file") = "grgm900c.cof")
      .def_static("Earth", &Body::Earth, py::arg("n") = 0, py::arg("m") = 0,
                  py::arg("gravity_file") = "EGM96.cof")
      .def_static("Mars", &Body::Mars, py::arg("n") = 0, py::arg("m") = 0,
                  py::arg("gravity_file") = "GMM1.cof")
      .def_static("Venus", &Body::Venus, py::arg("n") = 0, py::arg("m") = 0,
                  py::arg("gravity_file") = "MGN75HSAAP.cof")
      .def_static("Sun", &Body::Sun)
      .def_readonly("id", &Body::id)
      .def_readonly("name", &Body::name)
      .def_readonly("GM", &Body::GM)
      .def_readonly("R", &Body::R)
      .def_readonly("fixed_frame", &Body::fixed_frame)
      .def_readonly("inertial_frame", &Body::inertial_frame)
      .def_readonly("use_gravity_field", &Body::use_gravity_field)
      .def_readonly("gravity_field", &Body::gravity_field);

  // GravityField
  py::class_<GravityField<double>>(m, "GravityField")
      .def(py::init<>())
      .def_readonly("n_max", &GravityField<double>::n_max)
      .def_readonly("m_max", &GravityField<double>::m_max)
      .def_readonly("n", &GravityField<double>::n)
      .def_readonly("m", &GravityField<double>::m)
      .def_readonly("GM", &GravityField<double>::GM)
      .def_readonly("R", &GravityField<double>::R)
      .def_readonly("CS", &GravityField<double>::CS);
}

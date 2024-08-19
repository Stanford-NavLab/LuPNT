// lupnt
#include <lupnt/core/constants.h>
#include <lupnt/core/definitions.h>
#include <lupnt/dynamics/dynamics.h>
#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/body.h>

// pybind11
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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

void init_dynamics(py::module &m) {
  py::class_<IDynamics, PyIDyn<>> i_dynamics(m, "IDynamics");
  py::class_<IOrbitDynamics, IDynamics, PyIOrbDyn<>> i_orbit_dynamics(m, "IOrbitDynamics");
  py::class_<IAnalyticalOrbitDynamics, IOrbitDynamics, PyIAnOrbDyn<>> i_analytical_orbit_dynamics(
      m, "IAnalyticalOrbitDynamics");
  py::class_<KeplerianDynamics, IAnalyticalOrbitDynamics, PyIAnOrbDyn<KeplerianDynamics>>
      keplerian_dynamics(m, "KeplerianDynamics");

  keplerian_dynamics.def(py::init<double>(), py::arg("GM"))
      .def(
          "propagate",
          [](KeplerianDynamics &dyn, const Vec6d &x0, double t0, double tf,
             bool compute_stm) -> py::object {
            Vec6 x0_ = x0.cast<Real>();
            Real t0_ = Real(t0), tf_ = Real(tf);

            if (compute_stm) {
              Mat6d stm_out;  // Initialize stm_out for storing the state transition matrix
              Vec6 xf = dyn.Propagate(x0_, t0_, tf_, &stm_out).cast<double>();
              return py::make_tuple(xf.cast<double>(), stm_out.cast<double>());
            } else {
              Vec6 xf = dyn.Propagate(x0_, t0_, tf_, nullptr).cast<double>();
              return py::cast(xf);
            }
          },
          py::arg("x0"), py::arg("t0"), py::arg("tf"), py::arg("compute_stm") = false);

  // py::class_<KeplerianDynamics>(m, "KeplerianDynamics")
  //     .def(py::init<double>(), py::arg("GM"))
  //     .def(
  //         "propagate",
  //         [](KeplerianDynamics &dyn, const Vec6d &x0, double t0, double tf, Mat6d *stm =
  //         nullptr)
  //             -> Vec6d { return dyn.Propagate(x0.cast<Real>(), t0, tf, stm).cast<double>();
  //             },
  //         py::arg("x0"), py::arg("t0"), py::arg("tf"), py::arg("stm") = nullptr)
  //     .def(
  //         "propagate",
  //         [](KeplerianDynamics &dyn, const RowVec6d &x0, double t0, double tf, Mat6d *stm =
  //         nullptr)
  //             -> RowVec6d { return dyn.Propagate(x0.cast<Real>(), t0, tf,
  //             stm).cast<double>(); },
  // py::arg("x0"), py::arg("t0"), py::arg("tf"), py::arg("stm") = nullptr);
  // .def(
  // "propagate",
  // [](KeplerianDynamics &dyn, const Vec6d &x0, double t0, VecXd tf) -> MatX6 {
  //   return dyn.Propagate(x0.cast<Real>(), t0, tf.cast<Real>()).cast<double>();
  // },
  // py::arg("x0"), py::arg("t0"), py::arg("tf"), py::arg("stm") = nullptr);

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

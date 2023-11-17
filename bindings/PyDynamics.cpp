#include <lupnt/core/Constants.h>
#include <lupnt/dynamics/Dynamics.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include <Eigen/Dense>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

namespace py = pybind11;
using namespace LPT;

void init_dynamics(py::module &m) {
  // KeplerianDynamics
  py::class_<KeplerianDynamics>(m, "KeplerianDynamics")
      .def(py::init<const double>())
      .def("propagate",
           py::overload_cast<ClassicalOE &, ad::real>(
               &KeplerianDynamics::Propagate),
           py::arg("state"), py::arg("dt"))
      .def("propagate",
           py::overload_cast<QuasiNonsingularOE &, ad::real>(
               &KeplerianDynamics::Propagate),
           py::arg("state"), py::arg("dt"))
      .def("propagate",
           py::overload_cast<NonsingularOE &, ad::real>(
               &KeplerianDynamics::Propagate),
           py::arg("state"), py::arg("dt"))
      .def("propagate",
           py::overload_cast<EquinoctialOE &, ad::real>(
               &KeplerianDynamics::Propagate),
           py::arg("state"), py::arg("dt"))
      .def(
          "propagate_with_stm",
          [](KeplerianDynamics &dyn, ClassicalOE &state, ad::real dt) {
            Eigen::Matrix6d stm;
            dyn.PropagateWithStm(state, dt, stm);
            return stm;
          },
          py::arg("state"), py::arg("dt"), py::return_value_policy::move)
      .def(
          "propagate_with_stm",
          [](KeplerianDynamics &dyn, QuasiNonsingularOE &state, ad::real dt) {
            Eigen::Matrix6d stm;
            dyn.PropagateWithStm(state, dt, stm);
            return stm;
          },
          py::arg("state"), py::arg("dt"), py::return_value_policy::move)
      .def(
          "propagate_with_stm",
          [](KeplerianDynamics &dyn, NonsingularOE &state, ad::real dt) {
            Eigen::Matrix6d stm;
            dyn.PropagateWithStm(state, dt, stm);
            return stm;
          },
          py::arg("state"), py::arg("dt"), py::return_value_policy::move)
      .def(
          "propagate_with_stm",
          [](KeplerianDynamics &dyn, EquinoctialOE &state, ad::real dt) {
            Eigen::Matrix6d stm;
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
           py::overload_cast<OrbitState &, ad::real, ad::real, ad::real>(
               &NumericalDynamics::Propagate),
           py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("dt"))
      .def("propagate",
           py::overload_cast<ad::Vector6real &, ad::real, ad::real, ad::real>(
               &NumericalDynamics::Propagate),
           py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("dt"))
      .def(
          "propagate_with_stm",
          [](NumericalDynamics &dyn, CartesianOrbitState &state, ad::real t0,
             ad::real tf, ad::real dt) {
            Eigen::Matrix6d stm;
            dyn.PropagateWithStm(state, t0, tf, dt, stm);
            return stm;
          },
          py::arg("state"), py::arg("t0"), py::arg("tf"), py::arg("dt"),
          py::return_value_policy::move)
      .def(
          "propagate_with_stm",
          [](NumericalDynamics &dyn, ad::Vector6real &state, ad::real t0,
             ad::real tf, ad::real dt) {
            Eigen::Matrix6d stm;
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
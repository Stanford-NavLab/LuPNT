// lupnt
#include <lupnt/lupnt.h>
#include <pybind11/detail/common.h>

// pybind11
#include <string>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

using PyODE = std::function<VecXd(double, const VecXd &)>;

bool IsNotebook() {
  try {
    auto shell
        = py::module::import("IPython").attr("get_ipython")().attr("__class__").attr("__name__");
    if (shell.cast<std::string>() == "ZMQInteractiveShell") {
      return true;
    } else if (shell.cast<std::string>() == "TerminalInteractiveShell") {
      return false;
    } else {
      return false;
    }
  } catch (py::error_already_set &e) {
    return false;
  }
}

void InitDynamics(py::module &m) {
  // IntegratorType
  py::enum_<IntegratorType>(m, "IntegratorType")
      .value("RK4", IntegratorType::RK4)
      .value("RK8", IntegratorType::RK8)
      .value("RKF45", IntegratorType::RKF45)
      .value("PD45", IntegratorType::PD45)
      .export_values();

  // IntegratorParams
  py::class_<IntegratorParams>(m, "IntegratorParams")
      .def(py::init<>())
      .def(py::init<int, double, double>(), py::arg("max_iter") = 20, py::arg("abstol") = 1e-6,
           py::arg("reltol") = 1e-6)
      .def_readwrite("max_iter", &IntegratorParams::max_iter)
      .def_readwrite("abstol", &IntegratorParams::abstol)
      .def_readwrite("reltol", &IntegratorParams::reltol)
      .def(
          "set_terminate_if",
          [](IntegratorParams &p, py::function f) {
            p.terminate_if = [f](Real t, const VecX &x) -> bool {
              // Convert Eigen::VectorXd -> numpy array for Python
              py::gil_scoped_acquire gil;
              // Cast VecX -> VecXd (double) for safe conversion
              VecXd x_cast = x.cast<double>();
              py::array_t<double> arr(x_cast.size(), x_cast.data());
              return f(t.val(), arr).cast<bool>();
            };
          },
          py::arg("callback"), "Set a termination predicate f(t, x) -> bool");

  py::enum_<TerminationReason>(m, "TerminationReason")
      .value("ReachedTf", TerminationReason::ReachedTf)
      .value("UserCondition", TerminationReason::UserCondition);

  py::class_<TerminationInfo>(m, "TerminationInfo")
      .def_readonly("terminated", &TerminationInfo::terminated)
      .def_readonly("reason", &TerminationInfo::reason)
      .def_readonly("t_stop", &TerminationInfo::t_stop)
      .def_readonly("steps", &TerminationInfo::steps);

  // Dynamics
  py::class_<Dynamics>(m, "Dynamics")
      .def("set_print_progress", &Dynamics::SetPrintProgress)
      .def("propagate",
           py::overload_cast<const State &, Real, Real, const State *>(&Dynamics::Propagate),
           py::arg("x0"), py::arg("t0"), py::arg("tf"), py::arg("u") = py::none())
      .def(
          "propagate",
          [](Dynamics &dyn, const State &x0, Real t0, Real tf, const State *u = nullptr,
             bool stm = false) -> py::object {
            if (stm) {
              MatXd stm_out;
              State xf = dyn.Propagate(x0, t0, tf, u, &stm_out);
              return py::make_tuple(xf.cast<double>(), stm_out);
            } else {
              State xf2 = dyn.Propagate(x0, t0, tf, u);
              return py::cast(xf2.cast<double>());
            }
          },
          py::arg("x0"), py::arg("t0"), py::arg("tf"), py::arg("u"), py::arg("stm"))
      .def(
          "propagate",
          [](Dynamics &dyn, const State &x0, const VecX &tfs) -> MatX {
            MatX xfs = dyn.Propagate(x0, tfs);
            return xfs.cast<double>();
          },
          py::arg("x0"), py::arg("tfs"));

  // AnalyticalOrbitDynamics
  py::class_<AnalyticalDynamics, Dynamics>(m, "AnalyticalOrbitDynamics");

  // KeplerianDynamics
  py::class_<KeplerianDynamics<ClassicalOE>, AnalyticalDynamics>(m, "KeplerianDynamics")
      .def(py::init<double>(), py::arg("GM"));

  // NumericalOrbitDynamics
  py::class_<NumericalDynamics, Dynamics>(m, "NumericalOrbitDynamics")
      .def("set_ode_function",
           [](NumericalDynamics &dyn, PyODE odefunc) -> void {
             dyn.SetODE([odefunc](Real t, const VecX &x) -> VecX {
               return odefunc(t.val(), x.cast<double>()).cast<Real>();
             });
           })
      .def("get_time_step", [](NumericalDynamics &dyn) { return dyn.GetTimeStep().val(); })
      .def("set_time_step", [](NumericalDynamics &dyn, double dt) { dyn.SetTimeStep(dt); })
      .def("set_integrator",
           [](NumericalDynamics &dyn, IntegratorType integ) { dyn.SetIntegrator(integ); })
      .def("set_integrator_params", &NumericalDynamics::SetIntegratorParams)
      .def("compute_rates",
           [](NumericalDynamics &dyn, double t, const VecXd &x) -> VecXd {
             return dyn.ComputeRates(Real(t), x.cast<Real>().eval()).template cast<double>();
           })
      .def("propagate_with_info",
           [](NumericalDynamics &dyn, const State &x0, Real t0, Real tf) -> py::object {
             TerminationInfo info;
             State xf = dyn.PropagateEx(x0, t0, tf, &info);
             VecXd xf_cast = xf.cast<double>();
             return py::make_tuple(xf_cast, info);
           })
      .def("propagate_with_info",
           [](NumericalDynamics &dyn, const State &x0, Real t0, const VecX &tf) -> py::object {
             TerminationInfo info;
             MatX xf = dyn.PropagateEx(x0, t0, tf, &info);
             MatXd xf_cast = xf.cast<double>();
             return py::make_tuple(xf_cast, info);
           })
      .def("propagate_stm_with_info",
           [](NumericalDynamics &dyn, const State &x0, Real t0, Real tf) -> py::object {
             MatXd stm(6, 6);
             TerminationInfo info;
             State xf = dyn.PropagateExStm(x0, t0, tf, &stm, &info);
             VecXd xf_cast = xf.cast<double>();
             return py::make_tuple(xf_cast, stm, info);
           });

  // CartTwoBodyDynamics
  py::class_<CartesianTwoBodyDynamics, NumericalDynamics>(m, "CartesianTwoBodyDynamics")
      .def(py::init<double>(), py::arg("GM"));

  // JToCartTwoBodyDynamics
  py::class_<JToCartTwoBodyDynamics, NumericalDynamics>(m, "JToCartTwoBodyDynamics")
      .def(py::init<double, double, double>(), py::arg("GM"), py::arg("J2"), py::arg("R_body"));

  // J2KeplerianDynamics
  py::class_<J2KeplerianDynamics, NumericalDynamics>(m, "J2KeplerianDynamics")
      .def(py::init<double, double, double>(), py::arg("GM"), py::arg("J2"), py::arg("R_body"));

  // NBodyDynamics
  py::class_<NBodyDynamics, NumericalDynamics>(m, "NBodyDynamics")
      .def(py::init<>())
      .def("add_body", &NBodyDynamics::AddBody, py::arg("body"))
      .def("get_bodies", &NBodyDynamics::GetBodies)
      .def("set_frame", &NBodyDynamics::SetFrame, py::arg("frame"))
      .def("set_srp_coeff", py::overload_cast<Real, Real, Real>(&NBodyDynamics::SetSrpCoefficient),
           py::arg("CR"), py::arg("area"), py::arg("mass"))
      .def("set_drag_coeff", py::overload_cast<Real>(&NBodyDynamics::SetDragCoeff),
           py::arg("bcoeff"))
      .def("get_srp_coeff", [](NBodyDynamics &dyn) { return dyn.GetSrpCoeff().val(); })
      .def("get_drag_coeff", [](NBodyDynamics &dyn) { return dyn.GetDragCoeff().val(); })
      .def("get_use_relativity", &NBodyDynamics::GetUseRelativity)
      .def("set_use_relativity", &NBodyDynamics::SetUseRelativity, py::arg("use_relativity"))
      .def("get_units", &NBodyDynamics::GetUnits)
      .def("set_units", &NBodyDynamics::SetUnits, py::arg("units"))
      .def("set_autodiff", &NBodyDynamics::SetAutodiff)
      .def(
          "compute_accelerations",
          [](NBodyDynamics &dyn, double t, const VecXd &x, bool decompose_gravity) -> py::dict {
            std::map<std::string, Vec3> acc
                = dyn.ComputeAccelerations(Real(t), x.cast<Real>().eval(), decompose_gravity);
            py::dict out;
            for (const auto &kv : acc) {
              out[py::str(kv.first)] = kv.second.cast<double>().eval();
            }
            return out;
          },
          py::arg("t"), py::arg("x"), py::arg("decompose_gravity") = false,
          "Return a dict mapping force-term name to its Vec3 acceleration "
          "contribution at state x and epoch t (configured frame/units). When "
          "decompose_gravity is True, gravity-field bodies are resolved into "
          "individual spherical-harmonic terms (e.g. 'MOON_J2', 'MOON_C22').");

  // Body
  py::class_<Body>(m, "Body")
      .def(py::init<>())
      .def_static("Moon", py::overload_cast<int, int, std::string>(&Body::Moon), py::arg("n") = 0,
                  py::arg("m") = 0, py::arg("gravity_file") = "grgm900c.cof")
      .def_static("Moon", py::overload_cast<const UnitSystem &, int, int, std::string>(&Body::Moon),
                  py::arg("units"), py::arg("n") = 0, py::arg("m") = 0,
                  py::arg("gravity_file") = "grgm900c.cof")
      .def_static("Earth", py::overload_cast<int, int, std::string>(&Body::Earth), py::arg("n") = 0,
                  py::arg("m") = 0, py::arg("gravity_file") = "EGM96.cof")
      .def_static("Earth",
                  py::overload_cast<const UnitSystem &, int, int, std::string>(&Body::Earth),
                  py::arg("units"), py::arg("n") = 0, py::arg("m") = 0,
                  py::arg("gravity_file") = "EGM96.cof")
      .def_static("Mars", py::overload_cast<int, int, std::string>(&Body::Mars), py::arg("n") = 0,
                  py::arg("m") = 0, py::arg("gravity_file") = "GMM1.cof")
      .def_static("Mars", py::overload_cast<const UnitSystem &, int, int, std::string>(&Body::Mars),
                  py::arg("units"), py::arg("n") = 0, py::arg("m") = 0,
                  py::arg("gravity_file") = "GMM1.cof")
      .def_static("Venus", py::overload_cast<int, int, std::string>(&Body::Venus), py::arg("n") = 0,
                  py::arg("m") = 0, py::arg("gravity_file") = "MGN75HSAAP.cof")
      .def_static("Venus",
                  py::overload_cast<const UnitSystem &, int, int, std::string>(&Body::Venus),
                  py::arg("units"), py::arg("n") = 0, py::arg("m") = 0,
                  py::arg("gravity_file") = "MGN75HSAAP.cof")
      .def_static("Sun", py::overload_cast<>(&Body::Sun))
      .def_static("Sun", py::overload_cast<const UnitSystem &>(&Body::Sun), py::arg("units"))
      .def_static("Jupiter", py::overload_cast<>(&Body::Jupiter))
      .def_static("Jupiter", py::overload_cast<const UnitSystem &>(&Body::Jupiter),
                  py::arg("units"))
      .def_static("Saturn", py::overload_cast<>(&Body::Saturn))
      .def_static("Saturn", py::overload_cast<const UnitSystem &>(&Body::Saturn), py::arg("units"))
      .def_readonly("id", &Body::id)
      .def_readonly("name", &Body::name)
      .def_readonly("GM", &Body::GM)
      .def_readonly("R", &Body::R)
      .def_readonly("units", &Body::units)
      .def_readonly("fixed_frame", &Body::fixed_frame)
      .def_readonly("inertial_frame", &Body::inertial_frame)
      .def_readonly("use_gravity_field", &Body::use_gravity_field)
      .def_readonly("gravity_field", &Body::gravity_field);

  m.def("create_body", py::overload_cast<BodyId, int, int>(&CreateBody), py::arg("id"),
        py::arg("n") = 0, py::arg("m") = 0);
  m.def("create_body", py::overload_cast<BodyId, const UnitSystem &, int, int>(&CreateBody),
        py::arg("id"), py::arg("units"), py::arg("n") = 0, py::arg("m") = 0);

  // GravityField
  py::class_<GravityField<double>>(m, "GravityField")
      .def(py::init<>())
      .def_readonly("n_max", &GravityField<double>::n_max)
      .def_readonly("m_max", &GravityField<double>::m_max)
      .def_readonly("n", &GravityField<double>::n)
      .def_readonly("m", &GravityField<double>::m)
      .def_readonly("GM", &GravityField<double>::GM)
      .def_readonly("R", &GravityField<double>::R)
      .def_readonly("units", &GravityField<double>::units)
      .def_readonly("CS", &GravityField<double>::CS);
}

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
  py::enum_<IntegratorType>(m, "IntegratorType",
                            "Selects which numerical integrator NumericalOrbitDynamics uses.")
      .value("RK4", IntegratorType::RK4, "Classical 4th-order fixed-step Runge-Kutta.")
      .value("RK8", IntegratorType::RK8, "8th-order fixed-step Runge-Kutta.")
      .value("RKF45", IntegratorType::RKF45, "Runge-Kutta-Fehlberg 4(5) adaptive-step.")
      .value("PD45", IntegratorType::PD45, "Dormand-Prince 4(5) adaptive-step.")
      .export_values();

  // IntegratorParams
  py::class_<IntegratorParams>(
      m, "IntegratorParams",
      "Integrator tolerances, iteration limit, and optional early-termination predicate.")
      .def(py::init<>())
      .def(py::init<int, double, double>(), py::arg("max_iter") = 20, py::arg("abstol") = 1e-6,
           py::arg("reltol") = 1e-6,
           "Construct with max iterations, absolute and relative tolerances.")
      .def_readwrite("max_iter", &IntegratorParams::max_iter,
                     "Maximum adaptive-step iterations per step.")
      .def_readwrite("abstol", &IntegratorParams::abstol, "Absolute error tolerance.")
      .def_readwrite("reltol", &IntegratorParams::reltol, "Relative error tolerance.")
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

  py::enum_<TerminationReason>(m, "TerminationReason", "Why a PropagateEx call stopped.")
      .value("ReachedTf", TerminationReason::ReachedTf,
             "Propagation reached the requested final time.")
      .value("UserCondition", TerminationReason::UserCondition,
             "User-supplied terminate_if predicate fired.");

  py::class_<TerminationInfo>(
      m, "TerminationInfo",
      "Outcome of a PropagateEx call: whether/why it stopped, stop time, and step count.")
      .def_readonly("terminated", &TerminationInfo::terminated,
                    "True if propagation stopped before tf (user condition or max iterations).")
      .def_readonly("reason", &TerminationInfo::reason,
                    "Termination reason (see TerminationReason).")
      .def_readonly("t_stop", &TerminationInfo::t_stop, "Time actually reached [s].")
      .def_readonly("steps", &TerminationInfo::steps, "Number of integration steps taken.");

  // Dynamics
  py::class_<Dynamics>(m, "Dynamics",
                       "Base interface for state-propagation models (orbit, attitude, clock, ...).")
      .def("set_print_progress", &Dynamics::SetPrintProgress,
           "Enable or disable a progress bar during multi-step propagation.")
      .def("propagate",
           py::overload_cast<const State &, Real, Real, const State *>(&Dynamics::Propagate),
           py::arg("x0"), py::arg("t0"), py::arg("tf"), py::arg("u") = py::none(),
           "Propagate state x0 from epoch t0 to tf [s]; u is an optional control input.")
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
          py::arg("x0"), py::arg("t0"), py::arg("tf"), py::arg("u"), py::arg("stm"),
          "Propagate x0 from t0 to tf [s]; if stm=True, also return the state transition matrix.")
      .def(
          "propagate",
          [](Dynamics &dyn, const State &x0, const VecX &tfs) -> MatX {
            MatX xfs = dyn.Propagate(x0, tfs);
            return xfs.cast<double>();
          },
          py::arg("x0"), py::arg("tfs"),
          "Propagate x0 to each output epoch in tfs [s], returning states as rows.")
      .def(
          "propagate_stm",
          [](Dynamics &dyn, const VecXd &x0, double t0, double tf) -> py::tuple {
            State s(x0.cast<Real>().eval());
            MatXd stm;
            State xf = dyn.Propagate(s, Real(t0), Real(tf), nullptr, &stm);
            return py::make_tuple(xf.cast<double>().eval(), stm);
          },
          py::arg("x0"), py::arg("t0"), py::arg("tf"),
          "Propagate a numpy state vector [r; v] from t0 to tf, returning (state_tf, STM) as "
          "numpy arrays -- convenient for a filter's predict step authored in Python.");

  // AnalyticalOrbitDynamics
  py::class_<AnalyticalDynamics, Dynamics>(
      m, "AnalyticalOrbitDynamics", "Base class for analytical (closed-form) orbit propagation.");

  // KeplerianDynamics
  py::class_<KeplerianDynamics<ClassicalOE>, AnalyticalDynamics>(
      m, "KeplerianDynamics", "Analytical Keplerian propagation of classical orbital elements.")
      .def(py::init<double>(), py::arg("GM"), "Construct with central-body GM [m^3/s^2].");

  // NumericalOrbitDynamics
  py::class_<NumericalDynamics, Dynamics>(
      m, "NumericalOrbitDynamics",
      "Numerically integrated orbit dynamics driven by an ODE right-hand side and an integrator.")
      .def(
          "set_ode_function",
          [](NumericalDynamics &dyn, PyODE odefunc) -> void {
            dyn.SetODE([odefunc](Real t, const VecX &x) -> VecX {
              return odefunc(t.val(), x.cast<double>()).cast<Real>();
            });
          },
          "Set the ODE right-hand side f(t, x) -> dx/dt.")
      .def(
          "get_time_step", [](NumericalDynamics &dyn) { return dyn.GetTimeStep().val(); },
          "Return the integration time step [s].")
      .def(
          "set_time_step", [](NumericalDynamics &dyn, double dt) { dyn.SetTimeStep(dt); },
          "Set the integration time step [s].")
      .def(
          "set_integrator",
          [](NumericalDynamics &dyn, IntegratorType integ) { dyn.SetIntegrator(integ); },
          "Select the integrator type (see IntegratorType).")
      .def("set_integrator_params", &NumericalDynamics::SetIntegratorParams,
           "Set integrator tolerances, iteration limit, and termination predicate.")
      .def(
          "compute_rates",
          [](NumericalDynamics &dyn, double t, const VecXd &x) -> VecXd {
            return dyn.ComputeRates(Real(t), x.cast<Real>().eval()).template cast<double>();
          },
          "Evaluate the state derivative dx/dt at epoch t [s] and state x.")
      .def(
          "propagate_with_info",
          [](NumericalDynamics &dyn, const State &x0, Real t0, Real tf) -> py::object {
            TerminationInfo info;
            State xf = dyn.PropagateEx(x0, t0, tf, &info);
            VecXd xf_cast = xf.cast<double>();
            return py::make_tuple(xf_cast, info);
          },
          "Propagate x0 from t0 to tf [s], returning (state_tf, TerminationInfo).")
      .def(
          "propagate_with_info",
          [](NumericalDynamics &dyn, const State &x0, Real t0, const VecX &tf) -> py::object {
            TerminationInfo info;
            MatX xf = dyn.PropagateEx(x0, t0, tf, &info);
            MatXd xf_cast = xf.cast<double>();
            return py::make_tuple(xf_cast, info);
          },
          "Propagate x0 from t0 to each epoch in tf [s], returning (states, TerminationInfo).")
      .def(
          "propagate_stm_with_info",
          [](NumericalDynamics &dyn, const State &x0, Real t0, Real tf) -> py::object {
            MatXd stm(6, 6);
            TerminationInfo info;
            State xf = dyn.PropagateExStm(x0, t0, tf, &stm, &info);
            VecXd xf_cast = xf.cast<double>();
            return py::make_tuple(xf_cast, stm, info);
          },
          "Propagate x0 from t0 to tf [s], returning (state_tf, STM, TerminationInfo).");

  // CartTwoBodyDynamics
  py::class_<CartesianTwoBodyDynamics, NumericalDynamics>(
      m, "CartesianTwoBodyDynamics", "Numerical Cartesian two-body (point-mass) orbit dynamics.")
      .def(py::init<double>(), py::arg("GM"), "Construct with central-body GM [m^3/s^2].");

  // JToCartTwoBodyDynamics
  py::class_<JToCartTwoBodyDynamics, NumericalDynamics>(
      m, "JToCartTwoBodyDynamics",
      "Numerical Cartesian two-body dynamics with J2 oblateness perturbation.")
      .def(py::init<double, double, double>(), py::arg("GM"), py::arg("J2"), py::arg("R_body"),
           "Construct with GM [m^3/s^2], J2 [-], and body reference radius [m].");

  // J2KeplerianDynamics
  py::class_<J2KeplerianDynamics, NumericalDynamics>(
      m, "J2KeplerianDynamics", "Numerical Keplerian dynamics with secular J2 perturbation.")
      .def(py::init<double, double, double>(), py::arg("GM"), py::arg("J2"), py::arg("R_body"),
           "Construct with GM [m^3/s^2], J2 [-], and body reference radius [m].");

  // NBodyDynamics
  py::class_<NBodyDynamics, NumericalDynamics>(
      m, "NBodyDynamics",
      "Numerical point-mass, gravity-field, SRP, drag, and relativity orbit dynamics.")
      .def(py::init<>())
      .def("add_body", &NBodyDynamics::AddBody, py::arg("body"),
           "Add a central or perturbing body to the force model.")
      .def("get_bodies", &NBodyDynamics::GetBodies,
           "Return the bodies included in the force model.")
      .def("set_frame", &NBodyDynamics::SetFrame, py::arg("frame"),
           "Set the frame used for numerical integration.")
      .def("set_srp_coeff", py::overload_cast<Real, Real, Real>(&NBodyDynamics::SetSrpCoefficient),
           py::arg("CR"), py::arg("area"), py::arg("mass"),
           "Set the SRP coefficient from reflectivity CR [-], area [m^2], and mass [kg].")
      .def("set_drag_coeff", py::overload_cast<Real>(&NBodyDynamics::SetDragCoeff),
           py::arg("bcoeff"), "Set the drag ballistic coefficient directly [m^2/kg].")
      .def(
          "get_srp_coeff", [](NBodyDynamics &dyn) { return dyn.GetSrpCoeff().val(); },
          "Return the SRP ballistic coefficient [m^2/kg].")
      .def(
          "get_drag_coeff", [](NBodyDynamics &dyn) { return dyn.GetDragCoeff().val(); },
          "Return the drag ballistic coefficient [m^2/kg].")
      .def("get_use_relativity", &NBodyDynamics::GetUseRelativity,
           "Return whether post-Newtonian relativistic corrections are enabled.")
      .def("set_use_relativity", &NBodyDynamics::SetUseRelativity, py::arg("use_relativity"),
           "Enable or disable post-Newtonian relativistic corrections.")
      .def("get_units", &NBodyDynamics::GetUnits,
           "Return the active distance/time/mass unit system.")
      .def("set_units", &NBodyDynamics::SetUnits, py::arg("units"),
           "Set the active unit system for states, constants, and acceleration output.")
      .def("set_autodiff", &NBodyDynamics::SetAutodiff,
           "Enable or disable automatic-differentiation rates.")
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
  py::class_<Body>(
      m, "Body",
      "A celestial body: gravity parameter, radius, rotation rate, frames, optional gravity field.")
      .def(py::init<>())
      .def_static("Moon", py::overload_cast<int, int, std::string>(&Body::Moon), py::arg("n") = 0,
                  py::arg("m") = 0, py::arg("gravity_file") = "grgm900c.cof",
                  "Create the Moon, optionally loading an n x m spherical-harmonic gravity field.")
      .def_static(
          "Moon", py::overload_cast<const UnitSystem &, int, int, std::string>(&Body::Moon),
          py::arg("units"), py::arg("n") = 0, py::arg("m") = 0,
          py::arg("gravity_file") = "grgm900c.cof",
          "Create the Moon in the given unit system, optionally with an n x m gravity field.")
      .def_static("Earth", py::overload_cast<int, int, std::string>(&Body::Earth), py::arg("n") = 0,
                  py::arg("m") = 0, py::arg("gravity_file") = "EGM96.cof",
                  "Create Earth, optionally loading an n x m spherical-harmonic gravity field.")
      .def_static("Earth",
                  py::overload_cast<const UnitSystem &, int, int, std::string>(&Body::Earth),
                  py::arg("units"), py::arg("n") = 0, py::arg("m") = 0,
                  py::arg("gravity_file") = "EGM96.cof",
                  "Create Earth in the given unit system, optionally with an n x m gravity field.")
      .def_static("Mars", py::overload_cast<int, int, std::string>(&Body::Mars), py::arg("n") = 0,
                  py::arg("m") = 0, py::arg("gravity_file") = "GMM1.cof",
                  "Create Mars, optionally loading an n x m spherical-harmonic gravity field.")
      .def_static("Mars", py::overload_cast<const UnitSystem &, int, int, std::string>(&Body::Mars),
                  py::arg("units"), py::arg("n") = 0, py::arg("m") = 0,
                  py::arg("gravity_file") = "GMM1.cof",
                  "Create Mars in the given unit system, optionally with an n x m gravity field.")
      .def_static("Venus", py::overload_cast<int, int, std::string>(&Body::Venus), py::arg("n") = 0,
                  py::arg("m") = 0, py::arg("gravity_file") = "MGN75HSAAP.cof",
                  "Create Venus, optionally loading an n x m spherical-harmonic gravity field.")
      .def_static("Venus",
                  py::overload_cast<const UnitSystem &, int, int, std::string>(&Body::Venus),
                  py::arg("units"), py::arg("n") = 0, py::arg("m") = 0,
                  py::arg("gravity_file") = "MGN75HSAAP.cof",
                  "Create Venus in the given unit system, optionally with an n x m gravity field.")
      .def_static("Sun", py::overload_cast<>(&Body::Sun), "Create the Sun as a point-mass body.")
      .def_static("Sun", py::overload_cast<const UnitSystem &>(&Body::Sun), py::arg("units"),
                  "Create the Sun as a point-mass body in the given unit system.")
      .def_static("Jupiter", py::overload_cast<>(&Body::Jupiter),
                  "Create Jupiter as a point-mass body.")
      .def_static("Jupiter", py::overload_cast<const UnitSystem &>(&Body::Jupiter),
                  py::arg("units"), "Create Jupiter as a point-mass body in the given unit system.")
      .def_static("Saturn", py::overload_cast<>(&Body::Saturn),
                  "Create Saturn as a point-mass body.")
      .def_static("Saturn", py::overload_cast<const UnitSystem &>(&Body::Saturn), py::arg("units"),
                  "Create Saturn as a point-mass body in the given unit system.")
      .def_readonly("id", &Body::id, "Body identifier.")
      .def_readonly("name", &Body::name, "Body display name.")
      .def_readonly("GM", &Body::GM, "Gravitational parameter [m^3/s^2].")
      .def_readonly("R", &Body::R, "Reference radius [m].")
      .def_readonly("units", &Body::units, "Unit system of the body's constants.")
      .def_readonly("fixed_frame", &Body::fixed_frame, "Body-fixed reference frame.")
      .def_readonly("inertial_frame", &Body::inertial_frame, "Inertial reference frame.")
      .def_readonly("use_gravity_field", &Body::use_gravity_field,
                    "True if a spherical-harmonic gravity field is loaded.")
      .def_readonly("gravity_field", &Body::gravity_field,
                    "Spherical-harmonic gravity field coefficients and metadata.");

  m.def("create_body", py::overload_cast<BodyId, int, int>(&CreateBody), py::arg("id"),
        py::arg("n") = 0, py::arg("m") = 0,
        "Create a body from its identifier, optionally with an n x m gravity field.");
  m.def("create_body", py::overload_cast<BodyId, const UnitSystem &, int, int>(&CreateBody),
        py::arg("id"), py::arg("units"), py::arg("n") = 0, py::arg("m") = 0,
        "Create a body in the given unit system, optionally with an n x m gravity field.");

  // GravityField
  py::class_<GravityField<double>>(m, "GravityField",
                                   "Spherical-harmonic gravity field coefficients and metadata.")
      .def(py::init<>())
      .def_readonly("n_max", &GravityField<double>::n_max, "Maximum degree available.")
      .def_readonly("m_max", &GravityField<double>::m_max, "Maximum order available.")
      .def_readonly("n", &GravityField<double>::n, "Degree retained.")
      .def_readonly("m", &GravityField<double>::m, "Order retained.")
      .def_readonly("GM", &GravityField<double>::GM, "Gravitational parameter [m^3/s^2].")
      .def_readonly("R", &GravityField<double>::R, "Reference radius [m].")
      .def_readonly("units", &GravityField<double>::units, "Unit system of the field constants.")
      .def_readonly("CS", &GravityField<double>::CS,
                    "Unnormalized spherical-harmonic coefficients (C on/above diagonal, S below).");
}

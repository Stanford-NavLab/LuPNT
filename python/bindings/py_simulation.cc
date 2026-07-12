#include <lupnt/agents/agent.h>
#include <lupnt/applications/application.h>
#include <lupnt/lupnt.h>
#include <lupnt/simulations/simulation.h>
#include <lupnt/simulations/world.h>

#include <memory>
#include <string>

#include "py_pybind11.h"
#include "py_yaml.h"  // NOLINT - YAML::Node type caster (accepts dicts)

namespace py = pybind11;
using namespace lupnt;

// Trampoline so an Application subclass can be authored in pure Python: the C++ virtuals
// dispatch to Python overrides. `Step` is pure-virtual (every app must implement it);
// `Setup`/`Log` have base implementations (call `super().setup()` to keep the base's
// frequency-based Step scheduling). See `register_application` below for the factory hook.
class PyApplication : public Application {
public:
  using Application::Application;
  void Setup() override { PYBIND11_OVERRIDE_NAME(void, Application, "setup", Setup); }
  // Pass time as a plain Python float (t.val()), not a Real/autodiff object, so Python
  // authors can do ordinary arithmetic on it.
  void Step(Real t) override {
    py::gil_scoped_acquire gil;
    py::function f = py::get_override(static_cast<const Application*>(this), "step");
    if (!f) throw std::runtime_error("Application subclass must implement step(self, t)");
    f(t.val());
  }
  void Log(Real t) override {
    py::gil_scoped_acquire gil;
    py::function f = py::get_override(static_cast<const Application*>(this), "log");
    if (f)
      f(t.val());
    else
      Application::Log(t);
  }
};

// Trampoline so an Agent (a physical platform) can be authored in pure Python: the pure-virtual
// GetStateAt(t) dispatches to a Python `get_state_at(self, t)` returning a numpy [r; v] 6-vector;
// Setup/Step/Log optionally dispatch to Python. Use for a truth trajectory defined in Python
// (e.g. an analytic ephemeris). See `register_agent` below.
class PyAgent : public Agent {
public:
  using Agent::Agent;
  Cart6 GetStateAt(Real t) const override {
    py::gil_scoped_acquire gil;
    py::function f = py::get_override(static_cast<const Agent*>(this), "get_state_at");
    if (!f) throw std::runtime_error("Agent subclass must implement get_state_at(self, t)");
    VecXd v = f(t.val()).cast<VecXd>();
    LUPNT_CHECK(v.size() == 6, "Agent.get_state_at must return a 6-vector [r; v]", "PyAgent");
    Vec6 x6 = v.cast<Real>();
    return Cart6(x6);
  }
  void Setup() override { PYBIND11_OVERRIDE_NAME(void, Agent, "setup", Setup); }
  void Step(Real t) override {
    py::gil_scoped_acquire gil;
    py::function f = py::get_override(static_cast<const Agent*>(this), "step");
    if (f)
      f(t.val());
    else
      Agent::Step(t);
  }
  void Log(Real t) override {
    py::gil_scoped_acquire gil;
    py::function f = py::get_override(static_cast<const Agent*>(this), "log");
    if (f)
      f(t.val());
    else
      Agent::Log(t);
  }
};

// Trampoline so a Measurement model can be authored in pure Python: the pure-virtual
// Compute(x) dispatches to a Python `compute(self, x)` that returns ``(z, H, R)`` --- the
// predicted observable, its Jacobian dh/dx, and the noise covariance --- as numpy arrays.
// Applying the model to a truth state (from an agent/device) and adding noise generates an
// observation; applying it to the filter state predicts one. Clone() shares the Python object
// (models are read-only) so a filter's CreateFunction() closure keeps dispatching to Python.
class PyMeasurement : public Measurement {
public:
  using Measurement::Measurement;
  MeasData Compute(const State& x, MatXd* H = nullptr) const override {
    py::gil_scoped_acquire gil;
    py::function f = py::get_override(this, "compute");
    if (!f) throw std::runtime_error("Measurement subclass must implement compute(self, x)");
    py::tuple t = py::reinterpret_steal<py::tuple>(f(x.cast<double>().eval()).release());
    MeasData md;
    md.value = t[0].cast<VecXd>();
    if (H) *H = t[1].cast<MatXd>();
    md.covariance = t[2].cast<MatXd>();
    return md;
  }
  Ptr<Measurement> Clone() const override {
    py::gil_scoped_acquire gil;
    auto holder = std::make_shared<py::object>(py::cast(this));
    return std::shared_ptr<Measurement>(holder, const_cast<PyMeasurement*>(this));
  }
};

// Recursively convert a YAML::Node to a native Python object (dict / list / int / float /
// bool / str) for handing an app's config block to a Python-authored Application. Scalars are
// typed via yaml-cpp (int, then float -- which correctly parses scientific notation like
// 1e-8 -- then bool, else str), avoiding the YAML 1.1 round-trip that pyyaml mis-reads.
static py::object YamlToPy(const YAML::Node& node) {
  switch (node.Type()) {
    case YAML::NodeType::Map: {
      py::dict d;
      for (const auto& kv : node) d[py::str(kv.first.as<std::string>())] = YamlToPy(kv.second);
      return d;
    }
    case YAML::NodeType::Sequence: {
      py::list l;
      for (const auto& e : node) l.append(YamlToPy(e));
      return std::move(l);
    }
    case YAML::NodeType::Scalar: {
      long long i;
      if (YAML::convert<long long>::decode(node, i)) return py::int_(i);
      double dbl;
      if (YAML::convert<double>::decode(node, dbl)) return py::float_(dbl);
      bool b;
      if (YAML::convert<bool>::decode(node, b)) return py::bool_(b);
      return py::str(node.Scalar());
    }
    default: return py::none();
  }
}

// Framework-level agent-based simulation API: World, Application, Agent, Simulation.
// Application/Agent are registered as polymorphic bases with shared_ptr holders so
// that Agent::get_application() automatically downcasts to concrete Application
// subclasses (e.g. GroundStationManagerApp) registered elsewhere.
void InitSimulation(py::module& m) {
  // ---- World: shared read-only environment ----
  py::class_<World, std::shared_ptr<World>>(
      m, "World",
      "Shared, read-only physical environment for a Simulation: epoch, reference frame, force "
      "model, optional terrain/DEM, and a truth facade (get_state_at).")
      .def(
          "get_epoch", [](const World& w) { return w.GetEpoch().val(); },
          "TDB epoch [s past J2000] of simulation time t = 0")
      .def("get_frame", &World::GetFrame, "Reference frame of the shared environment")
      .def("has_force_model", &World::HasForceModel,
           "Whether a force_model block was provided (so an NBodyDynamics can be built)")
      .def("get_gm", &World::GetGM, "Central-body gravitational parameter [m^3/s^2]")
      .def(
          "gravity", [](const World& w, const Vec3d& r) { return w.Gravity(r); }, py::arg("r"),
          "Point-mass central-body gravitational acceleration at r (world frame) [m/s^2]")
      .def("has_terrain", &World::HasTerrain, "Whether a dem (terrain) block was provided")
      .def("get_elevation", &World::GetElevation, py::arg("east_m"), py::arg("north_m"),
           "Terrain elevation at a local ENU offset from the site center [m]")
      .def("enu_to_world", &World::EnuToWorld, py::arg("east_m"), py::arg("north_m"),
           py::arg("up_m"), "Convert a local ENU offset to a world-frame position [m]")
      .def(
          "r_enu_to_world", [](const World& w) { return w.REnuToWorld(); },
          "Rotation from the local ENU tangent frame to the world frame")
      .def(
          "site_center_world", [](const World& w) { return w.SiteCenterWorld(); },
          "Site center (local ENU origin) position in the world frame [m]")
      .def("site_lat_deg", &World::SiteLatDeg, "Site latitude [deg]")
      .def("site_lon_deg", &World::SiteLonDeg, "Site east longitude [deg]")
      .def(
          "dem_x", [](const World& w) { return w.GetDem().x(); },
          "Terrain DEM grid x coordinates (native projected meters)")
      .def(
          "dem_y", [](const World& w) { return w.GetDem().y(); },
          "Terrain DEM grid y coordinates (native projected meters)")
      .def(
          "dem_elevation", [](const World& w) { return w.GetDem().elevation(); },
          "Terrain DEM elevation grid [m]")
      .def(
          "dem_center",
          [](const World& w) { return Vec2d(w.GetDem().center_x(), w.GetDem().center_y()); },
          "DEM tile center [x, y] in native projected meters (the ENU East/North origin)")
      .def(
          "get_state_at",
          [](const World& w, const std::string& name, double t) {
            return w.GetStateAt(name, t).cast<double>().eval();
          },
          py::arg("name"), py::arg("t"),
          "Truth state [r; v] of agent `name` at simulation time t [s], in the world frame");

  // ---- Application: polymorphic base (subclassable from Python via PyApplication) ----
  py::class_<Application, PyApplication, std::shared_ptr<Application>>(
      m, "Application",
      "Polymorphic base for a scenario application (estimator/logic) hosted on an Agent. "
      "Subclass in Python and implement step(self, t) (and optionally setup/log).")
      .def(py::init<>(), "Construct an empty Application (base for a Python subclass).")
      .def("get_name", &Application::GetName, "This application's name")
      .def("set_name", &Application::SetName, py::arg("name"), "Set this application's name")
      .def(
          "get_frequency", [](const Application& a) { return a.GetFrequency().val(); },
          "Step() call frequency [Hz]")
      .def(
          "set_frequency", [](Application& a, double f) { a.SetFrequency(f); },
          py::arg("frequency"),
          "Set the Step() call frequency [Hz] used by setup() to schedule steps")
      .def("get_agent", &Application::GetAgent, py::return_value_policy::reference,
           "The Agent that hosts this application (set by the simulation before Setup()).")
      .def("setup", &Application::Setup,
           "Base Setup: schedules Step() at get_frequency() Hz. Call via super().setup() from a "
           "Python subclass to keep that scheduling.")
      .def("log", &Application::Log, py::arg("t"),
           "Base Log at simulation time t [s] (emits a debug message; override in a subclass).");

  // ---- Measurement: polymorphic model base (subclassable from Python via PyMeasurement) ----
  // Subclass and implement ``compute(self, x) -> (z, H, R)`` (numpy). ``evaluate(x)`` runs the
  // model through the C++ base (the same path a Filter uses), so the same Python class both
  // generates observations (apply to a truth state, add noise) and predicts them (apply to the
  // filter state).
  py::class_<Measurement, PyMeasurement, std::shared_ptr<Measurement>>(
      m, "Measurement",
      "Polymorphic base for a measurement model. Subclass in Python and implement "
      "compute(self, x) -> (z, H = dh/dx, R) to both generate and predict observations.")
      .def(py::init<>(), "Construct an empty Measurement (base for a Python subclass).")
      .def(
          "evaluate",
          [](const Measurement& meas, const VecXd& x) {
            State s(x.cast<Real>().eval());
            MatXd H;
            MeasData md = meas.Compute(s, &H);
            return py::make_tuple(md.value, H, md.covariance);
          },
          py::arg("x"),
          "Evaluate the measurement model at state x (numpy [r; v; ...]); returns "
          "(z, H = dh/dx, R). Dispatches to a Python subclass's compute(self, x).");

  // ---- Agent: polymorphic base (subclassable from Python via PyAgent) ----
  py::class_<Agent, PyAgent, std::shared_ptr<Agent>>(
      m, "Agent",
      "Polymorphic base for a physical platform (spacecraft, rover, lander) that hosts "
      "Applications. Subclass in Python and implement get_state_at(self, t) for a truth "
      "trajectory.")
      .def(py::init<>(), "Construct an empty Agent (base for a Python subclass).")
      .def("get_name", &Agent::GetName, "This agent's name")
      .def("set_name", &Agent::SetName, py::arg("name"), "Set this agent's name")
      .def("get_application", &Agent::GetApplication,
           "The primary Application hosted by this agent (the first one; downcasts to the "
           "concrete app type)")
      .def("get_applications", &Agent::GetApplications,
           "All Applications hosted by this agent, in order (for multi-application agents)")
      .def("get_application_by_name", &Agent::GetApplicationByName, py::arg("name"),
           "The hosted Application whose name ends with `name` (e.g. \"LanderNavApp\"), or None")
      .def(
          "get_state_at",
          [](const Agent& a, double t) { return a.GetStateAt(t).cast<double>().eval(); },
          py::arg("t"), "Cartesian state [r; v] of this agent at simulation time t [s]")
      .def("get_world", &Agent::GetWorld, py::return_value_policy::reference,
           "The shared World environment (use world.get_state_at(name, t) to read any agent's "
           "truth state).");

  // ---- Simulation: holds agents, the world, and the event queue ----
  py::class_<Simulation>(
      m, "Simulation",
      "Agent-based simulation: holds the agents, the shared World, and the timed event queue; "
      "run() drives the event loop to completion.")
      .def(py::init([](const std::string& config_path) {
             YAML::Node node = YAML::LoadFile(config_path);
             Config cfg(node);
             return std::make_unique<Simulation>(cfg);
           }),
           py::arg("config_path"),
           "Build a Simulation from a YAML scenario file (agents + world + event queue).")
      .def(py::init([](YAML::Node node) {
             Config cfg(node);
             return std::make_unique<Simulation>(cfg);
           }),
           py::arg("config"), "Build a Simulation from a config dict / YAML node.")
      .def("run", &Simulation::Run, "Run the event loop to completion.")
      .def("get_agent", &Simulation::GetAgent, py::arg("name"),
           py::return_value_policy::reference_internal, "Get the agent registered under `name`")
      .def("get_world", &Simulation::GetWorld, py::return_value_policy::reference_internal,
           "The shared read-only World environment (None if no world: block was defined)")
      .def(
          "get_duration", [](Simulation& s) { return s.GetDuration().val(); },
          "Total scenario duration [s]")
      .def(
          "get_time", [](Simulation& s) { return s.GetTime().val(); },
          "Current simulation time [s]");

  // ---- Author agents/apps in Python: register a Python subclass with the asset factory ----
  // After `register_application("MyApp", MyApp)`, a config `application: {class: MyApp, ...}` is
  // instantiated by calling `MyApp(config_dict)` (the `application:` block, as a dict). The
  // returned Python object is held by the C++ simulation; `agent.get_application()` hands the
  // SAME Python object back, so results stored on `self` are readable after `sim.run()`.
  m.def(
      "register_application",
      [](const std::string& name, py::object cls) {
        AssetFactory<Application, Config&>::Register(
            name, [cls](Config& config) -> std::shared_ptr<Application> {
              py::gil_scoped_acquire gil;
              py::object cfg = YamlToPy(config);
              py::object obj = cls(cfg);
              // Keep the Python instance alive for the C++ simulation's lifetime, and tie the
              // returned shared_ptr's lifetime to that Python object (so trampoline dispatch to
              // the Python overrides keeps working). The Application is aliased out of it. The
              // custom deleter re-acquires the GIL so the py::object is freed safely wherever the
              // C++ shared_ptr is released.
              auto holder = std::shared_ptr<py::object>(new py::object(obj), [](py::object* p) {
                py::gil_scoped_acquire gil;
                delete p;
              });
              Application* app = obj.cast<Application*>();
              return std::shared_ptr<Application>(holder, app);
            });
      },
      py::arg("name"), py::arg("cls"),
      "Register a Python Application subclass under `name` so a config `class: name` builds it. "
      "The class is constructed as `cls(config_dict)`; subclass Step(t) (and optionally "
      "Setup()/Log(t)), read truth via get_agent().get_world().get_state_at(agent, t), and store "
      "results on self.");

  // ---- Author agents in Python: register a Python Agent subclass with the asset factory ----
  m.def(
      "register_agent",
      [](const std::string& name, py::object cls) {
        AssetFactory<Agent, Config&>::Register(
            name, [cls](Config& config) -> std::shared_ptr<Agent> {
              py::gil_scoped_acquire gil;
              py::object obj = cls(YamlToPy(config));
              auto holder = std::shared_ptr<py::object>(new py::object(obj), [](py::object* p) {
                py::gil_scoped_acquire g;
                delete p;
              });
              return std::shared_ptr<Agent>(holder, obj.cast<Agent*>());
            });
      },
      py::arg("name"), py::arg("cls"),
      "Register a Python Agent subclass under `name` so a config `class: name` builds it. The "
      "class is constructed as `cls(config_dict)` and must implement get_state_at(self, t) "
      "returning a numpy [r; v] 6-vector (e.g. an analytic truth trajectory).");
}

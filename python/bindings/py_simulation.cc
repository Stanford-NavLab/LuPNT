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
  py::class_<World, std::shared_ptr<World>>(m, "World")
      .def(
          "get_epoch", [](const World& w) { return w.GetEpoch().val(); },
          "TDB epoch [s past J2000] of simulation time t = 0")
      .def("get_frame", &World::GetFrame, "Reference frame of the shared environment")
      .def("has_force_model", &World::HasForceModel)
      .def("get_gm", &World::GetGM, "Central-body gravitational parameter [m^3/s^2]")
      .def(
          "gravity", [](const World& w, const Vec3d& r) { return w.Gravity(r); }, py::arg("r"),
          "Point-mass central-body gravitational acceleration at r (world frame) [m/s^2]")
      .def("has_terrain", &World::HasTerrain)
      .def("get_elevation", &World::GetElevation, py::arg("east_m"), py::arg("north_m"),
           "Terrain elevation at a local ENU offset from the site center [m]")
      .def("enu_to_world", &World::EnuToWorld, py::arg("east_m"), py::arg("north_m"),
           py::arg("up_m"), "Convert a local ENU offset to a world-frame position [m]")
      .def(
          "r_enu_to_world", [](const World& w) { return w.REnuToWorld(); },
          "Rotation from the local ENU tangent frame to the world frame")
      .def("site_center_world", [](const World& w) { return w.SiteCenterWorld(); })
      .def("site_lat_deg", &World::SiteLatDeg)
      .def("site_lon_deg", &World::SiteLonDeg)
      .def("dem_x", [](const World& w) { return w.GetDem().x(); })
      .def("dem_y", [](const World& w) { return w.GetDem().y(); })
      .def("dem_elevation", [](const World& w) { return w.GetDem().elevation(); })
      .def("dem_center",
           [](const World& w) { return Vec2d(w.GetDem().center_x(), w.GetDem().center_y()); })
      .def(
          "get_state_at",
          [](const World& w, const std::string& name, double t) {
            return w.GetStateAt(name, t).cast<double>().eval();
          },
          py::arg("name"), py::arg("t"),
          "Truth state [r; v] of agent `name` at simulation time t [s], in the world frame");

  // ---- Application: polymorphic base (subclassable from Python via PyApplication) ----
  py::class_<Application, PyApplication, std::shared_ptr<Application>>(m, "Application")
      .def(py::init<>())
      .def("get_name", &Application::GetName)
      .def("set_name", &Application::SetName, py::arg("name"))
      .def("get_frequency", [](const Application& a) { return a.GetFrequency().val(); })
      .def(
          "set_frequency", [](Application& a, double f) { a.SetFrequency(f); },
          py::arg("frequency"))
      .def("get_agent", &Application::GetAgent, py::return_value_policy::reference,
           "The Agent that hosts this application (set by the simulation before Setup()).")
      .def("setup", &Application::Setup,
           "Base Setup: schedules Step() at get_frequency() Hz. Call via super().setup() from a "
           "Python subclass to keep that scheduling.")
      .def("log", &Application::Log, py::arg("t"));

  // ---- Agent: polymorphic base ----
  py::class_<Agent, std::shared_ptr<Agent>>(m, "Agent")
      .def("get_name", &Agent::GetName)
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
  py::class_<Simulation>(m, "Simulation")
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
           py::return_value_policy::reference_internal)
      .def("get_world", &Simulation::GetWorld, py::return_value_policy::reference_internal)
      .def("get_duration", [](Simulation& s) { return s.GetDuration().val(); })
      .def("get_time", [](Simulation& s) { return s.GetTime().val(); });

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
}

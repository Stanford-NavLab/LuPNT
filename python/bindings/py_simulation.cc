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

  // ---- Application: polymorphic base ----
  py::class_<Application, std::shared_ptr<Application>>(m, "Application")
      .def("get_name", &Application::GetName)
      .def("get_frequency", [](const Application& a) { return a.GetFrequency().val(); });

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
          py::arg("t"), "Cartesian state [r; v] of this agent at simulation time t [s]");

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
}

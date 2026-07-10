/**
 * @file py_lander_nav.cc
 * @brief Python bindings for the agent-based lunar-lander navigation application
 *        (`lupnt::LanderNavApp`, `lupnt/applications/lander/lander_nav_app.h`).
 *
 * The lander scenario runs through the generic agent framework:
 *   sim = pnt.Simulation("configs/lander_nav.yaml")
 *   app = sim.get_agent("Lander").get_application()   # -> LanderNavApp
 *   app.set_reference_trajectory_enu(traj)            # optional guidance-law truth path
 *   sim.run()
 * The app precomputes the descent truth, crater map, and relay orbits, synthesizes the IMU +
 * pseudorange + altimeter + crater-bearing measurements from the shared World, runs the MEKF,
 * and exposes the per-epoch truth/estimate/covariance series through the accessors below.
 *
 * `LanderNavApp` is registered as an `Application` subclass (polymorphic base bound in
 * `py_simulation.cc`) so `Agent.get_application()` downcasts automatically.
 */
#include <lupnt/lupnt.h>

#include <memory>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitLanderNav(py::module& m) {
  py::class_<LanderNavApp, Application, std::shared_ptr<LanderNavApp>>(m, "LanderNavApp")
      .def("set_reference_trajectory_enu", &LanderNavApp::SetReferenceTrajectoryEnu,
           py::arg("ref_traj_enu"),
           "Supply an [N x 3] ENU reference (truth) trajectory [m] about the DEM site center "
           "(U = height above the site datum). Overrides the built-in smoothstep descent and "
           "sets N; call before sim.run().")
      .def("site_id", &LanderNavApp::site_id)
      .def("site_name", &LanderNavApp::site_name)
      .def("dem_x", &LanderNavApp::dem_x)
      .def("dem_y", &LanderNavApp::dem_y)
      .def("dem_elevation", &LanderNavApp::dem_elevation)
      .def("crater_enu", &LanderNavApp::crater_enu, "[M x 2] crater (East, North) [m]")
      .def("time_s", &LanderNavApp::time_series)
      .def("pos_err_enu", &LanderNavApp::pos_err_enu,
           "[N x 3] (truth - est) position error in ENU [m]")
      .def("pos_sigma_enu", &LanderNavApp::pos_sigma_enu)
      .def("pos_err_norm", &LanderNavApp::pos_err_norm)
      .def("vel_err_norm", &LanderNavApp::vel_err_norm)
      .def("clock_bias_err", &LanderNavApp::clock_bias_err)
      .def("clock_bias_sigma", &LanderNavApp::clock_bias_sigma)
      .def("n_visible_sat", &LanderNavApp::n_visible_sat)
      .def("n_craters", &LanderNavApp::n_craters)
      .def("accel_bias_err", &LanderNavApp::accel_bias_err)
      .def("accel_bias_sigma", &LanderNavApp::accel_bias_sigma)
      .def("gyro_bias_err", &LanderNavApp::gyro_bias_err)
      .def("gyro_bias_sigma", &LanderNavApp::gyro_bias_sigma)
      .def("att_err_deg", &LanderNavApp::att_err_deg,
           "[N x 3] attitude error (rotation vector) [deg]")
      .def("att_sigma_deg", &LanderNavApp::att_sigma_deg)
      .def("att_quat_est", &LanderNavApp::att_quat_est,
           "[N x 4] estimated body-to-nav attitude quaternion [w,x,y,z] (MEKF)")
      .def("traj_enu_truth", &LanderNavApp::traj_enu_truth, "[N x 3] truth (East, North, Up) [m]")
      .def("traj_enu_est", &LanderNavApp::traj_enu_est, "[N x 3] estimated (East, North, Up) [m]")
      .def("alt_truth", &LanderNavApp::alt_truth, "[N] truth height above terrain [m]")
      .def("alt_est", &LanderNavApp::alt_est, "[N] estimated height above terrain [m]")
      .def("satellite_names", &LanderNavApp::satellite_names);
}

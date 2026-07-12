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
  // Guidance/control app: owns the descent truth trajectory (and the reference-trajectory hook).
  py::class_<LanderGncApp, Application, std::shared_ptr<LanderGncApp>>(
      m, "LanderGncApp",
      "Co-hosted lander guidance/control application that owns the descent truth trajectory "
      "(built-in smoothstep descent or a supplied ENU reference path).")
      .def("set_reference_trajectory_enu", &LanderGncApp::SetReferenceTrajectoryEnu,
           py::arg("ref_traj_enu"),
           "Supply an [N x 3] ENU reference (truth) trajectory [m] about the DEM site center "
           "(U = height above the site datum). Overrides the built-in smoothstep descent and "
           "sets N; call before sim.run().")
      .def("traj_enu_truth", &LanderGncApp::traj_enu_truth, "[N x 3] truth (East, North, Up) [m]")
      .def("alt_truth", &LanderGncApp::alt_truth, "[N] truth height above terrain [m]");

  py::class_<LanderNavApp, Application, std::shared_ptr<LanderNavApp>>(
      m, "LanderNavApp",
      "Strapdown-INS multiplicative EKF (MEKF) for a lunar lander on powered descent, fusing a "
      "full IMU, radar altimeter, crater-landmark bearings, and LunaNet (LANS) pseudoranges, with "
      "the IMU biases estimated online. Hosted on a Lander agent as its Application.")
      .def("site_id", &LanderNavApp::site_id, "Landing-site identifier")
      .def("site_name", &LanderNavApp::site_name, "Landing-site name")
      .def("dem_x", &LanderNavApp::dem_x,
           "Terrain DEM grid x coordinates (native projected meters)")
      .def("dem_y", &LanderNavApp::dem_y,
           "Terrain DEM grid y coordinates (native projected meters)")
      .def("dem_elevation", &LanderNavApp::dem_elevation, "Terrain DEM elevation grid [m]")
      .def("crater_enu", &LanderNavApp::crater_enu, "[M x 2] crater (East, North) [m]")
      .def("time_s", &LanderNavApp::time_series, "[N] epoch times [s]")
      .def("pos_err_enu", &LanderNavApp::pos_err_enu,
           "[N x 3] (truth - est) position error in ENU [m]")
      .def("pos_sigma_enu", &LanderNavApp::pos_sigma_enu,
           "[N x 3] 1-sigma position uncertainty in ENU [m]")
      .def("pos_err_norm", &LanderNavApp::pos_err_norm, "[N] 3D position error magnitude [m]")
      .def("vel_err_norm", &LanderNavApp::vel_err_norm, "[N] 3D velocity error magnitude [m/s]")
      .def("clock_bias_err", &LanderNavApp::clock_bias_err, "[N] clock-bias error [s]")
      .def("clock_bias_sigma", &LanderNavApp::clock_bias_sigma, "[N] clock-bias 1-sigma [s]")
      .def("n_visible_sat", &LanderNavApp::n_visible_sat,
           "[N] number of visible LunaNet satellites")
      .def("n_craters", &LanderNavApp::n_craters, "[N] number of tracked crater landmarks")
      .def("accel_bias_err", &LanderNavApp::accel_bias_err,
           "[N x 3] (truth - est) accelerometer bias [m/s^2]")
      .def("accel_bias_sigma", &LanderNavApp::accel_bias_sigma,
           "[N x 3] accelerometer-bias 1-sigma [m/s^2]")
      .def("gyro_bias_err", &LanderNavApp::gyro_bias_err,
           "[N x 3] (truth - est) gyroscope bias [rad/s]")
      .def("gyro_bias_sigma", &LanderNavApp::gyro_bias_sigma,
           "[N x 3] gyroscope-bias 1-sigma [rad/s]")
      .def("att_err_deg", &LanderNavApp::att_err_deg,
           "[N x 3] attitude error (rotation vector) [deg]")
      .def("att_sigma_deg", &LanderNavApp::att_sigma_deg, "[N x 3] attitude 1-sigma [deg]")
      .def("att_quat_est", &LanderNavApp::att_quat_est,
           "[N x 4] estimated body-to-nav attitude quaternion [w,x,y,z] (MEKF)")
      .def("traj_enu_truth", &LanderNavApp::traj_enu_truth, "[N x 3] truth (East, North, Up) [m]")
      .def("traj_enu_est", &LanderNavApp::traj_enu_est, "[N x 3] estimated (East, North, Up) [m]")
      .def("alt_truth", &LanderNavApp::alt_truth, "[N] truth height above terrain [m]")
      .def("alt_est", &LanderNavApp::alt_est, "[N] estimated height above terrain [m]")
      .def("satellite_names", &LanderNavApp::satellite_names,
           "Names of the LunaNet relay satellites");
}

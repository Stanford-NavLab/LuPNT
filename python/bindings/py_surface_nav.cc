/**
 * @file py_surface_nav.cc
 * @brief Python bindings for the agent-based lunar surface-rover navigation application
 *        (`lupnt::SurfaceRoverNavApp`, `lupnt/applications/rover/surface_rover_nav_app.h`).
 *
 * The rover scenario now runs through the generic agent framework:
 *   sim = pnt.Simulation("configs/surface_rover_nav.yaml"); sim.run()
 *   app = sim.get_agent("Rover").get_application()   # -> SurfaceRoverNavApp
 * The app precomputes the driving-arc truth, synthesizes the IMU + LCRNS pseudorange +
 * DEM-constraint measurements from the shared World, runs the strapdown-INS error-state EKF,
 * and exposes the per-epoch truth/estimate/covariance series through the accessors below (all
 * `VecXd`/`MatXd`/`VecXi`, so pybind returns NumPy arrays).
 *
 * `SurfaceRoverNavApp` is registered as an `Application` subclass (polymorphic base bound in
 * `py_simulation.cc`) so `Agent.get_application()` downcasts automatically.
 */
#include <lupnt/lupnt.h>

#include <memory>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitSurfaceNav(py::module& m) {
  py::class_<SurfaceRoverNavApp, Application, std::shared_ptr<SurfaceRoverNavApp>>(
      m, "SurfaceRoverNavApp",
      "Strapdown-INS error-state EKF for a lunar surface rover, fusing a full IMU, LCRNS (LANS) "
      "pseudoranges, and a DEM altitude constraint, with the IMU biases estimated online. Hosted "
      "on a Rover agent as its Application.")
      .def("site_id", &SurfaceRoverNavApp::site_id, "Rover-site identifier")
      .def("site_name", &SurfaceRoverNavApp::site_name, "Rover-site name")
      .def("dem_x", &SurfaceRoverNavApp::dem_x,
           "Terrain DEM grid x coordinates (native projected meters)")
      .def("dem_y", &SurfaceRoverNavApp::dem_y,
           "Terrain DEM grid y coordinates (native projected meters)")
      .def("dem_elevation", &SurfaceRoverNavApp::dem_elevation, "Terrain DEM elevation grid [m]")
      .def("time_s", &SurfaceRoverNavApp::time_series, "[N] epoch times [s]")
      .def("pos_err_enu", &SurfaceRoverNavApp::pos_err_enu,
           "[N x 3] (truth - est) position error in ENU [m]")
      .def("pos_sigma_enu", &SurfaceRoverNavApp::pos_sigma_enu,
           "[N x 3] 1-sigma position uncertainty in ENU [m]")
      .def("pos_err_norm", &SurfaceRoverNavApp::pos_err_norm, "[N] 3D position error magnitude [m]")
      .def("clock_bias_err", &SurfaceRoverNavApp::clock_bias_err, "[N] clock-bias error [s]")
      .def("clock_bias_sigma", &SurfaceRoverNavApp::clock_bias_sigma, "[N] clock-bias 1-sigma [s]")
      .def("n_visible", &SurfaceRoverNavApp::n_visible, "[N] number of visible LCRNS satellites")
      .def("accel_bias_err", &SurfaceRoverNavApp::accel_bias_err,
           "[N x 3] (truth - est) accelerometer bias [m/s^2]")
      .def("accel_bias_sigma", &SurfaceRoverNavApp::accel_bias_sigma,
           "[N x 3] accelerometer-bias 1-sigma [m/s^2]")
      .def("gyro_bias_err", &SurfaceRoverNavApp::gyro_bias_err,
           "[N x 3] (truth - est) gyroscope bias [rad/s]")
      .def("gyro_bias_sigma", &SurfaceRoverNavApp::gyro_bias_sigma,
           "[N x 3] gyroscope-bias 1-sigma [rad/s]")
      .def("att_err_deg", &SurfaceRoverNavApp::att_err_deg,
           "[N x 3] attitude error (rotation vector) [deg]")
      .def("att_sigma_deg", &SurfaceRoverNavApp::att_sigma_deg, "[N x 3] attitude 1-sigma [deg]")
      .def("rover_track_enu_truth", &SurfaceRoverNavApp::rover_track_enu_truth,
           "[N x 2] truth rover track (East, North) [m]")
      .def("rover_track_enu_est", &SurfaceRoverNavApp::rover_track_enu_est,
           "[N x 2] estimated rover track (East, North) [m]")
      .def("rover_alt_truth", &SurfaceRoverNavApp::rover_alt_truth,
           "[N] truth rover Up / elevation [m]")
      .def("satellite_names", &SurfaceRoverNavApp::satellite_names,
           "Names of the LCRNS relay satellites");
}

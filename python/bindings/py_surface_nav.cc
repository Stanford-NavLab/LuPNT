/**
 * @file py_surface_nav.cc
 * @brief Python bindings for `lupnt::RunSurfaceNav`
 *        (`lupnt/simulations/surface_nav/surface_nav_simulation.h`) -- a lunar
 *        surface-rover strapdown-INS navigation simulation fusing a full IMU
 *        (accelerometer + gyroscope, Kalibr noise model), LCRNS (LANS)
 *        pseudoranges, and a NASA-DEM altitude constraint in an error-state EKF
 *        that estimates the IMU biases online.
 *
 * Like `py_isl_odts.cc`, the config/result structs are plain `double`-valued so
 * they bind directly via `def_readwrite`/`def_readonly`; `Vec3d`/`MatXd`/`VecXd`
 * use pybind11's Eigen<->NumPy casters and the `std::vector<...>` members rely on
 * pybind11/stl.h.
 */
#include <lupnt/lupnt.h>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitSurfaceNav(py::module& m) {
  // ---- LcrnsSatConfig -------------------------------------------------------
  py::class_<LcrnsSatConfig>(m, "LcrnsSatConfig")
      .def(py::init<>())
      .def_readwrite("name", &LcrnsSatConfig::name)
      .def_readwrite("r0_m", &LcrnsSatConfig::r0_m,
                     "Initial position [m] in Frame.MOON_CI at start_epoch_utc")
      .def_readwrite("v0_mps", &LcrnsSatConfig::v0_mps,
                     "Initial velocity [m/s] in Frame.MOON_CI at start_epoch_utc");

  // ---- SurfaceNavConfig -----------------------------------------------------
  py::class_<SurfaceNavConfig>(m, "SurfaceNavConfig")
      .def(py::init<>())
      .def_readwrite("seed", &SurfaceNavConfig::seed)
      .def_readwrite("start_epoch_utc", &SurfaceNavConfig::start_epoch_utc)
      .def_readwrite("duration_s", &SurfaceNavConfig::duration_s)
      .def_readwrite("dt_s", &SurfaceNavConfig::dt_s)
      .def_readwrite("site_lat_deg", &SurfaceNavConfig::site_lat_deg)
      .def_readwrite("site_lon_deg", &SurfaceNavConfig::site_lon_deg)
      .def_readwrite("dem_half_width_m", &SurfaceNavConfig::dem_half_width_m)
      .def_readwrite("dem_max_res_m", &SurfaceNavConfig::dem_max_res_m)
      .def_readwrite("rover_start_east_m", &SurfaceNavConfig::rover_start_east_m)
      .def_readwrite("rover_start_north_m", &SurfaceNavConfig::rover_start_north_m)
      .def_readwrite("rover_speed_mps", &SurfaceNavConfig::rover_speed_mps)
      .def_readwrite("rover_heading_deg", &SurfaceNavConfig::rover_heading_deg)
      .def_readwrite("rover_turn_rate_dps", &SurfaceNavConfig::rover_turn_rate_dps)
      .def_readwrite("rover_clock_bias_s", &SurfaceNavConfig::rover_clock_bias_s)
      .def_readwrite("rover_clock_drift_sps", &SurfaceNavConfig::rover_clock_drift_sps)
      // IMU (Kalibr noise model)
      .def_readwrite("accel_noise_density", &SurfaceNavConfig::accel_noise_density,
                     "Accelerometer white-noise density [m/s^2 / sqrt(Hz)]")
      .def_readwrite("accel_bias_rw", &SurfaceNavConfig::accel_bias_rw,
                     "Accelerometer bias random-walk density [m/s^3 / sqrt(Hz)]")
      .def_readwrite("gyro_noise_density", &SurfaceNavConfig::gyro_noise_density,
                     "Gyroscope white-noise density [rad/s / sqrt(Hz)]")
      .def_readwrite("gyro_bias_rw", &SurfaceNavConfig::gyro_bias_rw,
                     "Gyroscope bias random-walk density [rad/s^2 / sqrt(Hz)]")
      .def_readwrite("accel_bias0", &SurfaceNavConfig::accel_bias0,
                     "Truth initial accel bias 1-sigma per axis [m/s^2]")
      .def_readwrite("gyro_bias0", &SurfaceNavConfig::gyro_bias0,
                     "Truth initial gyro bias 1-sigma per axis [rad/s]")
      // LCRNS
      .def_readwrite("satellites", &SurfaceNavConfig::satellites)
      .def_readwrite("elevation_mask_deg", &SurfaceNavConfig::elevation_mask_deg)
      .def_readwrite("pseudorange_sigma_m", &SurfaceNavConfig::pseudorange_sigma_m)
      .def_readwrite("sise_m", &SurfaceNavConfig::sise_m)
      // Filter init
      .def_readwrite("init_pos_sigma_m", &SurfaceNavConfig::init_pos_sigma_m)
      .def_readwrite("init_vel_sigma_mps", &SurfaceNavConfig::init_vel_sigma_mps)
      .def_readwrite("init_att_sigma_deg", &SurfaceNavConfig::init_att_sigma_deg)
      .def_readwrite("init_accel_bias_sigma", &SurfaceNavConfig::init_accel_bias_sigma)
      .def_readwrite("init_gyro_bias_sigma", &SurfaceNavConfig::init_gyro_bias_sigma)
      .def_readwrite("init_clock_bias_sigma_s", &SurfaceNavConfig::init_clock_bias_sigma_s)
      .def_readwrite("init_clock_drift_sigma_sps", &SurfaceNavConfig::init_clock_drift_sigma_sps)
      .def_readwrite("dem_sigma_m", &SurfaceNavConfig::dem_sigma_m)
      .def_readwrite("filter_bias_rw_scale", &SurfaceNavConfig::filter_bias_rw_scale,
                     "Filter IMU bias random-walk inflation vs truth (>= 1) for EKF consistency")
      .def_readwrite("enable_dem_constraint", &SurfaceNavConfig::enable_dem_constraint);

  // ---- SurfaceNavResults ----------------------------------------------------
  py::class_<SurfaceNavResults>(m, "SurfaceNavResults")
      .def_readonly("site_id", &SurfaceNavResults::site_id)
      .def_readonly("site_name", &SurfaceNavResults::site_name)
      .def_readonly("site_lat_deg", &SurfaceNavResults::site_lat_deg)
      .def_readonly("site_lon_deg", &SurfaceNavResults::site_lon_deg)
      .def_readonly("dem_x", &SurfaceNavResults::dem_x)
      .def_readonly("dem_y", &SurfaceNavResults::dem_y)
      .def_readonly("dem_elevation", &SurfaceNavResults::dem_elevation)
      .def_readonly("dem_center_x", &SurfaceNavResults::dem_center_x)
      .def_readonly("dem_center_y", &SurfaceNavResults::dem_center_y)
      .def_readonly("time_s", &SurfaceNavResults::time_s)
      .def_readonly("pos_err_enu", &SurfaceNavResults::pos_err_enu,
                    "[N x 3] (truth - est) position error in ENU [m]")
      .def_readonly("pos_sigma_enu", &SurfaceNavResults::pos_sigma_enu,
                    "[N x 3] 1-sigma position uncertainty in ENU [m]")
      .def_readonly("pos_err_norm", &SurfaceNavResults::pos_err_norm)
      .def_readonly("clock_bias_err", &SurfaceNavResults::clock_bias_err)
      .def_readonly("clock_bias_sigma", &SurfaceNavResults::clock_bias_sigma)
      .def_readonly("n_visible", &SurfaceNavResults::n_visible)
      .def_readonly("accel_bias_err", &SurfaceNavResults::accel_bias_err,
                    "[N x 3] (truth - est) accelerometer bias [m/s^2]")
      .def_readonly("accel_bias_sigma", &SurfaceNavResults::accel_bias_sigma,
                    "[N x 3] accelerometer-bias 1-sigma [m/s^2]")
      .def_readonly("gyro_bias_err", &SurfaceNavResults::gyro_bias_err,
                    "[N x 3] (truth - est) gyroscope bias [rad/s]")
      .def_readonly("gyro_bias_sigma", &SurfaceNavResults::gyro_bias_sigma,
                    "[N x 3] gyroscope-bias 1-sigma [rad/s]")
      .def_readonly("att_err_deg", &SurfaceNavResults::att_err_deg,
                    "[N x 3] attitude error (rotation vector) [deg]")
      .def_readonly("att_sigma_deg", &SurfaceNavResults::att_sigma_deg,
                    "[N x 3] attitude 1-sigma [deg]")
      .def_readonly("rover_track_enu_truth", &SurfaceNavResults::rover_track_enu_truth)
      .def_readonly("rover_track_enu_est", &SurfaceNavResults::rover_track_enu_est)
      .def_readonly("rover_alt_truth", &SurfaceNavResults::rover_alt_truth)
      .def_readonly("satellite_names", &SurfaceNavResults::satellite_names);

  m.def("run_surface_nav", &RunSurfaceNav, py::arg("config"),
        "Run the surface-rover IMU + LCRNS + DEM strapdown-INS navigation simulation.");
}

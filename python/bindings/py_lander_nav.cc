/**
 * @file py_lander_nav.cc
 * @brief Python bindings for `lupnt::RunLanderNav`
 *        (`lupnt/simulations/lander_nav/lander_nav_simulation.h`) -- a lunar
 *        lander powered-descent navigation simulation. A `Lander` agent hosts a
 *        `LanderNavApp` error-state INS EKF fusing a full IMU (accelerometer +
 *        gyroscope, Kalibr noise model), a nadir radar altimeter, crater-landmark
 *        bearings (terrain-relative navigation), and LunaNet (LANS) pseudoranges,
 *        estimating the IMU biases online.
 *
 * Reuses `LcrnsSatConfig` (already bound by `InitSurfaceNav`). The config/result
 * structs are plain `double`/`int`-valued so they bind directly via
 * `def_readwrite`/`def_readonly`; `Vec3d`/`MatXd`/`VecXd` use pybind11's
 * Eigen<->NumPy casters and `std::vector<...>` members rely on pybind11/stl.h.
 */
#include <lupnt/lupnt.h>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitLanderNav(py::module& m) {
  // ---- LanderNavConfig ------------------------------------------------------
  py::class_<LanderNavConfig>(m, "LanderNavConfig")
      .def(py::init<>())
      .def_readwrite("seed", &LanderNavConfig::seed)
      .def_readwrite("start_epoch_utc", &LanderNavConfig::start_epoch_utc)
      .def_readwrite("duration_s", &LanderNavConfig::duration_s)
      .def_readwrite("dt_s", &LanderNavConfig::dt_s)
      // Site / DEM
      .def_readwrite("site_lat_deg", &LanderNavConfig::site_lat_deg)
      .def_readwrite("site_lon_deg", &LanderNavConfig::site_lon_deg)
      .def_readwrite("dem_half_width_m", &LanderNavConfig::dem_half_width_m)
      .def_readwrite("dem_max_res_m", &LanderNavConfig::dem_max_res_m)
      // Descent truth
      .def_readwrite("descent_start_east_m", &LanderNavConfig::descent_start_east_m,
                     "Start East offset (downrange) from tile center [m]")
      .def_readwrite("descent_start_north_m", &LanderNavConfig::descent_start_north_m,
                     "Start North offset (crossrange) from tile center [m]")
      .def_readwrite("descent_end_east_m", &LanderNavConfig::descent_end_east_m,
                     "Landing East offset [m]")
      .def_readwrite("descent_end_north_m", &LanderNavConfig::descent_end_north_m,
                     "Landing North offset [m]")
      .def_readwrite("descent_start_alt_m", &LanderNavConfig::descent_start_alt_m,
                     "Height above terrain at t0 [m]")
      .def_readwrite("descent_end_alt_m", &LanderNavConfig::descent_end_alt_m,
                     "Height above terrain at touchdown/hover [m]")
      .def_readwrite("descent_heading_deg", &LanderNavConfig::descent_heading_deg,
                     "Body-x heading in ENU (0=E, 90=N) [deg]")
      .def_readwrite("lander_clock_bias_s", &LanderNavConfig::lander_clock_bias_s)
      .def_readwrite("lander_clock_drift_sps", &LanderNavConfig::lander_clock_drift_sps)
      .def_readwrite("ref_traj_enu", &LanderNavConfig::ref_traj_enu,
                     "Optional [N x 3] ENU reference (truth) trajectory [m] about the tile center "
                     "(U = height above the site datum). Non-empty overrides the built-in "
                     "smoothstep descent and sets N; see pylupnt.lander_guidance.")
      // IMU (Kalibr noise model)
      .def_readwrite("accel_noise_density", &LanderNavConfig::accel_noise_density,
                     "Accelerometer white-noise density [m/s^2 / sqrt(Hz)]")
      .def_readwrite("accel_bias_rw", &LanderNavConfig::accel_bias_rw,
                     "Accelerometer bias random-walk density [m/s^3 / sqrt(Hz)]")
      .def_readwrite("gyro_noise_density", &LanderNavConfig::gyro_noise_density,
                     "Gyroscope white-noise density [rad/s / sqrt(Hz)]")
      .def_readwrite("gyro_bias_rw", &LanderNavConfig::gyro_bias_rw,
                     "Gyroscope bias random-walk density [rad/s^2 / sqrt(Hz)]")
      .def_readwrite("accel_bias0", &LanderNavConfig::accel_bias0,
                     "Truth initial accel bias 1-sigma per axis [m/s^2]")
      .def_readwrite("gyro_bias0", &LanderNavConfig::gyro_bias0,
                     "Truth initial gyro bias 1-sigma per axis [rad/s]")
      // Radar altimeter
      .def_readwrite("enable_altimeter", &LanderNavConfig::enable_altimeter)
      .def_readwrite("altimeter_sigma_m", &LanderNavConfig::altimeter_sigma_m,
                     "Altimeter noise 1-sigma [m]")
      .def_readwrite("altimeter_max_range_m", &LanderNavConfig::altimeter_max_range_m,
                     "Max altitude with a valid altimeter return [m]")
      // Crater-landmark camera
      .def_readwrite("enable_craters", &LanderNavConfig::enable_craters)
      .def_readwrite("n_craters", &LanderNavConfig::n_craters,
                     "Number of synthetic craters in the map")
      .def_readwrite("crater_field_radius_m", &LanderNavConfig::crater_field_radius_m,
                     "Craters scattered within this radius of the site [m]")
      .def_readwrite("camera_fov_deg", &LanderNavConfig::camera_fov_deg,
                     "Camera half-cone about nadir [deg]")
      .def_readwrite("max_craters_per_epoch", &LanderNavConfig::max_craters_per_epoch)
      .def_readwrite("crater_sigma_arcsec", &LanderNavConfig::crater_sigma_arcsec,
                     "Crater-bearing noise 1-sigma [arcsec]")
      // LunaNet (LCRNS)
      .def_readwrite("enable_lunanet", &LanderNavConfig::enable_lunanet)
      .def_readwrite("satellites", &LanderNavConfig::satellites)
      .def_readwrite("elevation_mask_deg", &LanderNavConfig::elevation_mask_deg)
      .def_readwrite("pseudorange_sigma_m", &LanderNavConfig::pseudorange_sigma_m)
      .def_readwrite("sise_m", &LanderNavConfig::sise_m)
      // Filter init & tuning
      .def_readwrite("init_pos_sigma_m", &LanderNavConfig::init_pos_sigma_m)
      .def_readwrite("init_vel_sigma_mps", &LanderNavConfig::init_vel_sigma_mps)
      .def_readwrite("init_att_sigma_deg", &LanderNavConfig::init_att_sigma_deg)
      .def_readwrite("init_accel_bias_sigma", &LanderNavConfig::init_accel_bias_sigma)
      .def_readwrite("init_gyro_bias_sigma", &LanderNavConfig::init_gyro_bias_sigma)
      .def_readwrite("init_clock_bias_sigma_s", &LanderNavConfig::init_clock_bias_sigma_s)
      .def_readwrite("init_clock_drift_sigma_sps", &LanderNavConfig::init_clock_drift_sigma_sps)
      .def_readwrite("filter_bias_rw_scale", &LanderNavConfig::filter_bias_rw_scale,
                     "Filter IMU bias random-walk inflation vs truth (>= 1) for EKF consistency");

  // ---- LanderNavResults -----------------------------------------------------
  py::class_<LanderNavResults>(m, "LanderNavResults")
      .def_readonly("site_id", &LanderNavResults::site_id)
      .def_readonly("site_name", &LanderNavResults::site_name)
      .def_readonly("site_lat_deg", &LanderNavResults::site_lat_deg)
      .def_readonly("site_lon_deg", &LanderNavResults::site_lon_deg)
      .def_readonly("dem_x", &LanderNavResults::dem_x)
      .def_readonly("dem_y", &LanderNavResults::dem_y)
      .def_readonly("dem_elevation", &LanderNavResults::dem_elevation)
      .def_readonly("dem_center_x", &LanderNavResults::dem_center_x)
      .def_readonly("dem_center_y", &LanderNavResults::dem_center_y)
      .def_readonly("crater_enu", &LanderNavResults::crater_enu,
                    "[M x 2] crater (East, North) about tile center [m]")
      .def_readonly("time_s", &LanderNavResults::time_s)
      .def_readonly("pos_err_enu", &LanderNavResults::pos_err_enu,
                    "[N x 3] (truth - est) position error in ENU [m]")
      .def_readonly("pos_sigma_enu", &LanderNavResults::pos_sigma_enu,
                    "[N x 3] 1-sigma position uncertainty in ENU [m]")
      .def_readonly("pos_err_norm", &LanderNavResults::pos_err_norm)
      .def_readonly("vel_err_norm", &LanderNavResults::vel_err_norm)
      .def_readonly("clock_bias_err", &LanderNavResults::clock_bias_err)
      .def_readonly("clock_bias_sigma", &LanderNavResults::clock_bias_sigma)
      .def_readonly("n_visible_sat", &LanderNavResults::n_visible_sat,
                    "[N] number of visible LunaNet satellites")
      .def_readonly("n_craters", &LanderNavResults::n_craters,
                    "[N] number of tracked crater landmarks")
      .def_readonly("accel_bias_err", &LanderNavResults::accel_bias_err,
                    "[N x 3] (truth - est) accelerometer bias [m/s^2]")
      .def_readonly("accel_bias_sigma", &LanderNavResults::accel_bias_sigma,
                    "[N x 3] accelerometer-bias 1-sigma [m/s^2]")
      .def_readonly("gyro_bias_err", &LanderNavResults::gyro_bias_err,
                    "[N x 3] (truth - est) gyroscope bias [rad/s]")
      .def_readonly("gyro_bias_sigma", &LanderNavResults::gyro_bias_sigma,
                    "[N x 3] gyroscope-bias 1-sigma [rad/s]")
      .def_readonly("att_err_deg", &LanderNavResults::att_err_deg,
                    "[N x 3] attitude error (rotation vector) [deg]")
      .def_readonly("att_sigma_deg", &LanderNavResults::att_sigma_deg,
                    "[N x 3] attitude 1-sigma [deg]")
      .def_readonly("att_quat_est", &LanderNavResults::att_quat_est,
                    "[N x 4] estimated body-to-nav attitude quaternion [w,x,y,z] (MEKF)")
      .def_readonly("traj_enu_truth", &LanderNavResults::traj_enu_truth,
                    "[N x 3] truth (East, North, Up) [m]")
      .def_readonly("traj_enu_est", &LanderNavResults::traj_enu_est,
                    "[N x 3] estimated (East, North, Up) [m]")
      .def_readonly("alt_truth", &LanderNavResults::alt_truth, "[N] truth height above terrain [m]")
      .def_readonly("alt_est", &LanderNavResults::alt_est, "[N] estimated height above terrain [m]")
      .def_readonly("satellite_names", &LanderNavResults::satellite_names);

  m.def("run_lander_nav", &RunLanderNav, py::arg("config"),
        "Run the lunar-lander IMU + altimeter + crater-bearing + LunaNet descent-navigation "
        "simulation (a Lander agent hosting a LanderNavApp error-state INS EKF).");
}

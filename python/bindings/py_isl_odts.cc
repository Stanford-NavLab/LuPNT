/**
 * @file py_isl_odts.cc
 * @brief Python bindings for `lupnt::IslOdtsSimulation`
 *        (`lupnt/simulations/isl_odts/isl_odts_simulation.h`) -- an onboard
 *        inter-satellite-link (ISL) orbit determination and timing system
 *        (ODTS) simulation for a constellation of N satellites, where one hub
 *        satellite runs a Schmidt Extended Kalman Filter over crosslinks to
 *        each of the other N-1.
 *
 * Config/result structs are plain `double`-valued (no `lupnt::Real`), so they
 * bind directly via `def_readwrite`/`def_readonly` with pybind11's built-in
 * Eigen<->NumPy conversion, following the `GnssReceiverParams` pattern in
 * `py_gnss.cc`. `IslOdtsConfig.satellites` (std::vector<IslOdtsSatelliteConfig>)
 * and `IslOdtsResults.satellite_names`/`.truth_states` (std::vector<...>) rely on
 * pybind11/stl.h for Python list<->std::vector conversion.
 */
#include <lupnt/lupnt.h>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitIslOdts(py::module& m) {
  // ---- IslOdtsSatelliteConfig -----------------------------------------------

  py::class_<IslOdtsSatelliteConfig>(m, "IslOdtsSatelliteConfig")
      .def(py::init<>())
      .def_readwrite("name", &IslOdtsSatelliteConfig::name)
      .def_readwrite("r0_m", &IslOdtsSatelliteConfig::r0_m,
                     "Initial position [m] in Frame.MOON_CI at start_epoch_utc")
      .def_readwrite("v0_mps", &IslOdtsSatelliteConfig::v0_mps,
                     "Initial velocity [m/s] in Frame.MOON_CI at start_epoch_utc")
      .def_readwrite("clock_bias_s", &IslOdtsSatelliteConfig::clock_bias_s)
      .def_readwrite("clock_drift_sps", &IslOdtsSatelliteConfig::clock_drift_sps);

  // ---- IslLinkBudgetConfig ---------------------------------------------------

  py::class_<IslLinkBudgetConfig>(m, "IslLinkBudgetConfig")
      .def(py::init<>())
      .def_readwrite("enabled", &IslLinkBudgetConfig::enabled)
      .def_readwrite("tx_power_dbw", &IslLinkBudgetConfig::tx_power_dbw)
      .def_readwrite("tx_gain_dbi", &IslLinkBudgetConfig::tx_gain_dbi)
      .def_readwrite("rx_gain_dbi", &IslLinkBudgetConfig::rx_gain_dbi)
      .def_readwrite("frequency_hz", &IslLinkBudgetConfig::frequency_hz)
      .def_readwrite("system_noise_temp_k", &IslLinkBudgetConfig::system_noise_temp_k);

  // ---- IslSurfaceStationConfig -----------------------------------------------

  py::class_<IslSurfaceStationConfig>(m, "IslSurfaceStationConfig")
      .def(py::init<>())
      .def_readwrite("enabled", &IslSurfaceStationConfig::enabled,
                     "If true, a lunar surface station serves the satellites one at a time "
                     "(round-robin) with a one-way pseudorange (absolute position + clock aiding)")
      .def_readwrite("name", &IslSurfaceStationConfig::name)
      .def_readwrite("latitude_deg", &IslSurfaceStationConfig::latitude_deg,
                     "Station geodetic latitude [deg], Moon principal-axis frame")
      .def_readwrite("longitude_deg", &IslSurfaceStationConfig::longitude_deg,
                     "Station geodetic longitude [deg]")
      .def_readwrite("altitude_m", &IslSurfaceStationConfig::altitude_m,
                     "Station altitude above the reference sphere [m]")
      .def_readwrite("pseudorange_sigma_m", &IslSurfaceStationConfig::pseudorange_sigma_m,
                     "One-way pseudorange measurement noise 1-sigma [m]")
      .def_readwrite("elevation_mask_deg", &IslSurfaceStationConfig::elevation_mask_deg,
                     "A satellite is served only above this topocentric elevation [deg]");

  // ---- IslOdtsConfig ----------------------------------------------------------

  py::class_<IslOdtsConfig>(m, "IslOdtsConfig")
      .def(py::init<>())
      .def_readwrite("seed", &IslOdtsConfig::seed)
      .def_readwrite("start_epoch_utc", &IslOdtsConfig::start_epoch_utc)
      .def_readwrite("duration_s", &IslOdtsConfig::duration_s)
      .def_readwrite("dt_s", &IslOdtsConfig::dt_s)
      .def_readwrite("integration_step_s", &IslOdtsConfig::integration_step_s)
      .def_readwrite("satellites", &IslOdtsConfig::satellites,
                     "List of >=2 IslOdtsSatelliteConfig, fully cross-linked; every satellite "
                     "runs its own onboard SchmidtEKF in parallel. satellites[0] is only the "
                     "reference used for the crosslink-geometry reporting arrays")
      .def_readwrite("link_budget", &IslOdtsConfig::link_budget)
      .def_readwrite("surface_station", &IslOdtsConfig::surface_station,
                     "Rotating lunar surface-station one-way pseudorange aiding (disabled by "
                     "default); see IslSurfaceStationConfig")
      .def_readwrite("consider_exchange_interval_s", &IslOdtsConfig::consider_exchange_interval_s,
                     "Interval [s] at which the parallel filters exchange own estimate "
                     "(mean + covariance) to refresh each other's consider blocks; <=0 disables")
      .def_readwrite("moon_gravity_degree_truth", &IslOdtsConfig::moon_gravity_degree_truth)
      .def_readwrite("moon_gravity_order_truth", &IslOdtsConfig::moon_gravity_order_truth)
      .def_readwrite("moon_gravity_degree_filter", &IslOdtsConfig::moon_gravity_degree_filter)
      .def_readwrite("moon_gravity_order_filter", &IslOdtsConfig::moon_gravity_order_filter)
      .def_readwrite("include_earth", &IslOdtsConfig::include_earth)
      .def_readwrite("include_sun", &IslOdtsConfig::include_sun)
      .def_readwrite("use_relativity", &IslOdtsConfig::use_relativity)
      .def_readwrite("range_sigma_m", &IslOdtsConfig::range_sigma_m)
      .def_readwrite("range_rate_sigma_mps", &IslOdtsConfig::range_rate_sigma_mps)
      .def_readwrite("initial_position_sigma_m", &IslOdtsConfig::initial_position_sigma_m)
      .def_readwrite("initial_velocity_sigma_mps", &IslOdtsConfig::initial_velocity_sigma_mps)
      .def_readwrite("initial_clock_bias_sigma_s", &IslOdtsConfig::initial_clock_bias_sigma_s)
      .def_readwrite("initial_clock_drift_sigma_sps", &IslOdtsConfig::initial_clock_drift_sigma_sps)
      .def_readwrite("consider_position_sigma_m", &IslOdtsConfig::consider_position_sigma_m)
      .def_readwrite("consider_velocity_sigma_mps", &IslOdtsConfig::consider_velocity_sigma_mps)
      .def_readwrite("consider_clock_bias_sigma_s", &IslOdtsConfig::consider_clock_bias_sigma_s)
      .def_readwrite("consider_clock_drift_sigma_sps",
                     &IslOdtsConfig::consider_clock_drift_sigma_sps)
      .def_readwrite("process_accel_sigma_mps2", &IslOdtsConfig::process_accel_sigma_mps2);

  // ---- IslOdtsResults ----------------------------------------------------------
  // Time series, row k <-> t_s[k]. truth_states[j] columns:
  // [r_x,r_y,r_z,v_x,v_y,v_z,clock_bias_s,clock_drift_sps]. est[j]/cov_diag[j] are
  // satellite j's filter, stacked in GLOBAL sat index order (block g = filter j's
  // estimate of satellite g), so est[j][:, 8*j:8*j+8] is filter j's own state. Per-link
  // matrices (range_* / cn0_dbhz) have one column per crosslink from satellites[0] to
  // satellites[i+1] (n_links = n_sat - 1). range_resid_m[j] holds filter j's crosslink
  // residuals to its neighbors (columns = global indices != j, ascending).

  py::class_<IslOdtsResults>(m, "IslOdtsResults")
      .def(py::init<>())
      .def_readonly("t_s", &IslOdtsResults::t_s)
      .def_readonly("satellite_names", &IslOdtsResults::satellite_names)
      .def_readonly("truth_states", &IslOdtsResults::truth_states)
      .def_readonly("est", &IslOdtsResults::est,
                    "List [n_sat] of [N x 8*n_sat]; est[j] = satellite j's filter estimate "
                    "in global sat order (own block at columns 8*j:8*j+8)")
      .def_readonly("cov_diag", &IslOdtsResults::cov_diag,
                    "List [n_sat] of [N x 8*n_sat] covariance diagonals, same layout as est")
      .def_readonly("range_true_m", &IslOdtsResults::range_true_m)
      .def_readonly("range_rate_true_mps", &IslOdtsResults::range_rate_true_mps)
      .def_readonly("range_obs_m", &IslOdtsResults::range_obs_m)
      .def_readonly("range_rate_obs_mps", &IslOdtsResults::range_rate_obs_mps)
      .def_readonly("range_resid_m", &IslOdtsResults::range_resid_m,
                    "List [n_sat] of [N x n_links]; range_resid_m[j] = filter j's pre-fit "
                    "crosslink range residuals to its neighbors")
      .def_readonly("cn0_dbhz", &IslOdtsResults::cn0_dbhz)
      .def_readonly("served_sat_idx", &IslOdtsResults::served_sat_idx,
                    "[N] global index of the satellite served by the station (-1 if none)")
      .def_readonly("station_pr_true_m", &IslOdtsResults::station_pr_true_m,
                    "[N] truth station pseudorange, NaN when the station is idle")
      .def_readonly("station_pr_obs_m", &IslOdtsResults::station_pr_obs_m,
                    "[N] noisy station pseudorange observation, NaN when idle")
      .def_readonly("station_pr_resid_m", &IslOdtsResults::station_pr_resid_m,
                    "[N] served filter's pre-fit station pseudorange residual, NaN when idle")
      .def_readonly("station_pos_mci", &IslOdtsResults::station_pos_mci,
                    "[N x 3] station inertial position [m], Frame.MOON_CI");

  // ---- IslOdtsSimulation -------------------------------------------------------

  py::class_<IslOdtsSimulation>(m, "IslOdtsSimulation")
      .def(py::init<IslOdtsConfig>(), py::arg("config"))
      .def("setup", &IslOdtsSimulation::Setup)
      .def("run", &IslOdtsSimulation::Run)
      .def("get_config", &IslOdtsSimulation::GetConfig, py::return_value_policy::reference_internal)
      .def("get_results", &IslOdtsSimulation::GetResults,
           py::return_value_policy::reference_internal);
}

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
#include <lupnt/applications/lunar_sat_odts/ground_odts_app.h>
#include <lupnt/applications/lunar_sat_odts/satellite_odts_app.h>
#include <lupnt/lupnt.h>

#include <memory>

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
      .def_readwrite("exchange_use_covariance_intersection",
                     &IslOdtsConfig::exchange_use_covariance_intersection,
                     "If True (default), fuse each broadcast neighbor estimate into the "
                     "consider block with Covariance Intersection (consistent under the "
                     "unknown inter-filter correlation); if False, overwrite the block and "
                     "zero its cross-covariance (optimistic, kept for comparison)")
      .def_readwrite("exchange_ci_weight", &IslOdtsConfig::exchange_ci_weight,
                     "Fixed CI weight w in [0,1] on the local consider block (1-w on the "
                     "broadcast); negative (default) selects w per block by minimizing the "
                     "fused-block covariance trace")
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
      .def_readwrite("process_accel_sigma_mps2", &IslOdtsConfig::process_accel_sigma_mps2)
      .def_readwrite("surface_stations", &IslOdtsConfig::surface_stations,
                     "List of surface-station beacons; if non-empty it overrides the single "
                     "surface_station (each station exchanges a one-way pseudorange with every "
                     "satellite above its elevation mask)")
      .def_readwrite("enable_two_way_time_transfer", &IslOdtsConfig::enable_two_way_time_transfer,
                     "If True, each crosslink also measures the range-equivalent clock-bias "
                     "difference C*(b_i - b_j), making the constellation's relative clocks "
                     "observable (two-way ranging alone is clock-free)")
      .def_readwrite("time_transfer_sigma_m", &IslOdtsConfig::time_transfer_sigma_m,
                     "Two-way time-transfer noise 1-sigma [m] (comparable to the ranging noise)")
      .def_readwrite("enable_two_way_frequency_transfer",
                     &IslOdtsConfig::enable_two_way_frequency_transfer,
                     "If True, each crosslink also measures the clock-DRIFT difference "
                     "C*(d_i - d_j) (rate companion to the time transfer), making relative "
                     "clock drift observable")
      .def_readwrite("frequency_transfer_sigma_mps", &IslOdtsConfig::frequency_transfer_sigma_mps,
                     "Two-way frequency-transfer noise 1-sigma [m/s]")
      .def_readwrite("enable_station_doppler", &IslOdtsConfig::enable_station_doppler,
                     "If True, the surface-station beacon links also provide one-way Doppler "
                     "(pseudorange-rate), adding velocity + clock-drift observability")
      .def_readwrite("station_doppler_sigma_mps", &IslOdtsConfig::station_doppler_sigma_mps,
                     "One-way station Doppler noise 1-sigma [m/s]")
      .def_readwrite("enable_centralized_ground_filter",
                     &IslOdtsConfig::enable_centralized_ground_filter,
                     "If True, also run a centralized EKF on the ground that estimates every "
                     "satellite's orbit+clock from the station pseudoranges only (no ISL)")
      .def_readwrite("central_process_accel_sigma_mps2",
                     &IslOdtsConfig::central_process_accel_sigma_mps2,
                     "Process-noise sigma [m/s^2] for the centralized ground filter (usually "
                     "much larger than the ISL filters' value; range-only ground tracking is "
                     "weakly observable and diverges if over-confident)")
      .def_readwrite("central_outlier_threshold", &IslOdtsConfig::central_outlier_threshold,
                     "Normalized-residual outlier-rejection threshold [sigma] for the "
                     "centralized filter (rejects large linearized residuals from fast "
                     "perilune passes to keep the EKF stable; large value disables)");

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
      .def_readonly("cov_own_full", &IslOdtsResults::cov_own_full,
                    "List [n_sat] of [N x 64]; row k is satellite j's own 8x8 covariance "
                    "row-major flattened (reshape to (N,8,8)) -- for NEES consistency checks")
      .def_readonly("range_true_m", &IslOdtsResults::range_true_m)
      .def_readonly("range_rate_true_mps", &IslOdtsResults::range_rate_true_mps)
      .def_readonly("range_obs_m", &IslOdtsResults::range_obs_m)
      .def_readonly("range_rate_obs_mps", &IslOdtsResults::range_rate_obs_mps)
      .def_readonly("range_resid_m", &IslOdtsResults::range_resid_m,
                    "List [n_sat] of [N x n_links]; range_resid_m[j] = filter j's pre-fit "
                    "crosslink range residuals to its neighbors")
      .def_readonly("cn0_dbhz", &IslOdtsResults::cn0_dbhz)
      .def_readonly("time_transfer_true_m", &IslOdtsResults::time_transfer_true_m,
                    "[N x n_links] two-way time-transfer truth C*(b_0 - b_{i+1}); NaN if disabled")
      .def_readonly("time_transfer_obs_m", &IslOdtsResults::time_transfer_obs_m,
                    "[N x n_links] noisy two-way time-transfer observation; NaN if disabled")
      .def_readonly("station_pos_mci", &IslOdtsResults::station_pos_mci,
                    "List [n_station] of [N x 3] station inertial positions [m], Frame.MOON_CI")
      .def_readonly("station_visible", &IslOdtsResults::station_visible,
                    "[N x n_sat] number of stations that see each satellite each epoch")
      .def_readonly("station_pr_true_m", &IslOdtsResults::station_pr_true_m,
                    "[N x n_sat] truth station pseudorange per satellite (first visible "
                    "station), NaN where no station sees it")
      .def_readonly("station_pr_obs_m", &IslOdtsResults::station_pr_obs_m,
                    "[N x n_sat] noisy station pseudorange per satellite, NaN if not visible")
      .def_readonly("station_pr_resid_m", &IslOdtsResults::station_pr_resid_m,
                    "[N x n_sat] onboard filter's pre-fit station pseudorange residual per "
                    "satellite, NaN if not visible")
      .def_readonly("est_central", &IslOdtsResults::est_central,
                    "[N x 8*n_sat] centralized ground filter estimate, [r,v,cb,cd] per "
                    "satellite in global order (zeros if the centralized filter is disabled)")
      .def_readonly("cov_central_full", &IslOdtsResults::cov_central_full,
                    "List [n_sat] of [N x 64]; centralized filter's own 8x8 covariance per "
                    "satellite, row-major flattened (reshape to (N,8,8))");

  // ---- IslOdtsCoordinatorApp ---------------------------------------------------
  // The agent-based coordinator hosted on an `IslOdtsManager` agent. Retrieve it from a
  // `pnt.Simulation` via `sim.get_agent("IslManager").get_application()` (downcasts here),
  // then read its time-series `IslOdtsResults`.
  py::class_<IslOdtsCoordinatorApp, Application, std::shared_ptr<IslOdtsCoordinatorApp>>(
      m, "IslOdtsCoordinatorApp")
      .def("get_config", &IslOdtsCoordinatorApp::GetConfig,
           py::return_value_policy::reference_internal)
      .def("get_results", &IslOdtsCoordinatorApp::GetResults,
           py::return_value_policy::reference_internal,
           "Time-series IslOdtsResults, identical layout to IslOdtsSimulation.get_results()");

  // ---- Distributed variant: per-satellite onboard app + ground-segment app ----
  py::class_<SatelliteOdtsApp, Application, std::shared_ptr<SatelliteOdtsApp>>(m,
                                                                               "SatelliteOdtsApp")
      .def("sat_name", &SatelliteOdtsApp::SatName)
      .def("neighbor_names", &SatelliteOdtsApp::NeighborNames)
      .def("time_grid", &SatelliteOdtsApp::TimeGrid)
      .def("truth_state", &SatelliteOdtsApp::TruthState, "[N x 8] own truth [r,v,cb,cd]")
      .def("own_estimate", &SatelliteOdtsApp::OwnEstimate, "[N x 8] own onboard estimate")
      .def("own_cov_diag", &SatelliteOdtsApp::OwnCovDiag, "[N x 8] own covariance diagonal")
      .def("own_cov_full", &SatelliteOdtsApp::OwnCovFull, "[N x 64] own 8x8 covariance (row-major)")
      .def("num_anchors_total", &SatelliteOdtsApp::NumAnchorsTotal);

  py::class_<GroundOdtsApp, Application, std::shared_ptr<GroundOdtsApp>>(m, "GroundOdtsApp")
      .def("satellite_names", &GroundOdtsApp::SatelliteNames)
      .def("time_grid", &GroundOdtsApp::TimeGrid)
      .def("est_central", &GroundOdtsApp::EstCentral, "[N x 8*n_sat] centralized ground estimate")
      .def("cov_central_full", &GroundOdtsApp::CovCentralFull,
           "n_sat x [N x 64] per-satellite 8x8 covariance (row-major)");
}

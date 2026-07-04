/**
 * @file py_isl_odts.cc
 * @brief Python bindings for `lupnt::IslOdtsSimulation`
 *        (`lupnt/simulations/IslOdts/isl_odts_simulation.h`) -- an onboard
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

  // ---- IslOdtsConfig ----------------------------------------------------------

  py::class_<IslOdtsConfig>(m, "IslOdtsConfig")
      .def(py::init<>())
      .def_readwrite("seed", &IslOdtsConfig::seed)
      .def_readwrite("start_epoch_utc", &IslOdtsConfig::start_epoch_utc)
      .def_readwrite("duration_s", &IslOdtsConfig::duration_s)
      .def_readwrite("dt_s", &IslOdtsConfig::dt_s)
      .def_readwrite("integration_step_s", &IslOdtsConfig::integration_step_s)
      .def_readwrite("satellites", &IslOdtsConfig::satellites,
                      "List of >=2 IslOdtsSatelliteConfig; satellites[0] is the hub satellite "
                      "whose onboard SchmidtEKF is run, satellites[1:] are the satellites it "
                      "links to (one crosslink each)")
      .def_readwrite("link_budget", &IslOdtsConfig::link_budget)
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
  // [r_x,r_y,r_z,v_x,v_y,v_z,clock_bias_s,clock_drift_sps]. est/cov_diag stack
  // [own(8), consider_1(8), ..., consider_{n_sat-1}(8)]. Per-link matrices (range_*,
  // cn0_dbhz) have one column per crosslink (n_links = n_sat - 1), in the same
  // order as satellite_names[1:].

  py::class_<IslOdtsResults>(m, "IslOdtsResults")
      .def(py::init<>())
      .def_readonly("t_s", &IslOdtsResults::t_s)
      .def_readonly("satellite_names", &IslOdtsResults::satellite_names)
      .def_readonly("truth_states", &IslOdtsResults::truth_states)
      .def_readonly("est", &IslOdtsResults::est)
      .def_readonly("cov_diag", &IslOdtsResults::cov_diag)
      .def_readonly("range_true_m", &IslOdtsResults::range_true_m)
      .def_readonly("range_rate_true_mps", &IslOdtsResults::range_rate_true_mps)
      .def_readonly("range_obs_m", &IslOdtsResults::range_obs_m)
      .def_readonly("range_rate_obs_mps", &IslOdtsResults::range_rate_obs_mps)
      .def_readonly("range_resid_m", &IslOdtsResults::range_resid_m)
      .def_readonly("range_rate_resid_mps", &IslOdtsResults::range_rate_resid_mps)
      .def_readonly("cn0_dbhz", &IslOdtsResults::cn0_dbhz);

  // ---- IslOdtsSimulation -------------------------------------------------------

  py::class_<IslOdtsSimulation>(m, "IslOdtsSimulation")
      .def(py::init<IslOdtsConfig>(), py::arg("config"))
      .def("setup", &IslOdtsSimulation::Setup)
      .def("run", &IslOdtsSimulation::Run)
      .def("get_config", &IslOdtsSimulation::GetConfig, py::return_value_policy::reference_internal)
      .def("get_results", &IslOdtsSimulation::GetResults,
           py::return_value_policy::reference_internal);
}

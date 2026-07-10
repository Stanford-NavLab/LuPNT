/**
 * @file py_gnss_odts.cc
 * @brief Python bindings for the lunar GNSS ODTS scenario -- the config/summary structs
 *        and free-function link precompute helpers in
 *        `lupnt/simulations/lunar_gnss_odts/lunar_gnss_odts_simulation.h`, plus the
 *        agent-based `lupnt::LunarGnssOdtsApp`
 *        (`lupnt/applications/lunar_gnss_odts/lunar_gnss_odts_app.h`). A lunar-orbiting
 *        receiver determines its own orbit + clock from cislunar GNSS sidelobe
 *        pseudorange, Doppler, and (optionally) TDCP measurements, run through a UDU EKF
 *        (or UDU stochastic-cloning EKF when TDCP is enabled). Drive it from a
 *        `pnt.Simulation` (a physical `Spacecraft` receiver hosts the app).
 *
 * Config/summary structs mix `std::filesystem::path` and plain scalar members.
 * Path-valued fields are exposed as plain strings, converted to
 * `std::filesystem::path` on the C++ side, following the pattern used for
 * `Sp3Loader`/`AntexLoader`/`GnssConstellation` in `py_gnss.cc`. All other
 * fields bind directly via `def_readwrite`/`def_readonly`.
 */
#include <lupnt/applications/lunar_gnss_odts/lunar_gnss_odts_app.h>
#include <lupnt/lupnt.h>

#include <filesystem>
#include <string>
#include <vector>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

namespace {
  std::vector<std::string> PathsToStrings(const std::vector<std::filesystem::path>& paths) {
    std::vector<std::string> out;
    out.reserve(paths.size());
    for (const auto& p : paths) out.push_back(p.string());
    return out;
  }

  std::vector<std::filesystem::path> StringsToPaths(const std::vector<std::string>& strs) {
    std::vector<std::filesystem::path> out;
    out.reserve(strs.size());
    for (const auto& s : strs) out.push_back(std::filesystem::path(s));
    return out;
  }
}  // namespace

void InitGnssOdts(py::module& m) {
  // ---- ReceiverAppConfig ------------------------------------------------------

  py::class_<ReceiverAppConfig>(m, "ReceiverAppConfig")
      .def(py::init<>())
      .def_readwrite("rate_hz", &ReceiverAppConfig::rate_hz,
                     "Receiver application call rate on the receiver's local clock [Hz]");

  // ---- DesignConfig -------------------------------------------------------------

  py::class_<DesignConfig>(m, "DesignConfig")
      .def(py::init<>())
      .def_property(
          "database_path", [](const DesignConfig& c) { return c.database_path.string(); },
          [](DesignConfig& c, const std::string& s) { c.database_path = std::filesystem::path(s); },
          "Path to a gnss_designs.yaml-style design database (ignored unless set)")
      .def_readwrite("name", &DesignConfig::name, "Design name looked up in the database")
      .def_readwrite("receiver_params", &DesignConfig::receiver_params)
      .def_readwrite("cn0_threshold_dbhz", &DesignConfig::cn0_threshold_dbhz)
      .def_readwrite("cn0_acquisition_threshold_dbhz",
                     &DesignConfig::cn0_acquisition_threshold_dbhz)
      .def_readwrite("cn0_tracking_threshold_dbhz", &DesignConfig::cn0_tracking_threshold_dbhz)
      .def_readwrite("apply_cn0_threshold", &DesignConfig::apply_cn0_threshold)
      .def_readwrite("setup_transmitters", &DesignConfig::setup_transmitters,
                     "Load GPS/Galileo transmit antenna gain patterns so link budget/CN0 "
                     "(including sidelobe reception) is modeled instead of assumed nominal")
      .def_readwrite("receiver_antenna_name", &DesignConfig::receiver_antenna_name,
                     "Receiver antenna gain-pattern name for link-budget/CN0 computation; "
                     "empty uses omni 0 dB")
      .def_readwrite("use_cn0_measurement_sigmas", &DesignConfig::use_cn0_measurement_sigmas)
      .def_readwrite("tx_yaw_dedicated", &DesignConfig::tx_yaw_dedicated,
                     "If True, the transmitter CN0 antenna-gain geometry uses the block-specific "
                     "dedicated eclipse yaw-steering law (GPS/Galileo) rather than the nominal "
                     "Sun-pointing frame; approximation (no maneuver windows), default nominal "
                     "is bit-identical");

  // ---- PlasmaDelayConfig ---------------------------------------------------------

  py::class_<PlasmaDelayConfig>(m, "PlasmaDelayConfig")
      .def(py::init<>())
      .def_readwrite("simulate_truth", &PlasmaDelayConfig::simulate_truth,
                     "Apply the precomputed ionosphere/plasmasphere delay table to truth "
                     "measurements")
      .def_readwrite("model_in_filter", &PlasmaDelayConfig::model_in_filter,
                     "Also apply the delay table to filter (estimated) measurements, instead "
                     "of leaving the delay unmodeled and absorbed by noise inflation")
      .def_readwrite("raytrace_step_size_km", &PlasmaDelayConfig::raytrace_step_size_km)
      .def_readwrite("raytrace_correction", &PlasmaDelayConfig::raytrace_correction)
      .def_readwrite("raytrace_fine_correction", &PlasmaDelayConfig::raytrace_fine_correction)
      .def_readwrite("raytrace_straight_ray", &PlasmaDelayConfig::raytrace_straight_ray)
      .def_readwrite("raytrace_compute_higher_order",
                     &PlasmaDelayConfig::raytrace_compute_higher_order)
      .def_readwrite("raytrace_use_adaptive_step", &PlasmaDelayConfig::raytrace_use_adaptive_step)
      .def_readwrite("raytrace_use_fortran_gcpm", &PlasmaDelayConfig::raytrace_use_fortran_gcpm)
      .def_readwrite("raytrace_cutoff_radius_re", &PlasmaDelayConfig::raytrace_cutoff_radius_re)
      .def_readwrite("raytrace_gradient_step_km", &PlasmaDelayConfig::raytrace_gradient_step_km)
      .def_readwrite("raytrace_correction_tolerance_m",
                     &PlasmaDelayConfig::raytrace_correction_tolerance_m)
      .def_readwrite("raytrace_kp", &PlasmaDelayConfig::raytrace_kp)
      .def_readwrite("raytrace_rz12", &PlasmaDelayConfig::raytrace_rz12)
      .def_readwrite("raytrace_integrator", &PlasmaDelayConfig::raytrace_integrator)
      .def_readwrite("raytrace_correction_method", &PlasmaDelayConfig::raytrace_correction_method)
      .def_readwrite("filter_pseudorange_noise_inflation_m",
                     &PlasmaDelayConfig::filter_pseudorange_noise_inflation_m)
      .def_readwrite("filter_doppler_noise_inflation_hz",
                     &PlasmaDelayConfig::filter_doppler_noise_inflation_hz);

  // ---- ConstellationSourceConfig --------------------------------------------------

  py::class_<ConstellationSourceConfig>(m, "ConstellationSourceConfig")
      .def(py::init<>())
      .def_property(
          "sp3_directory",
          [](const ConstellationSourceConfig& c) { return c.sp3_directory.string(); },
          [](ConstellationSourceConfig& c, const std::string& s) {
            c.sp3_directory = std::filesystem::path(s);
          },
          "Directory scanned for SP3 files when auto_select_sp3 is true")
      .def_readwrite("auto_select_sp3", &ConstellationSourceConfig::auto_select_sp3)
      .def_property(
          "sp3_files",
          [](const ConstellationSourceConfig& c) { return PathsToStrings(c.sp3_files); },
          [](ConstellationSourceConfig& c, const std::vector<std::string>& v) {
            c.sp3_files = StringsToPaths(v);
          },
          "Explicit SP3 file list (ignored when auto_select_sp3 is true)")
      .def_property(
          "antex_file", [](const ConstellationSourceConfig& c) { return c.antex_file.string(); },
          [](ConstellationSourceConfig& c, const std::string& s) {
            c.antex_file = std::filesystem::path(s);
          })
      .def_readwrite("use_all_gps", &ConstellationSourceConfig::use_all_gps)
      .def_readwrite("include_galileo", &ConstellationSourceConfig::include_galileo)
      .def_readwrite("gps_prns", &ConstellationSourceConfig::gps_prns)
      .def_readwrite("galileo_prns", &ConstellationSourceConfig::galileo_prns)
      .def_property(
          "brdc_directory",
          [](const ConstellationSourceConfig& c) { return c.brdc_directory.string(); },
          [](ConstellationSourceConfig& c, const std::string& s) {
            c.brdc_directory = std::filesystem::path(s);
          },
          "Directory scanned for RINEX-nav (BRDC) files when use_broadcast_ephemeris is set")
      .def_property(
          "brdc_files",
          [](const ConstellationSourceConfig& c) { return PathsToStrings(c.brdc_files); },
          [](ConstellationSourceConfig& c, const std::vector<std::string>& v) {
            c.brdc_files = StringsToPaths(v);
          },
          "Explicit RINEX-nav file list (overrides brdc_directory scan)")
      .def_readwrite("use_broadcast_ephemeris", &ConstellationSourceConfig::use_broadcast_ephemeris,
                     "Feed the filter (receiver) model the broadcast transmitter ephemeris "
                     "while truth keeps precise SP3; injects the debiased broadcast error")
      .def_readwrite("debias_broadcast_clock", &ConstellationSourceConfig::debias_broadcast_clock,
                     "Remove the per-constellation median broadcast-minus-precise clock offset")
      .def_readwrite("debias_qzss_radial", &ConstellationSourceConfig::debias_qzss_radial,
                     "Remove the per-QZSS-satellite median radial broadcast-minus-precise orbit "
                     "offset");

  // ---- LunarGnssODTSConfig -------------------------------------------------------

  py::class_<LunarGnssODTSConfig>(m, "LunarGnssODTSConfig")
      .def(py::init<>())
      .def_readwrite("seed", &LunarGnssODTSConfig::seed)
      .def_readwrite("monte_carlo_runs", &LunarGnssODTSConfig::monte_carlo_runs)
      .def_readwrite("duration_s", &LunarGnssODTSConfig::duration_s)
      .def_readwrite("dt_s", &LunarGnssODTSConfig::dt_s)
      .def_readwrite("ephemeris_dt_s", &LunarGnssODTSConfig::ephemeris_dt_s)
      .def_property(
          "output_dir", [](const LunarGnssODTSConfig& c) { return c.output_dir.string(); },
          [](LunarGnssODTSConfig& c, const std::string& s) {
            c.output_dir = std::filesystem::path(s);
          })
      .def_property(
          "links_file", [](const LunarGnssODTSConfig& c) { return c.links_file.string(); },
          [](LunarGnssODTSConfig& c, const std::string& s) {
            c.links_file = std::filesystem::path(s);
          },
          "Output CSV for precomputed link geometry/CN0 (Precompute stage)")
      .def_property(
          "delays_file", [](const LunarGnssODTSConfig& c) { return c.delays_file.string(); },
          [](LunarGnssODTSConfig& c, const std::string& s) {
            c.delays_file = std::filesystem::path(s);
          },
          "Input CSV of precomputed plasma delays (consumed by Run when "
          "plasma.simulate_truth is true)")
      .def_readwrite("receiver_app", &LunarGnssODTSConfig::receiver_app)
      .def_readwrite("design", &LunarGnssODTSConfig::design)
      .def_readwrite("plasma", &LunarGnssODTSConfig::plasma)
      .def_readwrite("constellation", &LunarGnssODTSConfig::constellation)
      .def_readwrite("start_epoch_utc", &LunarGnssODTSConfig::start_epoch_utc)
      .def_readwrite("receiver_a_m", &LunarGnssODTSConfig::receiver_a_m)
      .def_readwrite("receiver_ecc", &LunarGnssODTSConfig::receiver_ecc)
      .def_readwrite("receiver_inc_rad", &LunarGnssODTSConfig::receiver_inc_rad)
      .def_readwrite("receiver_raan_rad", &LunarGnssODTSConfig::receiver_raan_rad)
      .def_readwrite("receiver_argp_rad", &LunarGnssODTSConfig::receiver_argp_rad)
      .def_readwrite("receiver_mean_anomaly_rad", &LunarGnssODTSConfig::receiver_mean_anomaly_rad)
      .def_readwrite("clock_bias_s", &LunarGnssODTSConfig::clock_bias_s)
      .def_readwrite("clock_drift_sps", &LunarGnssODTSConfig::clock_drift_sps)
      .def_readwrite("clock_model", &LunarGnssODTSConfig::clock_model,
                     "Receiver clock model: OCXO, USO, CSAC, MINI_RAFS, RAFS, DSAC")
      .def_readwrite("use_three_state_clock_truth",
                     &LunarGnssODTSConfig::use_three_state_clock_truth,
                     "Generate the truth clock with a 3-state [bias,drift,drift-rate] model "
                     "while the filter stays 2-state")
      .def_readwrite("clock_drift_rate_sps2", &LunarGnssODTSConfig::clock_drift_rate_sps2,
                     "Initial truth clock drift-rate [s/s^2] (3-state truth clock only)")
      .def_readwrite("moon_gravity_degree_truth", &LunarGnssODTSConfig::moon_gravity_degree_truth)
      .def_readwrite("moon_gravity_order_truth", &LunarGnssODTSConfig::moon_gravity_order_truth)
      .def_readwrite("moon_gravity_degree_filter", &LunarGnssODTSConfig::moon_gravity_degree_filter)
      .def_readwrite("moon_gravity_order_filter", &LunarGnssODTSConfig::moon_gravity_order_filter)
      .def_readwrite("moon_gravity_degree_constellation",
                     &LunarGnssODTSConfig::moon_gravity_degree_constellation)
      .def_readwrite("moon_gravity_order_constellation",
                     &LunarGnssODTSConfig::moon_gravity_order_constellation)
      .def_readwrite("include_earth", &LunarGnssODTSConfig::include_earth)
      .def_readwrite("include_sun", &LunarGnssODTSConfig::include_sun)
      .def_readwrite("use_relativity", &LunarGnssODTSConfig::use_relativity,
                     "Apply Moon-centered relativistic clock-rate correction in truth and "
                     "filter propagation (JointOrbitClockDynamics)")
      .def_readwrite("use_srp_truth", &LunarGnssODTSConfig::use_srp_truth)
      .def_readwrite("use_srp_filter", &LunarGnssODTSConfig::use_srp_filter)
      .def_readwrite("srp_coeff_truth_m2_kg", &LunarGnssODTSConfig::srp_coeff_truth_m2_kg)
      .def_readwrite("srp_coeff_filter_m2_kg", &LunarGnssODTSConfig::srp_coeff_filter_m2_kg)
      .def_readwrite("pseudorange_sigma_m", &LunarGnssODTSConfig::pseudorange_sigma_m)
      .def_readwrite("doppler_sigma_hz", &LunarGnssODTSConfig::doppler_sigma_hz)
      .def_readwrite("use_pseudorange", &LunarGnssODTSConfig::use_pseudorange)
      .def_readwrite("use_doppler", &LunarGnssODTSConfig::use_doppler)
      .def_readwrite("use_tdcp", &LunarGnssODTSConfig::use_tdcp,
                     "Add time-differenced carrier phase (TDCP) rows and switch the filter to "
                     "a UDU stochastic-cloning EKF carrying [x_k, x_(k-1)]")
      .def_readwrite("carrier_phase_sigma_m", &LunarGnssODTSConfig::carrier_phase_sigma_m)
      .def_readwrite("tdcp_sigma_m", &LunarGnssODTSConfig::tdcp_sigma_m)
      .def_readwrite("use_ionosphere_free", &LunarGnssODTSConfig::use_ionosphere_free,
                     "Process pseudorange as the dual-frequency (L1+L5) ionosphere-free "
                     "combination; TDCP stays single-frequency (L1)")
      .def_readwrite("filter_pseudorange_noise_inflation_m",
                     &LunarGnssODTSConfig::filter_pseudorange_noise_inflation_m,
                     "Filter-only pseudorange noise inflation [m], added in quadrature on top "
                     "of the C/N0-derived (or IF-combined) sigma")
      .def_readwrite("filter_tdcp_noise_inflation_m",
                     &LunarGnssODTSConfig::filter_tdcp_noise_inflation_m,
                     "Filter-only TDCP noise inflation [m], added in quadrature on top of the "
                     "C/N0-derived carrier-phase sigma")
      .def_readwrite("pseudorange_min_tangent_altitude_m",
                     &LunarGnssODTSConfig::pseudorange_min_tangent_altitude_m,
                     "Reject pseudorange links whose LOS tangent point passes below this Earth "
                     "altitude [m]")
      .def_readwrite("tdcp_min_tangent_altitude_m",
                     &LunarGnssODTSConfig::tdcp_min_tangent_altitude_m,
                     "Reject TDCP links whose LOS tangent point passes below this Earth "
                     "altitude [m]")
      .def_readwrite("estimate_srp_coefficient", &LunarGnssODTSConfig::estimate_srp_coefficient)
      .def_readwrite("initial_srp_coeff_m2_kg", &LunarGnssODTSConfig::initial_srp_coeff_m2_kg)
      .def_readwrite("initial_position_sigma_m", &LunarGnssODTSConfig::initial_position_sigma_m)
      .def_readwrite("initial_velocity_sigma_mps", &LunarGnssODTSConfig::initial_velocity_sigma_mps)
      .def_readwrite("initial_clock_bias_sigma_s", &LunarGnssODTSConfig::initial_clock_bias_sigma_s)
      .def_readwrite("initial_clock_drift_sigma_sps",
                     &LunarGnssODTSConfig::initial_clock_drift_sigma_sps)
      .def_readwrite("initial_srp_coeff_sigma_m2_kg",
                     &LunarGnssODTSConfig::initial_srp_coeff_sigma_m2_kg)
      .def_readwrite("process_accel_sigma_mps2", &LunarGnssODTSConfig::process_accel_sigma_mps2)
      .def_readwrite("process_srp_coeff_sigma_m2_kg_sqrt_s",
                     &LunarGnssODTSConfig::process_srp_coeff_sigma_m2_kg_sqrt_s)
      .def_readwrite("integration_step_s", &LunarGnssODTSConfig::integration_step_s)
      .def_readwrite("precompute_progress_interval_s",
                     &LunarGnssODTSConfig::precompute_progress_interval_s,
                     "Minimum wall-clock seconds between Stage 1 precompute progress prints")
      .def_readwrite("run_progress_interval_epochs",
                     &LunarGnssODTSConfig::run_progress_interval_epochs,
                     "Epochs between live EKF progress prints during Run; 0 = automatic")
      .def_readwrite("debug_print_matrix_epochs", &LunarGnssODTSConfig::debug_print_matrix_epochs,
                     "Print STM and measurement-Jacobian diagnostics for the first N filter epochs")
      .def_readwrite("debug_print_matrix_max_rows",
                     &LunarGnssODTSConfig::debug_print_matrix_max_rows,
                     "Maximum measurement rows included in matrix diagnostics")
      .def_readwrite("precompute_num_threads", &LunarGnssODTSConfig::precompute_num_threads,
                     "Threads for the Stage 1 constellation loop (0 = all cores)");

  // ---- LunarGnssODTSSummary ------------------------------------------------------

  py::class_<LunarGnssODTSSummary>(m, "LunarGnssODTSSummary")
      .def(py::init<>())
      .def_readonly("monte_carlo_index", &LunarGnssODTSSummary::monte_carlo_index)
      .def_readonly("num_epochs", &LunarGnssODTSSummary::num_epochs)
      .def_readonly("final_position_error_m", &LunarGnssODTSSummary::final_position_error_m)
      .def_readonly("final_velocity_error_mps", &LunarGnssODTSSummary::final_velocity_error_mps)
      .def_readonly("final_clock_bias_error_m", &LunarGnssODTSSummary::final_clock_bias_error_m)
      .def_readonly("final_clock_drift_error_mps",
                    &LunarGnssODTSSummary::final_clock_drift_error_mps)
      .def_readonly("final_srp_coeff_error_m2_kg",
                    &LunarGnssODTSSummary::final_srp_coeff_error_m2_kg)
      .def_readonly("rms_position_error_m", &LunarGnssODTSSummary::rms_position_error_m)
      .def_readonly("rms_velocity_error_mps", &LunarGnssODTSSummary::rms_velocity_error_mps);

  // ---- LunarGnssOdtsApp (agent-based) --------------------------------------------
  // The app hosted on the physical `Spacecraft` receiver. Retrieve it from a
  // `pnt.Simulation` via `sim.get_agent("gnss_manager").get_application()` (downcasts here),
  // then read its per-seed `LunarGnssODTSSummary` list; the full time series live in the
  // trajectory_mc<N>.csv / summary.csv outputs under config.output_dir.
  py::class_<LunarGnssOdtsApp, Application, std::shared_ptr<LunarGnssOdtsApp>>(m,
                                                                               "LunarGnssOdtsApp")
      .def("precompute", &LunarGnssOdtsApp::Precompute,
           "Build receiver truth trajectory and GNSS link/CN0 geometry, writing "
           "config.links_file (optional; the EKF Step computes links in-memory otherwise)")
      .def("get_config", &LunarGnssOdtsApp::GetConfig, py::return_value_policy::reference_internal)
      .def("get_summaries", &LunarGnssOdtsApp::GetSummaries,
           py::return_value_policy::reference_internal,
           "Per-seed LunarGnssODTSSummary list; the full time series live in the "
           "trajectory_mc<N>.csv outputs");

  m.def("lunar_gnss_odts_precompute_epoch_count", &LunarGnssODTSPrecomputeEpochCount,
        py::arg("config"), "Return the number of receiver epochs used by GNSS link precompute");
  m.def("lunar_gnss_odts_link_cache_valid", &LunarGnssODTSLinkCacheValid, py::arg("config"),
        "True when config.links_file exists and its metadata matches the config fingerprint");
  m.def("finalize_lunar_gnss_odts_link_cache", &FinalizeLunarGnssODTSLinkCache, py::arg("config"),
        "Write the config fingerprint metadata for an externally assembled links CSV");
  m.def("precompute_lunar_gnss_odts_links_range", &PrecomputeLunarGnssODTSLinksRange,
        py::arg("config"), py::arg("epoch_begin"), py::arg("epoch_end"),
        "Build GNSS link/CN0 geometry for the half-open epoch range [epoch_begin, epoch_end), "
        "writing config.links_file");
}

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

  py::class_<ReceiverAppConfig>(m, "ReceiverAppConfig", "Receiver application scheduling config")
      .def(py::init<>())
      .def_readwrite("rate_hz", &ReceiverAppConfig::rate_hz,
                     "Receiver application call rate on the receiver's local clock [Hz]");

  // ---- DesignConfig -------------------------------------------------------------

  py::class_<DesignConfig>(m, "DesignConfig",
                           "Receiver/antenna design and C/N0-threshold config for the link budget")
      .def(py::init<>())
      .def_property(
          "database_path", [](const DesignConfig& c) { return c.database_path.string(); },
          [](DesignConfig& c, const std::string& s) { c.database_path = std::filesystem::path(s); },
          "Path to a gnss_designs.yaml-style design database (ignored unless set)")
      .def_readwrite("name", &DesignConfig::name, "Design name looked up in the database")
      .def_readwrite("receiver_params", &DesignConfig::receiver_params,
                     "GNSS receiver hardware parameters used in the link budget")
      .def_readwrite("cn0_threshold_dbhz", &DesignConfig::cn0_threshold_dbhz,
                     "Deprecated single C/N0 threshold [dBHz]; used only if the acquisition/"
                     "tracking thresholds are <0")
      .def_readwrite("cn0_acquisition_threshold_dbhz",
                     &DesignConfig::cn0_acquisition_threshold_dbhz,
                     "C/N0 required to acquire a new satellite [dBHz]")
      .def_readwrite("cn0_tracking_threshold_dbhz", &DesignConfig::cn0_tracking_threshold_dbhz,
                     "C/N0 required to keep an existing lock [dBHz]")
      .def_readwrite("apply_cn0_threshold", &DesignConfig::apply_cn0_threshold,
                     "Gate visible satellites by the C/N0 acquisition/tracking thresholds")
      .def_readwrite("setup_transmitters", &DesignConfig::setup_transmitters,
                     "Load GPS/Galileo transmit antenna gain patterns so link budget/CN0 "
                     "(including sidelobe reception) is modeled instead of assumed nominal")
      .def_readwrite("receiver_antenna_name", &DesignConfig::receiver_antenna_name,
                     "Receiver antenna gain-pattern name for link-budget/CN0 computation; "
                     "empty uses omni 0 dB")
      .def_readwrite("use_cn0_measurement_sigmas", &DesignConfig::use_cn0_measurement_sigmas,
                     "Derive per-link measurement sigmas from C/N0 instead of the fixed config "
                     "sigmas")
      .def_readwrite("tx_yaw_dedicated", &DesignConfig::tx_yaw_dedicated,
                     "If True, the transmitter CN0 antenna-gain geometry uses the block-specific "
                     "dedicated eclipse yaw-steering law (GPS/Galileo) rather than the nominal "
                     "Sun-pointing frame; approximation (no maneuver windows), default nominal "
                     "is bit-identical");

  // ---- PlasmaDelayConfig ---------------------------------------------------------

  py::class_<PlasmaDelayConfig>(
      m, "PlasmaDelayConfig",
      "Ionosphere/plasmasphere delay config: truth/filter application flags plus the "
      "ray-trace parameters used to build the delay table")
      .def(py::init<>())
      .def_readwrite("simulate_truth", &PlasmaDelayConfig::simulate_truth,
                     "Apply the precomputed ionosphere/plasmasphere delay table to truth "
                     "measurements")
      .def_readwrite("model_in_filter", &PlasmaDelayConfig::model_in_filter,
                     "Also apply the delay table to filter (estimated) measurements, instead "
                     "of leaving the delay unmodeled and absorbed by noise inflation")
      .def_readwrite("raytrace_step_size_km", &PlasmaDelayConfig::raytrace_step_size_km,
                     "Ray-trace integration step size [km]")
      .def_readwrite("raytrace_correction", &PlasmaDelayConfig::raytrace_correction,
                     "Enable the ray-path bending correction")
      .def_readwrite("raytrace_fine_correction", &PlasmaDelayConfig::raytrace_fine_correction,
                     "Enable the fine ray-path bending correction")
      .def_readwrite("raytrace_straight_ray", &PlasmaDelayConfig::raytrace_straight_ray,
                     "Trace a straight-line ray instead of a bent ray path")
      .def_readwrite("raytrace_compute_higher_order",
                     &PlasmaDelayConfig::raytrace_compute_higher_order,
                     "Include higher-order (beyond first-order) plasma delay terms")
      .def_readwrite("raytrace_use_adaptive_step", &PlasmaDelayConfig::raytrace_use_adaptive_step,
                     "Use adaptive step-size control in the ray-trace integrator")
      .def_readwrite("raytrace_use_fortran_gcpm", &PlasmaDelayConfig::raytrace_use_fortran_gcpm,
                     "Use the Fortran GCPM plasmasphere density model")
      .def_readwrite("raytrace_cutoff_radius_re", &PlasmaDelayConfig::raytrace_cutoff_radius_re,
                     "Outer ray-trace cutoff radius [Earth radii]")
      .def_readwrite("raytrace_gradient_step_km", &PlasmaDelayConfig::raytrace_gradient_step_km,
                     "Finite-difference step for density gradients [km]")
      .def_readwrite("raytrace_correction_tolerance_m",
                     &PlasmaDelayConfig::raytrace_correction_tolerance_m,
                     "Ray-path correction convergence tolerance [m]")
      .def_readwrite("raytrace_kp", &PlasmaDelayConfig::raytrace_kp,
                     "Kp geomagnetic activity index")
      .def_readwrite("raytrace_rz12", &PlasmaDelayConfig::raytrace_rz12,
                     "IRI R12 sunspot index (>0 uses value; -1/-2 historical)")
      .def_readwrite("raytrace_integrator", &PlasmaDelayConfig::raytrace_integrator,
                     "Ray-trace ODE integrator name (e.g. 'RK4')")
      .def_readwrite("raytrace_correction_method", &PlasmaDelayConfig::raytrace_correction_method,
                     "Ray-path correction optimizer name (e.g. 'neldermead')")
      .def_readwrite("filter_pseudorange_noise_inflation_m",
                     &PlasmaDelayConfig::filter_pseudorange_noise_inflation_m,
                     "Filter-only pseudorange noise inflation [m] absorbing unmodeled plasma delay")
      .def_readwrite("filter_doppler_noise_inflation_hz",
                     &PlasmaDelayConfig::filter_doppler_noise_inflation_hz,
                     "Filter-only Doppler noise inflation [Hz] absorbing unmodeled plasma delay");

  // ---- ConstellationSourceConfig --------------------------------------------------

  py::class_<ConstellationSourceConfig>(
      m, "ConstellationSourceConfig",
      "Truth GNSS constellation source: precise SP3+BRDC or numerically propagated almanac, "
      "with satellite selection and broadcast-error options")
      .def(py::init<>())
      .def_property(
          "sp3_directory",
          [](const ConstellationSourceConfig& c) { return c.sp3_directory.string(); },
          [](ConstellationSourceConfig& c, const std::string& s) {
            c.sp3_directory = std::filesystem::path(s);
          },
          "Directory scanned for SP3 files when auto_select_sp3 is true")
      .def_readwrite("auto_select_sp3", &ConstellationSourceConfig::auto_select_sp3,
                     "Auto-select the SP3 files covering the epoch window from sp3_directory")
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
          },
          "ANTEX (.atx) transmitter antenna phase-center/gain file")
      .def_readwrite("use_all_gps", &ConstellationSourceConfig::use_all_gps,
                     "Include all available GPS satellites (ignore gps_prns)")
      .def_readwrite("include_galileo", &ConstellationSourceConfig::include_galileo,
                     "Include Galileo satellites")
      .def_readwrite("gps_prns", &ConstellationSourceConfig::gps_prns,
                     "Explicit GPS PRN list (used when use_all_gps is false)")
      .def_readwrite("galileo_prns", &ConstellationSourceConfig::galileo_prns,
                     "Explicit Galileo PRN list")
      .def_readwrite("include_qzss", &ConstellationSourceConfig::include_qzss,
                     "Include QZSS satellites (L1/L5)")
      .def_readwrite("qzss_prns", &ConstellationSourceConfig::qzss_prns, "Explicit QZSS PRN list")
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
                     "offset")
      .def_readwrite("source", &ConstellationSourceConfig::source,
                     "Truth constellation source: 'sp3_brdc' (default, precise SP3 + BRDC error, "
                     "past epochs only) or 'almanac' (seed Keplerian elements from a YUMA almanac "
                     "or BRDC, numerically propagate with J2 + Sun/Moon -> runs at future epochs)")
      .def_readwrite("almanac_seed", &ConstellationSourceConfig::almanac_seed,
                     "Seed for the future-epoch ('almanac') path: 'sp3' (default) seeds from the "
                     "latest available precise SP3 (real constellation, cm level) then numerically "
                     "propagates; 'yuma'/'brdc' seed from a coarse almanac / broadcast message "
                     "instead. An explicit almanac_file forces the YUMA seed.")
      .def_readwrite("propagation_model", &ConstellationSourceConfig::propagation_model,
                     "Force model for the future-epoch numerical propagation: 'full' (default, "
                     "Earth 8x8 + Sun + Moon) or 'j2' (Earth + J2 zonal only). Both use RK8, which "
                     "preserves the semi-major axis over multi-month spans.")
      .def_property(
          "almanac_file",
          [](const ConstellationSourceConfig& c) { return c.almanac_file.string(); },
          [](ConstellationSourceConfig& c, const std::string& s) {
            c.almanac_file = std::filesystem::path(s);
          },
          "YUMA almanac seed file for source=='almanac' (empty -> seed from the BRDC files)")
      .def_readwrite("synthetic_sise_radial_m", &ConstellationSourceConfig::synthetic_sise_radial_m,
                     "Almanac-mode modeled SISE: per-PRN radial orbit-error std [m]")
      .def_readwrite("synthetic_sise_along_m", &ConstellationSourceConfig::synthetic_sise_along_m,
                     "Almanac-mode modeled SISE: per-PRN along-track orbit-error std [m]")
      .def_readwrite("synthetic_sise_cross_m", &ConstellationSourceConfig::synthetic_sise_cross_m,
                     "Almanac-mode modeled SISE: per-PRN cross-track orbit-error std [m]")
      .def_readwrite("synthetic_sise_clock_m", &ConstellationSourceConfig::synthetic_sise_clock_m,
                     "Almanac-mode modeled SISE: per-PRN clock-error std [m]");

  // ---- LunarGnssODTSConfig -------------------------------------------------------

  py::class_<LunarGnssODTSConfig>(
      m, "LunarGnssODTSConfig",
      "Top-level lunar GNSS ODTS scenario config: Monte-Carlo/timing, receiver orbit+clock, "
      "dynamics fidelity, measurements, filter tuning, and nested design/plasma/constellation "
      "configs")
      .def(py::init<>())
      .def_readwrite("seed", &LunarGnssODTSConfig::seed, "Base RNG seed")
      .def_readwrite("monte_carlo_runs", &LunarGnssODTSConfig::monte_carlo_runs,
                     "Number of Monte-Carlo seeds run")
      .def_readwrite("duration_s", &LunarGnssODTSConfig::duration_s, "Scenario duration [s]")
      .def_readwrite("dt_s", &LunarGnssODTSConfig::dt_s, "Filter measurement/update step [s]")
      .def_readwrite("ephemeris_dt_s", &LunarGnssODTSConfig::ephemeris_dt_s,
                     "Transmitter ephemeris sampling step [s]")
      .def_property(
          "output_dir", [](const LunarGnssODTSConfig& c) { return c.output_dir.string(); },
          [](LunarGnssODTSConfig& c, const std::string& s) {
            c.output_dir = std::filesystem::path(s);
          },
          "Directory for trajectory_mc<N>.csv / summary.csv outputs")
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
      .def_readwrite("receiver_app", &LunarGnssODTSConfig::receiver_app,
                     "Nested receiver application scheduling config")
      .def_readwrite("design", &LunarGnssODTSConfig::design,
                     "Nested receiver/antenna design and C/N0 config")
      .def_readwrite("plasma", &LunarGnssODTSConfig::plasma,
                     "Nested ionosphere/plasmasphere delay config")
      .def_readwrite("constellation", &LunarGnssODTSConfig::constellation,
                     "Nested truth GNSS constellation source config")
      .def_readwrite("start_epoch_utc", &LunarGnssODTSConfig::start_epoch_utc,
                     "Scenario start epoch, ISO 8601 UTC (e.g. '2025-01-01T12:00:00')")
      .def_readwrite("receiver_a_m", &LunarGnssODTSConfig::receiver_a_m,
                     "Receiver orbit semi-major axis [m] (Moon-centered)")
      .def_readwrite("receiver_ecc", &LunarGnssODTSConfig::receiver_ecc,
                     "Receiver orbit eccentricity")
      .def_readwrite("receiver_inc_rad", &LunarGnssODTSConfig::receiver_inc_rad,
                     "Receiver orbit inclination [rad]")
      .def_readwrite("receiver_raan_rad", &LunarGnssODTSConfig::receiver_raan_rad,
                     "Receiver orbit right ascension of ascending node [rad]")
      .def_readwrite("receiver_argp_rad", &LunarGnssODTSConfig::receiver_argp_rad,
                     "Receiver orbit argument of periapsis [rad]")
      .def_readwrite("receiver_mean_anomaly_rad", &LunarGnssODTSConfig::receiver_mean_anomaly_rad,
                     "Receiver orbit initial mean anomaly [rad]")
      .def_readwrite("clock_bias_s", &LunarGnssODTSConfig::clock_bias_s,
                     "Initial receiver clock bias [s]")
      .def_readwrite("clock_drift_sps", &LunarGnssODTSConfig::clock_drift_sps,
                     "Initial receiver clock drift [s/s]")
      .def_readwrite("clock_model", &LunarGnssODTSConfig::clock_model,
                     "Receiver clock model: OCXO, USO, CSAC, MINI_RAFS, RAFS, DSAC")
      .def_readwrite("use_three_state_clock_truth",
                     &LunarGnssODTSConfig::use_three_state_clock_truth,
                     "Generate the truth clock with a 3-state [bias,drift,drift-rate] model "
                     "while the filter stays 2-state")
      .def_readwrite("clock_drift_rate_sps2", &LunarGnssODTSConfig::clock_drift_rate_sps2,
                     "Initial truth clock drift-rate [s/s^2] (3-state truth clock only)")
      .def_readwrite("use_three_state_clock_filter",
                     &LunarGnssODTSConfig::use_three_state_clock_filter,
                     "Estimate a 3-state [bias,drift,drift-rate] clock in the filter (matches the "
                     "manuscript); default is the 2-state [bias,drift] filter")
      .def_readwrite(
          "initial_clock_drift_rate_sigma_sps2",
          &LunarGnssODTSConfig::initial_clock_drift_rate_sigma_sps2,
          "Initial 1-sigma for the filter clock drift-rate [s/s^2] (3-state filter only)")
      .def_readwrite("moon_gravity_degree_truth", &LunarGnssODTSConfig::moon_gravity_degree_truth,
                     "Moon gravity field degree for receiver truth propagation")
      .def_readwrite("moon_gravity_order_truth", &LunarGnssODTSConfig::moon_gravity_order_truth,
                     "Moon gravity field order for receiver truth propagation")
      .def_readwrite("moon_gravity_degree_filter", &LunarGnssODTSConfig::moon_gravity_degree_filter,
                     "Moon gravity field degree for filter propagation")
      .def_readwrite("moon_gravity_order_filter", &LunarGnssODTSConfig::moon_gravity_order_filter,
                     "Moon gravity field order for filter propagation")
      .def_readwrite("moon_gravity_degree_constellation",
                     &LunarGnssODTSConfig::moon_gravity_degree_constellation,
                     "Moon gravity field degree for almanac-mode constellation propagation")
      .def_readwrite("moon_gravity_order_constellation",
                     &LunarGnssODTSConfig::moon_gravity_order_constellation,
                     "Moon gravity field order for almanac-mode constellation propagation")
      .def_readwrite("include_earth", &LunarGnssODTSConfig::include_earth,
                     "Include Earth third-body gravity in propagation")
      .def_readwrite("include_sun", &LunarGnssODTSConfig::include_sun,
                     "Include Sun third-body gravity in propagation")
      .def_readwrite("use_relativity", &LunarGnssODTSConfig::use_relativity,
                     "Apply Moon-centered relativistic clock-rate correction in truth and "
                     "filter propagation (JointOrbitClockDynamics)")
      .def_readwrite("use_srp_truth", &LunarGnssODTSConfig::use_srp_truth,
                     "Include solar radiation pressure in truth propagation")
      .def_readwrite("use_srp_filter", &LunarGnssODTSConfig::use_srp_filter,
                     "Include solar radiation pressure in filter propagation")
      .def_readwrite("srp_coeff_truth_m2_kg", &LunarGnssODTSConfig::srp_coeff_truth_m2_kg,
                     "Truth SRP coefficient (Cr*A/m) [m^2/kg]")
      .def_readwrite("srp_coeff_filter_m2_kg", &LunarGnssODTSConfig::srp_coeff_filter_m2_kg,
                     "Filter SRP coefficient (Cr*A/m) [m^2/kg]")
      .def_readwrite("pseudorange_sigma_m", &LunarGnssODTSConfig::pseudorange_sigma_m,
                     "Pseudorange measurement noise std [m] (when not C/N0-derived)")
      .def_readwrite("doppler_sigma_hz", &LunarGnssODTSConfig::doppler_sigma_hz,
                     "Doppler measurement noise std [Hz] (when not C/N0-derived)")
      .def_readwrite("use_pseudorange", &LunarGnssODTSConfig::use_pseudorange,
                     "Process pseudorange measurements")
      .def_readwrite("use_doppler", &LunarGnssODTSConfig::use_doppler,
                     "Process Doppler measurements")
      .def_readwrite("use_tdcp", &LunarGnssODTSConfig::use_tdcp,
                     "Add time-differenced carrier phase (TDCP) rows and switch the filter to "
                     "a UDU stochastic-cloning EKF carrying [x_k, x_(k-1)]")
      .def_readwrite("carrier_phase_sigma_m", &LunarGnssODTSConfig::carrier_phase_sigma_m,
                     "Carrier-phase measurement noise std [m]")
      .def_readwrite("tdcp_sigma_m", &LunarGnssODTSConfig::tdcp_sigma_m,
                     "Shared (truth+filter) additive TDCP noise std [m]; 0 uses the carrier-phase "
                     "C/N0 floor")
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
      .def_readwrite("estimate_srp_coefficient", &LunarGnssODTSConfig::estimate_srp_coefficient,
                     "Estimate the SRP coefficient as a filter state")
      .def_readwrite("initial_srp_coeff_m2_kg", &LunarGnssODTSConfig::initial_srp_coeff_m2_kg,
                     "Initial filter SRP coefficient estimate (Cr*A/m) [m^2/kg]")
      .def_readwrite("initial_position_sigma_m", &LunarGnssODTSConfig::initial_position_sigma_m,
                     "Initial position uncertainty std [m]")
      .def_readwrite("initial_velocity_sigma_mps", &LunarGnssODTSConfig::initial_velocity_sigma_mps,
                     "Initial velocity uncertainty std [m/s]")
      .def_readwrite("initial_clock_bias_sigma_s", &LunarGnssODTSConfig::initial_clock_bias_sigma_s,
                     "Initial clock-bias uncertainty std [s]")
      .def_readwrite("initial_clock_drift_sigma_sps",
                     &LunarGnssODTSConfig::initial_clock_drift_sigma_sps,
                     "Initial clock-drift uncertainty std [s/s]")
      .def_readwrite("initial_srp_coeff_sigma_m2_kg",
                     &LunarGnssODTSConfig::initial_srp_coeff_sigma_m2_kg,
                     "Initial SRP-coefficient uncertainty std [m^2/kg]")
      .def_readwrite("process_accel_sigma_mps2", &LunarGnssODTSConfig::process_accel_sigma_mps2,
                     "Process-noise acceleration std [m/s^2]")
      .def_readwrite("process_srp_coeff_sigma_m2_kg_sqrt_s",
                     &LunarGnssODTSConfig::process_srp_coeff_sigma_m2_kg_sqrt_s,
                     "SRP-coefficient random-walk process-noise std [m^2/kg/sqrt(s)] (numerical "
                     "floor, not a physical model)")
      .def_readwrite("integration_step_s", &LunarGnssODTSConfig::integration_step_s,
                     "Numerical propagation integration step [s]")
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

  py::class_<LunarGnssODTSSummary>(m, "LunarGnssODTSSummary",
                                   "Per-seed ODTS run summary (final and RMS estimation errors)")
      .def(py::init<>())
      .def_readonly("monte_carlo_index", &LunarGnssODTSSummary::monte_carlo_index,
                    "Monte-Carlo seed index")
      .def_readonly("num_epochs", &LunarGnssODTSSummary::num_epochs, "Number of filter epochs")
      .def_readonly("final_position_error_m", &LunarGnssODTSSummary::final_position_error_m,
                    "Final-epoch position error [m]")
      .def_readonly("final_velocity_error_mps", &LunarGnssODTSSummary::final_velocity_error_mps,
                    "Final-epoch velocity error [m/s]")
      .def_readonly("final_clock_bias_error_m", &LunarGnssODTSSummary::final_clock_bias_error_m,
                    "Final-epoch clock-bias error [m]")
      .def_readonly("final_clock_drift_error_mps",
                    &LunarGnssODTSSummary::final_clock_drift_error_mps,
                    "Final-epoch clock-drift error [m/s]")
      .def_readonly("final_srp_coeff_error_m2_kg",
                    &LunarGnssODTSSummary::final_srp_coeff_error_m2_kg,
                    "Final-epoch SRP-coefficient error [m^2/kg] (NaN if not estimated)")
      .def_readonly("rms_position_error_m", &LunarGnssODTSSummary::rms_position_error_m,
                    "RMS position error over the run [m]")
      .def_readonly("rms_velocity_error_mps", &LunarGnssODTSSummary::rms_velocity_error_mps,
                    "RMS velocity error over the run [m/s]");

  // ---- LunarGnssOdtsApp (agent-based) --------------------------------------------
  // The app hosted on the physical `Spacecraft` receiver. Retrieve it from a
  // `pnt.Simulation` via `sim.get_agent("gnss_manager").get_application()` (downcasts here),
  // then read its per-seed `LunarGnssODTSSummary` list; the full time series live in the
  // trajectory_mc<N>.csv / summary.csv outputs under config.output_dir.
  py::class_<LunarGnssOdtsApp, Application, std::shared_ptr<LunarGnssOdtsApp>>(
      m, "LunarGnssOdtsApp",
      "Coordinator app for the lunar GNSS ODTS scenario, hosted on a Spacecraft receiver; its "
      "scheduled Step runs the whole Monte-Carlo EKF body and writes the CSV outputs")
      .def("precompute", &LunarGnssOdtsApp::Precompute,
           "Build receiver truth trajectory and GNSS link/CN0 geometry, writing "
           "config.links_file (optional; the EKF Step computes links in-memory otherwise)")
      .def("get_config", &LunarGnssOdtsApp::GetConfig, py::return_value_policy::reference_internal,
           "The resolved LunarGnssODTSConfig driving this app")
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

#pragma once

#include <filesystem>
#include <limits>
#include <string>
#include <vector>

#include "lupnt/agents/gnss_constellation.h"
#include "lupnt/simulations/simulation.h"

namespace lupnt {

  struct ReceiverAppConfig {
    double rate_hz = 1.0;
  };

  struct DesignConfig {
    std::filesystem::path database_path = "gnss_designs.yaml";
    std::string name = "lunar_gnss_baseline";
    lupnt::GnssReceiverParams receiver_params;
    double cn0_threshold_dbhz = 15.0;  // deprecated single threshold; used if acq/track are <0
    double cn0_acquisition_threshold_dbhz = -1.0;  // C/N0 to acquire a new satellite [dBHz]
    double cn0_tracking_threshold_dbhz = -1.0;     // C/N0 to keep an existing lock [dBHz]
    bool apply_cn0_threshold = false;
    bool setup_transmitters = false;
    bool use_cn0_measurement_sigmas = false;
    // If true, the transmitter C/N0 antenna-gain geometry uses the block-specific *dedicated*
    // eclipse yaw-steering law (GPS/Galileo) instead of the canonical nominal Sun-pointing
    // frame. Approximation (no maneuver-window state machine); default nominal is bit-identical.
    bool tx_yaw_dedicated = false;
    // Receiver antenna gain pattern for link-budget/C/N0 computation. Empty means omni (0 dB).
    std::string receiver_antenna_name = "moongpsr";
  };

  struct PlasmaDelayConfig {
    bool simulate_truth = false;
    bool model_in_filter = false;
    double raytrace_step_size_km = 100.0;
    bool raytrace_correction = false;
    bool raytrace_fine_correction = false;
    bool raytrace_straight_ray = true;
    bool raytrace_compute_higher_order = true;
    bool raytrace_use_adaptive_step = true;
    bool raytrace_use_fortran_gcpm = true;
    double raytrace_cutoff_radius_re = 4.0;
    double raytrace_gradient_step_km = 1.0;
    double raytrace_correction_tolerance_m = 100.0;
    double raytrace_kp = 3.0;
    double raytrace_rz12 = -1.0;  // IRI R12 sunspot index (>0 uses value; -1/-2 historical)
    std::string raytrace_integrator = "RK4";
    std::string raytrace_correction_method = "neldermead";
    double filter_pseudorange_noise_inflation_m = 0.0;
    double filter_doppler_noise_inflation_hz = 0.0;
  };

  struct ConstellationSourceConfig {
    std::filesystem::path sp3_directory = "../../data/LuPNT_data/ephemeris/gnsslibpy/sp3";
    bool auto_select_sp3 = true;
    std::vector<std::filesystem::path> sp3_files;
    std::filesystem::path antex_file = "../../data/LuPNT_data/gnss/igs20.atx";
    bool use_all_gps = true;
    bool include_galileo = false;
    std::vector<int> gps_prns;
    std::vector<int> galileo_prns;

    // Broadcast (RINEX-nav / BRDC) transmitter ephemeris. When `use_broadcast_ephemeris`
    // is true, the truth measurements keep the precise SP3 transmitter states, but the
    // *filter* (receiver) measurement model is fed the broadcast transmitter position and
    // clock -- i.e. the debiased broadcast-minus-precise error is injected as an unmodeled
    // measurement error, exactly as a real receiver would experience it. Directory is
    // scanned for RINEX-nav files (`*_MN.rnx`) covering the epoch; explicit `brdc_files`
    // override the scan.
    std::filesystem::path brdc_directory;
    std::vector<std::filesystem::path> brdc_files;
    bool use_broadcast_ephemeris = false;
    // Remove the per-constellation systematic clock offset (median over all satellites and
    // epochs of precise-minus-broadcast) from the injected clock error -- see Montenbruck &
    // Steigenberger, "Performance evaluation of the CNAV broadcast ephemeris", J. Nav. 2018.
    bool debias_broadcast_clock = true;
    // Remove the per-satellite median radial (precise-minus-broadcast) orbit offset for QZSS.
    bool debias_qzss_radial = true;
  };

  struct LunarGnssODTSConfig {
    int seed = 42;
    int monte_carlo_runs = 1;
    double duration_s = 3600.0;
    double dt_s = 60.0;
    double ephemeris_dt_s = 60.0;
    std::filesystem::path output_dir = "output/gnss_filtering";
    std::filesystem::path links_file;
    std::filesystem::path delays_file;
    ReceiverAppConfig receiver_app;
    DesignConfig design;
    PlasmaDelayConfig plasma;
    ConstellationSourceConfig constellation;

    std::string start_epoch_utc = "2025-01-01T12:00:00";
    double receiver_a_m = 6541.4e3;
    double receiver_ecc = 0.6;
    double receiver_inc_rad = 65.5 * lupnt::RAD;
    double receiver_raan_rad = 60.0 * lupnt::RAD;
    double receiver_argp_rad = 90.0 * lupnt::RAD;
    double receiver_mean_anomaly_rad = 0.0;
    double clock_bias_s = 0.0;
    double clock_drift_sps = 0.0;
    // Receiver clock model (truth + filter): OCXO, USO, CSAC, MINI_RAFS, RAFS, DSAC.
    std::string clock_model = "OCXO";
    // Generate the *truth* clock with a three-state model [bias, drift, drift-rate] while the
    // filter keeps its two-state [bias, drift] model -- an unmodeled-clock-dynamics scenario.
    bool use_three_state_clock_truth = false;
    // Initial truth clock drift-rate [s/s^2] (only used when use_three_state_clock_truth).
    double clock_drift_rate_sps2 = 0.0;

    int moon_gravity_degree_truth = 20;
    int moon_gravity_order_truth = 20;
    int moon_gravity_degree_filter = 12;
    int moon_gravity_order_filter = 12;
    int moon_gravity_degree_constellation = 12;
    int moon_gravity_order_constellation = 12;
    bool include_earth = true;
    bool include_sun = true;
    bool use_relativity = true;
    bool use_srp_truth = false;
    bool use_srp_filter = false;
    double srp_coeff_truth_m2_kg = 2.0e-3;
    double srp_coeff_filter_m2_kg = 2.0e-3;

    double pseudorange_sigma_m = 2.0;
    double doppler_sigma_hz = 0.1;
    bool use_pseudorange = true;
    bool use_doppler = true;
    bool use_tdcp = false;
    double carrier_phase_sigma_m = 0.01;
    // Shared (truth+filter) additive TDCP sigma; default 0 so truth TDCP noise is the raw
    // carrier-phase C/N0 floor. The filter-only inflation is filter_tdcp_noise_inflation_m.
    double tdcp_sigma_m = 0.0;

    // Dual-frequency ionosphere-free (L1+L5) pseudorange combination (TDCP stays L1).
    bool use_ionosphere_free = false;
    // Filter-only measurement-noise inflation added in quadrature on top of the
    // C/N0-derived sigmas (truth observations use the C/N0 sigmas only).
    double filter_pseudorange_noise_inflation_m = 0.0;
    double filter_tdcp_noise_inflation_m = 0.0;
    // Reject links whose line-of-sight tangent point passes below this Earth altitude
    // [m]; applied per measurement type (pseudorange vs TDCP).
    double pseudorange_min_tangent_altitude_m = 0.0;
    double tdcp_min_tangent_altitude_m = 0.0;

    bool estimate_srp_coefficient = false;
    double initial_srp_coeff_m2_kg = 2.0e-3;
    double initial_position_sigma_m = 100.0;
    double initial_velocity_sigma_mps = 0.1;
    double initial_clock_bias_sigma_s = 1.0e-6;
    double initial_clock_drift_sigma_sps = 1.0e-9;
    double initial_srp_coeff_sigma_m2_kg = 1.0e-3;
    double process_accel_sigma_mps2 = 1.0e-7;
    // Small random-walk floor (not a physical model): keeps the estimated SRP-coefficient
    // variance from shrinking toward zero, which would ill-condition the filter covariance.
    double process_srp_coeff_sigma_m2_kg_sqrt_s = 1.0e-8;
    double integration_step_s = 20.0;
    // Minimum wall-clock seconds between precompute progress prints (Stage 1 link geometry).
    double precompute_progress_interval_s = 3.0;
    // Epochs between live EKF progress prints during Run(). 0 keeps the legacy automatic
    // cadence of about 50 updates over the full run.
    int run_progress_interval_epochs = 0;
    // Print STM and measurement-Jacobian diagnostics for the first N filter epochs. 0 disables.
    int debug_print_matrix_epochs = 0;
    // Maximum measurement rows to include in each diagnostic matrix print.
    int debug_print_matrix_max_rows = 8;
    // Threads for the Stage 1 constellation loop (0 = OpenMP default / all cores; capped to
    // the number of constellation-frequency sets being built).
    int precompute_num_threads = 0;
  };

  struct LunarGnssODTSSummary {
    int monte_carlo_index = 0;
    int num_epochs = 0;
    double final_position_error_m = 0.0;
    double final_velocity_error_mps = 0.0;
    double final_clock_bias_error_m = 0.0;
    double final_clock_drift_error_mps = 0.0;
    double final_srp_coeff_error_m2_kg = std::numeric_limits<double>::quiet_NaN();
    double rms_position_error_m = 0.0;
    double rms_velocity_error_mps = 0.0;
  };

  LunarGnssODTSConfig LoadLunarGnssODTSConfig(const std::filesystem::path& path);
  /// @brief Parse a `LunarGnssODTSConfig` from an already-loaded YAML node (the shared body
  /// of `LoadLunarGnssODTSConfig`). `root` must carry the same section layout as the
  /// scenario file (`simulation:`, `pipeline:`, `truth:`, `constellation:`, `plasma:`,
  /// `dynamics:`, `measurements:`, `filter:`, `receiver_app:`, `design:`); relative paths are
  /// resolved against `base_dir`. Used by the agent-based `LunarGnssOdtsApp` to read its
  /// config from an `application:` block.
  LunarGnssODTSConfig ParseLunarGnssODTSConfig(const Config& root,
                                               const std::filesystem::path& base_dir);
  /// @brief Parse a `plasma:` block (simulate_truth / model_in_filter / raytrace params) into
  /// `cfg`. Shared by the app-block parser and the `world:`-level plasma environment so the
  /// ionosphere/plasmasphere model can be declared once under `world:` (its natural home as a
  /// shared truth property). A null node leaves `cfg` at its defaults.
  void ParsePlasmaDelayConfig(const Config& plasma, PlasmaDelayConfig& cfg);
  /// @brief Resolve a struct-built config for running: auto-select the SP3 products covering
  /// the epoch window when `constellation.auto_select_sp3` is set and no explicit files are
  /// given (matches the struct/Python config path of the former `LunarGnssODTSSimulation`).
  void ResolveLunarGnssODTSConfigForRun(LunarGnssODTSConfig& cfg);
  int LunarGnssODTSPrecomputeEpochCount(const LunarGnssODTSConfig& config);
  bool LunarGnssODTSLinkCacheValid(const LunarGnssODTSConfig& config);
  void FinalizeLunarGnssODTSLinkCache(const LunarGnssODTSConfig& config);
  void PrecomputeLunarGnssODTSLinks(const LunarGnssODTSConfig& config);
  void PrecomputeLunarGnssODTSLinksRange(const LunarGnssODTSConfig& config, int epoch_begin,
                                         int epoch_end);
  /// @brief Run the GNSS ODTS Monte-Carlo body. When `receiver_truth` is non-null it is used
  /// as the receiver's truth trajectory (the physical `Spacecraft` agent's self-propagated
  /// grid, agent-driven path); otherwise the trajectory is built internally from the config
  /// (legacy struct / standalone path). Both must be sampled on the receiver time grid
  /// returned by `LunarGnssODTSReceiverElapsedTimes`.
  std::vector<LunarGnssODTSSummary> RunLunarGnssODTSMonteCarlo(
      const LunarGnssODTSConfig& config, const std::vector<State>* receiver_truth = nullptr);

  /// @brief Elapsed times [s] of the receiver's ODTS epochs (from the receiver-app schedule),
  /// so a host agent can sample its truth on exactly the grid the engine expects.
  VecXd LunarGnssODTSReceiverElapsedTimes(const LunarGnssODTSConfig& config);

}  // namespace lupnt

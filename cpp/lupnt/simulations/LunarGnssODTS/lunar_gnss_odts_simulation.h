#pragma once

#include <filesystem>
#include <limits>
#include <string>
#include <vector>

#include "lupnt/agents/gnss_constellation.h"
#include "lupnt/core/simulation.h"

namespace lupnt {

  struct ReceiverAppConfig {
    double rate_hz = 1.0;
  };

  struct DesignConfig {
    std::filesystem::path database_path = "gnss_designs.yaml";
    std::string name = "lunar_gnss_baseline";
    lupnt::GnssReceiverParams receiver_params;
    double cn0_threshold_dbhz = 15.0;
    bool apply_cn0_threshold = false;
    bool setup_transmitters = false;
    bool use_cn0_measurement_sigmas = false;
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
    double tdcp_sigma_m = 0.02;

    bool estimate_srp_coefficient = false;
    double initial_srp_coeff_m2_kg = 2.0e-3;
    double initial_position_sigma_m = 100.0;
    double initial_velocity_sigma_mps = 0.1;
    double initial_clock_bias_sigma_s = 1.0e-6;
    double initial_clock_drift_sigma_sps = 1.0e-9;
    double initial_srp_coeff_sigma_m2_kg = 1.0e-3;
    double process_accel_sigma_mps2 = 1.0e-7;
    double process_clock_bias_sigma_s_sqrt_s = 1.0e-11;
    double process_clock_drift_sigma_sps_sqrt_s = 1.0e-13;
    double process_srp_coeff_sigma_m2_kg_sqrt_s = 1.0e-10;
    double integration_step_s = 20.0;
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
  void PrecomputeLunarGnssODTSLinks(const LunarGnssODTSConfig& config);
  std::vector<LunarGnssODTSSummary> RunLunarGnssODTSMonteCarlo(const LunarGnssODTSConfig& config);

  class LunarGnssODTSSimulation : public Simulation {
  public:
    explicit LunarGnssODTSSimulation(std::filesystem::path config_path);
    explicit LunarGnssODTSSimulation(LunarGnssODTSConfig config);

    void Setup() override;
    void Precompute() override;
    void Run() override;

    const LunarGnssODTSConfig& GetConfig() const { return config_; }
    const std::vector<LunarGnssODTSSummary>& GetSummaries() const { return summaries_; }

  private:
    std::filesystem::path config_path_;
    LunarGnssODTSConfig config_;
    bool setup_complete_ = false;
    std::vector<LunarGnssODTSSummary> summaries_;
  };

}  // namespace lupnt

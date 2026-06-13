#pragma once

#include <lupnt/lupnt.h>

#include <set>

#include "src/measurement_storage.h"
#include "src/simulation_config.h"

namespace filtering_sim {
  using namespace lupnt;

  using FreqTypePair = std::pair<MeasFreq, MeasType>;

  struct MeasFlag {
    bool use_meas;
    int min_alt;
    int max_alt;
    int process_inv;
  };

  struct FilterConfig {
    // Filter parameters
    bool dmc;  // Dynamic model compensation
    VecXd dmc_beta;
    ProcessNoiseAlgorithm process_noise_algorithm;
    ProcessNoise process_noise;
    MatXd Q_a;
    bool use_adaptive_process_noise = false;

    // Outlier rejection
    bool outlier_rejection;
    double outlier_threshold;

    // SRP
    bool use_srp;
    double srp_initial_error_ratio;
    double true_srp_b_coeff;  // [m^2/kg] True SRP coefficient

    // Initial uncertainties
    double sigma_r0_sat;     // [m] Initial position/bias error
    double sigma_v0_sat;     // [m/s] Initial velocity/drift error
    double sigma_clk0_sat;   // [m] Initial clock offset error
    double sigma_clkd0_sat;  // [m/s] Initial clock drift error

    // Filter dimensions
    int N_t_filter;
    int N_sat_filter;
    int N_clk;

    // Measurements
    std::vector<MeasFreq> ambiguity_meas_freqs = {};
    int n_ambiguity_signal_types = 0;

    bool use_pr;
    bool use_tdcp;
    bool use_graphic;
    bool use_carrier;
    bool use_ionofree_pr;
    bool use_ionofree_carrier;
    bool use_ionofree_tdcp;
    bool use_ionofree_any;

    // Integers
    bool est_integers;
    double tau_ambiguity;
    double q_ambiguity;

    // Clock models
    ClockModel clk_model_sat;
    Ptr<Clock> clock_;
    Ptr<ClockDynamics> clock_dynamics_;
    Ptr<Clock> clock_truth_;
    Ptr<ClockDynamics> clock_dynamics_truth_;
    double inflate_q2 = 1.0;  // Inflation factor for q2 in clock noise

    // Cycle slip ratio
    double cycle_slip_ratio;
    double cycle_slip_magnitude;
    double cycle_slip_detect_sigma;
    double cycle_slip_max_cn0;
    bool fix_cycleslip;

    // Process noise
    double sigma_acc;  // [m/s^2/sqrt(Hz)] Acceleration process noise

    // Noise Models
    double sigma_ure;            // [m] Inflated standard deviation for user range error
    double sigma_iono;           // [m] Inflated standard deviation for ionospheric error
    double sigma_tdcp;           // [m] Inflated standard deviation for TDCP error
    double sigma_graphic;        // [m] Inflated standard deviation for graphic error
    double sigma_tdcp_ionofree;  // [m] Inflated standard deviation for iono-free TDCP error

    double sigma_inflation_code;
    double sigma_inflation_carrier;
    double sigma_inflation_tdcp;
    double sigma_inflation_graphic;
    double sigma_inflation_ionofree_code;
    double sigma_inflation_ionofree_carrier;

    // Measurement Processing
    std::map<FreqTypePair, MeasFlag> meas_flags;

    // Verbose settings
    int verbose_filter_dynamics;
    int verbose_filter_predict;
    int verbose_filter_measurement;
    int verbose_filter_update;
    int verbose_filter_smoothing;

    // Autodiff
    bool use_meas_autodiff;  // Use autodiff for measurement processing

    // UDU decomposition
    bool use_udu;             // Use UDU decomposition for covariance
    bool udu_predict_single;  // Run the measurement model for each individual measurement (vs
                              // correct with linear models)

    // Smoothing
    int N_smooth_iter;       // Number of smoothing iterations
    double smooth_fraction;  // Fraction of time steps to smooth over
    bool run_smoother = false;
    int smooth_start_tidx;  // Time index to start smoothing from

    // Recompute filter results
    bool recompute = false;

    // Dynamics model
    bool joint_orbit_clock_dynamics;  // Use joint orbit-clock dynamics model

    // Measurements to use
    std::set<std::string> measurements;
  };

  struct FilterState {
    MatXd rva_true_sat;
    MatXd rva_est_sat;
    MatXd rva_sigma_sat;
    MatXd clk_true_sat;
    MatXd clk_est_sat;
    MatXd clk_sigma_sat;
    VecXd srp_true_sat;
    VecXd srp_est_sat;
    VecXd srp_sigma_sat;
    MatXd ambiguity_true_sat;      // True ambiguities
    MatXd ambiguity_est_sat;       // Estimated ambiguities
    MatXd ambiguity_sigma_sat;     // Ambiguity sigmas
    MatXd ambiguity_is_estimated;  // Flags for ambiguities being estimated (1: estimated, 0: not
                                   // estimated)
    // Measurement logging
    VecXd num_meas;

    // RTN
    MatXd rtn_err_pos;
    MatXd rtn_err_vel;
    MatXd rtn_err_pos_sigma;
    MatXd rtn_err_vel_sigma;

    MatXd P_est;

    // Time
    double t0_tai;
    VecXd ts_filter;
    VecXd accumulated_relativistic_drift;  // Accumulated relativistic drift over time

    int N_ekf_sat_full;
    int N_ekf_sat_est;  // Total number of states in the EKF
    int N_rva;
    int N_clk;
    int N_srp;
    int N_ambiguity_full;  // Integer ambiguity (L1, L5 for each satellite)
    int N_ambiguity_est;   // Estimated ambiguity (only those being estimated, changes over time)
    std::map<FreqPrnPair, int> freq_prn_to_full_state_index;  // Map (freq, prn) to state index
    std::map<int, FreqPrnPair> full_state_index_to_freq_prn;  // Map state index to (freq, prn)

    int N_t;
  };

  bool AddMeasFlags(const YAML::Node& meas_config, const std::string& meas_type_str,
                    MeasType meas_type, std::map<FreqTypePair, MeasFlag>& meas_flags);
  FilterConfig SetupFilterConfig(const SimulationConfig& sim_config, const VecXd& tspan, int seed,
                                 bool verbose = false);
  FilterState InitializeFilterState(const FilterConfig& cfg, const SimulationConfig& sim_config,
                                    MatXd sat_rv, VecXd tspan, double t0_tai, int seed,
                                    bool verbose = false);

  Ptr<NBodyDynamics> SetupDynamics(bool use_ad, bool use_srp, int sph_order, double srp_coeff,
                                   double dt_prop);

  JointState CreateJointState(const SimulationConfig& sim_config, const FilterConfig& filter_cfg,
                              FilterState& state, bool verbose = false);

}  // namespace filtering_sim

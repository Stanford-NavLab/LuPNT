#include "src/state_meas_manager.h"

#include <H5Lpublic.h>
#include <lupnt/lupnt.h>

#include "lupnt/core/definitions.h"
#include "lupnt/filters/ekf.h"
#include "lupnt/states/state.h"
#include "src/delays.h"
#include "src/dynamics_models.h"
#include "src/gnssmeas_loader.h"
#include "src/udu_filter.h"
#include "src/udu_utils.h"

namespace filtering_sim {
  using namespace lupnt;

  StateMeasManager::StateMeasManager(const SimulationConfig& sim_config,
                                     const GnssMeasLoader& gnss_meas_loader, int seed,
                                     const std::filesystem::path& output_dir, bool verbose) {
    // Create Results File
    double dt = sim_config.dt;
    std::string case_name = sim_config.case_name;
    std::string mc_name = "mc" + std::to_string(seed);

    filtering_results_path_.resize(sim_config.N_smooth_iter);
    smoothing_results_path_.resize(sim_config.N_smooth_iter);

    for (int n_iter = 0; n_iter < sim_config.N_smooth_iter; n_iter++) {
      std::string filter_file_name = "filter_iter_" + std::to_string(n_iter) + ".h5";
      std::string smooth_file_name = "smooth_iter_" + std::to_string(n_iter) + ".h5";
      filtering_results_path_[n_iter] = output_dir / "filter" / sim_config.norbit_dt_str / case_name
                                        / mc_name / filter_file_name;
      smoothing_results_path_[n_iter] = output_dir / "filter" / sim_config.norbit_dt_str / case_name
                                        / mc_name / smooth_file_name;
    }

    bool recompute_filter = sim_config.recompute_filter;
    std::cout << "Recompute Filter: " << (recompute_filter ? "Yes" : "No") << std::endl;

    std::filesystem::create_directories(
        filtering_results_path_[sim_config.N_smooth_iter - 1].parent_path());
    if (std::filesystem::exists(smoothing_results_path_[sim_config.N_smooth_iter - 1])
        && !recompute_filter) {
      std::cout << "Filter and Smoothing results already exist at "
                << smoothing_results_path_[sim_config.N_smooth_iter - 1]
                << ". Skipping filter execution." << std::endl;
      run_filter_ = false;
      return;
    }

    run_filter_ = true;

    meas_storage_ = MeasStorage();
    sim_config_ = sim_config;

    tspan_ = gnss_meas_loader.get_tspan();
    t_tais_ = gnss_meas_loader.get_tais();
    MatX6d rx_posvel_mci = gnss_meas_loader.get_receiver_posvel_mci();
    first_measurement_ = true;
    double t0_tai = gnss_meas_loader.get_t0_tai();

    filter_cfg_ = SetupFilterConfig(sim_config, tspan_, seed, verbose);
    filter_state_ = InitializeFilterState(filter_cfg_, sim_config_, rx_posvel_mci, tspan_, t0_tai,
                                          seed, verbose);
    joint_state_ = CreateJointState(sim_config_, filter_cfg_, filter_state_, verbose);
    SetupGnssMeasurement(gnss_meas_loader, seed);
    InitializeAmbiguityStates(verbose);
    std::cout << "State and Measurement Manager initialized." << std::endl;
  }

  StateMeasManager::~StateMeasManager() {}

  void StateMeasManager::SetupGnssMeasurement(const GnssMeasLoader& gnssmeas_loader, int seed) {
    // Random variables
    RandomEngine::SetSeed(seed);
    auto rng = RandomEngine::Get();
    std::uniform_real_distribution<double> real_dist(0.0, 1.0);
    std::uniform_int_distribution<int> int_dist(-filter_cfg_.cycle_slip_magnitude,
                                                filter_cfg_.cycle_slip_magnitude);

    VecXd tspan = gnssmeas_loader.get_tspan();
    int lent = filter_state_.N_t;
    auto meas_config = sim_config_.config["measurements"];

    std::vector<GnssMeas> all_meas = gnssmeas_loader.get_data();

    auto pbar = Logger::GetProgressBar(lent, "Generating GNSS Measurements", "Main");
    std::atomic<std::size_t> counter{0};
    int total_meas = 0;

    for (int k = 0; k < lent; k++) {
      std::vector<int> meas_types;
      std::vector<double> y_meas_vec;
      std::vector<double> Rdiag_meas_vec;

      VecXd clk_true_sat = filter_state_.clk_true_sat.row(k);

      // Collect measurements at time index k
      std::vector<int> data_db_indices = gnssmeas_loader.get_dataidx_for_tidx(k);
      int n_data_tidx = data_db_indices.size();

      // Judget the previous, current, and next time index tuples
      std::vector<ConstFreqPrnTuple> prev_tuples = {};
      if (k > 0) {
        prev_tuples = gnssmeas_loader.get_tuples_for_tidx(k - 1);
      }
      std::vector<ConstFreqPrnTuple> next_tuples = {};
      if (k < lent - 1) {
        next_tuples = gnssmeas_loader.get_tuples_for_tidx(k + 1);
      }

      // insert empty measurement if no data
      if (data_db_indices.empty()) {
        continue;
      }

      bool acquired_this_step = false;
      bool final_track_this_step = false;

      // Process each measurement
      std::map<FreqPrnPair, int>
          accumulated_slips;  // Map to track accumulated slips per unique PRN

      for (int i = 0; i < n_data_tidx; i++) {
        int data_idx = data_db_indices[i];
        const GnssMeas& meas = all_meas[data_idx];
        std::string signal = meas.signal;
        GnssConst gnss_const = StringToGnssConst(meas.gnss_const);
        MeasFreq meas_freq = StringToMeasFreq(signal);
        int unique_prn = GetUniquePrn(meas.prn, gnss_const);

        FreqPrnPair freq_prn_key = std::make_pair(meas_freq, unique_prn);

        // Judge if this signal is acquired this step
        ConstFreqPrnTuple curr_tuple = std::make_tuple(meas.gnss_const, meas.signal, meas.prn);
        if (std::find(prev_tuples.begin(), prev_tuples.end(), curr_tuple) == prev_tuples.end()) {
          // Not found in previous step, so newly acquired
          acquired_this_step = true;
        }
        if (std::find(next_tuples.begin(), next_tuples.end(), curr_tuple) == next_tuples.end()) {
          // Not found in next step, so final track
          final_track_this_step = true;
        }

        // Compute range in MCI frame
        Vec3d rx_pos_mci = meas.posvel_rx_mci.head<3>();
        Vec3d tx_pos_mci = meas.posvel_tx_mci.head<3>();
        double range_mci = (rx_pos_mci - tx_pos_mci).norm();
        Vec6d gnss_rv_ephem_mci
            = meas.posvel_tx_ephem_mci;  // Use ephemeris-based position-velocity
        double gnss_clock_ephem = meas.clockbias_ephem;

        double lamda = get_wavelength(meas_freq);

        // 1. Reset accumulation on new acquisition
        if (acquired_this_step) {
          accumulated_slips[freq_prn_key] = 0;
        }

        // Ensure the entry exists
        if (accumulated_slips.find(freq_prn_key) == accumulated_slips.end()) {
          accumulated_slips[freq_prn_key] = 0;
        }

        double shapiro_delay_m = meas.shapiro_delay;
        double relativistic_delay_m = 0.0;
        if (!filter_cfg_.joint_orbit_clock_dynamics) {
          // Add relativistic delay only if not using joint dynamics
          relativistic_delay_m = meas.relativistic_delay;
        }

        // Always generate pseudorange and carrier
        double pr = range_mci + (clk_true_sat(0) - meas.clockbias_tx) + meas.total_delay_m
                    + shapiro_delay_m + relativistic_delay_m;  // deteministic part only
        double pr_noise = SampleNormal(0.0, meas.sigma_range);
        double pr_y = pr + pr_noise;
        double carrier_range_true = range_mci + (clk_true_sat(0) - meas.clockbias_tx)
                                    + meas.total_carrier_delay_m + shapiro_delay_m
                                    + relativistic_delay_m;  // deteministic part only
        double carrier_noise = SampleNormal(0.0, meas.sigma_carrier);
        double carrier_range_y = carrier_range_true + carrier_noise;
        // 2. Determine NEW cycle slip for this epoch
        double cn0 = meas.cn0;
        int new_cycle_slip = 0;

        if (filter_cfg_.cycle_slip_ratio > 0.0 && cn0 <= filter_cfg_.cycle_slip_max_cn0) {
          if (acquired_this_step) {
            // No cycle slip on first acquisition
            new_cycle_slip = 0;
          } else {
            // Check if true cycle slip occurs
            double slip_rand = real_dist(rng);
            if (slip_rand < filter_cfg_.cycle_slip_ratio) {
              // Cycle slip occurs
              int slip_magnitude = int_dist(rng);
              if (slip_magnitude == 0) {
                slip_magnitude = 1;  // Ensure non-zero slip
              }
              new_cycle_slip = slip_magnitude;
              accumulated_slips[freq_prn_key] += new_cycle_slip;
            }
          }
        }

        // 3, Apply accumulated slips to carrier phase
        carrier_range_y += accumulated_slips[freq_prn_key] * lamda;

        int n_cycles = 0;
        double carrier_y_phase = 0.0;

        // Unmodeled delays
        double unmodeled_pr_delay = meas.total_delay_m + meas.ephem_range_error_m;
        double unmodeled_carrier_delay = meas.total_carrier_delay_m + meas.ephem_range_error_m;

        // time
        double tspan_meas = tspan_(k);
        double tai_meas = t_tais_(k);

        if (acquired_this_step) {
          // On first acquisition, reset carrier phase integer ambiguity
          n_cycles = static_cast<int>(std::floor(carrier_range_y / lamda));
          carrier_y_phase = carrier_range_y - n_cycles * lamda;
          meas_storage_.add_measurement(
              k, tspan_meas, tai_meas, meas_freq, unique_prn, pr_y, carrier_range_y,
              carrier_y_phase, n_cycles, accumulated_slips[freq_prn_key], unmodeled_pr_delay,
              unmodeled_carrier_delay, meas.sigma_range, meas.sigma_carrier, gnss_rv_ephem_mci,
              gnss_clock_ephem, acquired_this_step, final_track_this_step);
        } else {
          // Subsequent measurements
          int prev_index = meas_storage_.get_index(k - 1, meas_freq, unique_prn);
          assert(prev_index != -1);

          // Note that we use the OLD n_cycles (from prev_index).
          // Since carrier_range_y has jumped, but n_cycles has not.
          // the residuals (Measurement - Model) will correctly show the slip.
          n_cycles = meas_storage_.carrier_phase_integer[prev_index];
          carrier_y_phase = carrier_range_y - n_cycles * lamda;
          meas_storage_.add_measurement(
              k, tspan_meas, tai_meas, meas_freq, unique_prn, pr_y, carrier_range_y,
              carrier_y_phase, n_cycles, accumulated_slips[freq_prn_key], unmodeled_pr_delay,
              unmodeled_carrier_delay, meas.sigma_range, meas.sigma_carrier, gnss_rv_ephem_mci,
              gnss_clock_ephem, acquired_this_step, final_track_this_step);
        }
        total_meas++;
      }  // For each data index

      auto c = counter.fetch_add(1);
      if (c % 100) pbar->Update(c);
    }  // For each time step

    pbar->Finish();  // End of time steps

    // Initialize Ambiguity States
    std::cout << "Total Time Steps with GNSS Measurements: " << lent << std::endl;
    std::cout << "Total GNSS Measurements Generated: " << total_meas << std::endl;
    std::cout << "Measurement Setup Finished" << std::endl;

    return;
  }  // End of SetupGnssMeasurement

  /**
   * Initialize ambiguity states in the filter state
   * @param filter_config Filter configuration
   * @param filter_state Current filter state
   * @param meas_manager Measurement manager with measurement data
   * @param verbose Whether to print initialization info
   * @return Updated filter state with initialized ambiguity states
   */
  void StateMeasManager::InitializeAmbiguityStates(bool verbose) {
    FilterState state = filter_state_;
    int n_signal_types = filter_cfg_.n_ambiguity_signal_types;  // L1, L5, L1-L5

    // Determine number of ambiguities from measurements
    int N_prns = GetUniquePrnCount();

    int N_ambiguities = 0;

    if (filter_cfg_.est_integers) {
      N_ambiguities = n_signal_types * N_prns;
    }
    state.N_ambiguity_full = N_ambiguities;
    state.N_ekf_sat_full = state.N_rva + state.N_clk + state.N_srp + state.N_ambiguity_full;
    state.N_ambiguity_est = 0;  // For the first step, no ambiguities are estimated
    state.N_ekf_sat_est = state.N_rva + state.N_clk + state.N_srp + state.N_ambiguity_est;

    int integer_state_start_idx = state.N_rva + state.N_clk + state.N_srp;

    // Map the unique PRNs to state indice
    if (filter_cfg_.est_integers) {
      state.freq_prn_to_full_state_index.clear();
      state.full_state_index_to_freq_prn.clear();
      std::vector<int> unique_prns = GetAllUniquePrns();
      for (int i = 0; i < N_prns; i++) {
        int prn = unique_prns[i];
        for (int freq = 0; freq < n_signal_types; freq++) {
          MeasFreq meas_freq = filter_cfg_.ambiguity_meas_freqs[freq];
          state.freq_prn_to_full_state_index[{meas_freq, prn}]
              = integer_state_start_idx + n_signal_types * i + freq;
          state.full_state_index_to_freq_prn[integer_state_start_idx + n_signal_types * i + freq]
              = {meas_freq, prn};
        }
      }

      // Resize ambiguity state arrays
      state.ambiguity_true_sat.resize(state.N_t, N_ambiguities);
      state.ambiguity_est_sat.resize(state.N_t, N_ambiguities);
      state.ambiguity_sigma_sat.resize(state.N_t, N_ambiguities);
      state.ambiguity_is_estimated.resize(state.N_t, N_ambiguities);

      // Initialize true ambiguities
      for (int i = 0; i < filter_cfg_.N_t_filter; i++) {
        std::vector<int> prns = GetTrackedPrnsAtTimeIndex(i);
        for (size_t j = 0; j < prns.size(); j++) {
          for (int freq = 0; freq < n_signal_types; freq++) {
            MeasFreq meas_freq = filter_cfg_.ambiguity_meas_freqs[freq];
            int prn = prns[j];
            int state_index = state.freq_prn_to_full_state_index[{meas_freq, prn}];
            state.ambiguity_true_sat(i, state_index - integer_state_start_idx)
                = GetTrueAmbiguity(i, meas_freq, prn);
          }
        }
      }

      if (verbose) {
        std::cout << "Initialized " << N_ambiguities << " ambiguity states for " << N_prns
                  << " unique PRNs." << std::endl;
      }
    } else {
      if (verbose) {
        std::cout << "No ambiguity states to initialize." << std::endl;
      }
    }

    // Update the filter state
    filter_state_ = state;
  }

  std::vector<int> StateMeasManager::GetAllUniquePrns() const {
    std::vector<int> prns(meas_storage_.unique_prns.begin(), meas_storage_.unique_prns.end());
    return prns;
  }

  std::vector<int> StateMeasManager::GetTrackedPrnsAtTimeIndex(int tidx) const {
    std::vector<int> prns;
    auto it = meas_storage_.freqs_prns_at_tidx.find(tidx);
    if (it != meas_storage_.freqs_prns_at_tidx.end()) {
      for (const auto& pair : it->second) {
        prns.push_back(pair.second);
      }
    }
    return prns;
  }

  int StateMeasManager::GetTrueAmbiguity(int tidx, MeasFreq freq, int prn) const {
    // Placeholder: In actual implementation, retrieve from true state or measurement storage
    int index = meas_storage_.get_index(tidx, freq, prn);
    if (index != -1) {
      return meas_storage_.carrier_phase_integer[index];
    } else {
      return 0;  // Default if not found (or multi-freq combination)
    }
  }

  StateCovPair StateMeasManager::FullStateToEstState(const VecXd& x_full, const MatXd& P_full,
                                                     std::vector<int>& ambig_full_state_indices,
                                                     bool is_tdcp) {
    // Placeholder for converting full state to estimated state
    // Implement the logic based on tracked PRNs and filter configuration
    int N_sat_state = filter_state_.N_rva + filter_state_.N_clk + filter_state_.N_srp;
    int N_est_state = N_sat_state + filter_state_.N_ambiguity_est;
    int N_full_state = filter_state_.N_ekf_sat_full;
    VecXd x_est = VecXd::Zero(N_est_state);
    MatXd P_est = MatXd::Zero(N_est_state, N_est_state);
    if (is_tdcp) {
      x_est.resize(2 * N_est_state);
      P_est.resize(2 * N_est_state, 2 * N_est_state);
    }

    // First copy satellite states
    for (int i = 0; i < N_sat_state; i++) {
      x_est(i) = x_full(i);
      for (int j = 0; j < N_sat_state; j++) {
        P_est(i, j) = P_full(i, j);
      }

      if (is_tdcp) {
        x_est(i + N_est_state) = x_full(i + N_full_state);
        for (int j = 0; j < N_sat_state; j++) {
          P_est(i + N_est_state, j + N_est_state) = P_full(i + N_full_state, j + N_full_state);
          P_est(i + N_est_state, j) = P_full(i + N_full_state, j);
          P_est(i, j + N_est_state) = P_full(i, j + N_full_state);
        }
      }
    }

    // Then copy ambiguity states
    if (ambig_full_state_indices.empty()) {
      return {x_est, P_est};
    }

    int idx = N_sat_state;
    for (const auto& ambig_idx : ambig_full_state_indices) {
      // Skip if no ambiguity for this index, or multiple indices (dual frequency)
      x_est(idx) = x_full(ambig_idx);

      // Copy ambig-sat state covariances
      P_est(idx, idx) = P_full(ambig_idx, ambig_idx);
      for (int i = 0; i < N_sat_state; i++) {
        P_est(idx, i) = P_full(ambig_idx, i);
        P_est(i, idx) = P_full(i, ambig_idx);
      }

      if (is_tdcp) {
        x_est(idx + N_est_state) = x_full(ambig_idx + N_full_state);
        P_est(idx + N_est_state, idx + N_est_state)
            = P_full(ambig_idx + N_full_state, ambig_idx + N_full_state);
        for (int i = 0; i < N_sat_state; i++) {
          // bottom left corner (add N_est_state/N_full_state to row)
          P_est(idx + N_est_state, i) = P_full(ambig_idx + N_full_state, i);
          P_est(i + N_est_state, idx) = P_full(i + N_full_state, ambig_idx);
          // top right corner (add N_est_state/N_full_state to column)
          P_est(idx, i + N_est_state) = P_full(ambig_idx, i + N_full_state);
          P_est(i, idx + N_est_state) = P_full(i, ambig_idx + N_full_state);
          // bottom right corner (add N_est_state/N_full_state to both row and column)
          P_est(idx + N_est_state, i + N_est_state)
              = P_full(ambig_idx + N_full_state, i + N_full_state);
          P_est(i + N_est_state, idx + N_est_state)
              = P_full(i + N_full_state, ambig_idx + N_full_state);
        }
      }

      idx++;
    }

    return {x_est, P_est};
  }

  StateCovPair StateMeasManager::EstStateToFullState(const VecXd& x_est, const MatXd& P_est,
                                                     std::vector<int>& ambig_full_state_indices,
                                                     bool is_tdcp) {
    // Placeholder for converting estimated state to full state
    // Implement the logic based on tracked PRNs and filter configuration

    int N_sat_state = filter_state_.N_rva + filter_state_.N_clk + filter_state_.N_srp;
    int N_est_state = N_sat_state + filter_state_.N_ambiguity_est;
    int N_full_state = filter_state_.N_ekf_sat_full;

    VecXd x_full = VecXd::Zero(N_full_state);
    MatXd P_full = MatXd::Zero(N_full_state, N_full_state);

    if (is_tdcp) {
      x_full.resize(2 * N_full_state);
      P_full.resize(2 * N_full_state, 2 * N_full_state);
    }

    // Copy satellite states
    for (int i = 0; i < N_sat_state; i++) {
      x_full(i) = x_est(i);
      for (int j = 0; j < N_sat_state; j++) {
        P_full(i, j) = P_est(i, j);
      }

      if (is_tdcp) {
        x_full(i + N_full_state) = x_est(i + N_est_state);
        for (int j = 0; j < N_sat_state; j++) {
          P_full(i + N_full_state, j + N_full_state) = P_est(i + N_est_state, j + N_est_state);
          P_full(i + N_full_state, j) = P_est(i + N_est_state, j);
          P_full(i, j + N_full_state) = P_est(i, j + N_est_state);
        }
      }
    }

    if (ambig_full_state_indices.empty()) {
      return {x_full, P_full};
    }

    // Copy ambiguity states
    int idx = N_sat_state;
    for (const auto& ambig_idx : ambig_full_state_indices) {
      // Skip if no ambiguity for this index, or multiple indices (dual frequency)
      x_full(ambig_idx) = x_est(idx);
      // Copy ambig-sat state covariances
      P_full(ambig_idx, ambig_idx) = P_est(idx, idx);
      for (int i = 0; i < N_sat_state; i++) {
        P_full(ambig_idx, i) = P_est(idx, i);
        P_full(i, ambig_idx) = P_est(i, idx);
      }

      if (is_tdcp) {
        x_full(ambig_idx + N_full_state) = x_est(idx + N_est_state);
        P_full(ambig_idx + N_full_state, ambig_idx + N_full_state)
            = P_est(idx + N_est_state, idx + N_est_state);
        for (int i = 0; i < N_sat_state; i++) {
          // bottom left corner (add N_full_state/N_est_state to row)
          P_full(ambig_idx + N_full_state, i) = P_est(idx + N_est_state, i);
          P_full(i + N_full_state, ambig_idx) = P_est(i + N_est_state, idx);
          // top right corner (add N_full_state/N_est_state to column)
          P_full(ambig_idx, i + N_full_state) = P_est(idx, i + N_est_state);
          P_full(i, ambig_idx + N_full_state) = P_est(i, idx + N_est_state);
          // bottom right corner (add N_full_state/N_est_state to both row and column)
          P_full(ambig_idx + N_full_state, i + N_full_state)
              = P_est(idx + N_est_state, i + N_est_state);
          P_full(i + N_full_state, ambig_idx + N_full_state)
              = P_est(i + N_est_state, idx + N_est_state);
        }
      }

      idx++;
    }

    // For simplicity, return the input as is
    return {x_full, P_full};
  }

  int StateMeasManager::FullIndexToEstIndex(
      int full_index, const std::vector<int>& ambig_full_state_indices) const {
    int N_sat_state = filter_state_.N_rva + filter_state_.N_clk + filter_state_.N_srp;
    if (full_index < N_sat_state) {
      return full_index;  // Satellite state
    } else {
      // Ambiguity state
      auto it
          = std::find(ambig_full_state_indices.begin(), ambig_full_state_indices.end(), full_index);
      if (it != ambig_full_state_indices.end()) {
        return N_sat_state + std::distance(ambig_full_state_indices.begin(), it);
      } else {
        throw std::runtime_error("Full index not found in ambiguity state indices");
      }
    }
  }

  int StateMeasManager::EstIndexToFullIndex(
      int est_index, const std::vector<int>& ambig_full_state_indices) const {
    int N_sat_state = filter_state_.N_rva + filter_state_.N_clk + filter_state_.N_srp;
    if (est_index < N_sat_state) {
      return est_index;  // Satellite state
    } else {
      // Ambiguity state
      int ambig_idx = est_index - N_sat_state;
      if (ambig_idx >= 0 && ambig_idx < ambig_full_state_indices.size()) {
        return ambig_full_state_indices[ambig_idx];
      } else {
        throw std::runtime_error("Estimated index out of range for ambiguity states");
      }
    }
  }

  double StateMeasManager::GenerateNewAmbiguityEstimate(int meas_index) {
    double pr_meas = meas_storage_.pseudorange[meas_index];
    double carrier_meas = meas_storage_.carrier_phase[meas_index];
    MeasFreq meas_freq = meas_storage_.freq[meas_index];

    // Simple rounding method for new ambiguity estimate
    double lamda = get_wavelength(meas_freq);  // Use provided measurement frequency
    double N_float = (pr_meas - carrier_meas) / lamda - AMBIGUITY_OFFSET_N;
    return N_float;
  }

  double StateMeasManager::GenerateNewAmbiguityCovariance(int meas_index) {
    double sigma_pr = meas_storage_.sigma_pseudorange[meas_index];
    double sigma_carrier = meas_storage_.sigma_carrier_phase[meas_index];
    MeasFreq meas_freq = meas_storage_.freq[meas_index];
    double lamda = get_wavelength(meas_freq);  // Use provided measurement frequency
    double sigma_N = std::sqrt(sigma_pr * sigma_pr + sigma_carrier * sigma_carrier) / lamda;

    // PR - Carrier difference can be off by multiple wavelengths due to ionospheric delay
    // The ionospheric delay can be up to ~100 meters, and is doubled in the PR - Carrier difference
    double max_delay_m = 100.0;  // Assume max iono delay of 100 meters for initial covariance
    sigma_N = (2 * max_delay_m) / lamda + sigma_N;  // Large initial uncertainty
    double P_N = sigma_N * sigma_N;

    return P_N;
  }

  /******************************************************
   * Measurement Manipulation at Each Time Step
   *******************************************************/

  MeasData StateMeasManager::CreateEmptyMeasData(int tidx) {
    MeasData meas_data;
    meas_data.tidx = tidx;
    meas_data.tspan = meas_storage_.tspan_at_tidx[tidx];
    meas_data.t_tai = meas_storage_.tai_at_tidx[tidx];
    meas_data.n_meas = 0;
    meas_data.z = VecXd::Zero(0);
    meas_data.R_diag = VecXd::Zero(0);
    meas_data.N_integer_est = VecXd::Zero(0);
    meas_data.ambig_range_est = VecXd::Zero(0);
    meas_data.true_delays = VecXd::Zero(0);
    meas_data.ambig_est_state_index_for_allmeas = {};
    meas_data.ambig_full_state_indices = {};
    meas_data.prns = {};
    meas_data.freqs = {};
    meas_data.types = {};
    meas_data.freq_prn_pairs = {};
    meas_data.r_ephem_mci = {};
    meas_data.clock_ephem = {};
    meas_data.r_prev_ephem_mci = {};
    meas_data.clock_prev_ephem = {};

    return meas_data;
  }

  bool StateMeasManager::IsUseMeasurement(MeasFreq freq, MeasType type, double tangential_alt,
                                          int tidx) const {
    auto it = filter_cfg_.meas_flags.find({freq, type});
    if (it != filter_cfg_.meas_flags.end()) {
      const MeasFlag& flag = it->second;
      if (!flag.use_meas) {
        return false;
      }
      if (tangential_alt < flag.min_alt || tangential_alt > flag.max_alt) {
        return false;
      }
      if ((tidx % flag.process_inv) != 0) {
        return false;
      }
      return true;
    }
    // If no specific flag found, do not use the measurement
    return false;
  }

  void StateMeasManager::AddSingleFreqMeasurement(int tidx, State& x_est_full, VecMeasData& vmd) {
    // Information
    bool use_pr = filter_cfg_.use_pr;
    bool use_carrier = filter_cfg_.use_carrier;
    bool use_graphic = filter_cfg_.use_graphic;
    bool use_tdcp = filter_cfg_.use_tdcp;

    // State Information ----------------------------------------------
    Vec3 rx_pos_mci = x_est_full.segment<3>(0);  // Receiver position from state estimate
    double t_tai = t_tais_(tidx);
    double min_alt = 0;

    // Compute Minimum Satellite Altitude
    double sigma_ure = filter_cfg_.sigma_ure;
    double sigma_iono = filter_cfg_.sigma_iono;
    double sigma_tdcp = filter_cfg_.sigma_tdcp;
    double sigma_graphic = filter_cfg_.sigma_graphic;

    int n_signals = static_cast<int>(vmd.db_indices.size());
    vmd.use_signal_carrier_vec.resize(n_signals, false);

    for (int i = 0; i < n_signals; i++) {
      int dbidx = vmd.db_indices[i];
      double pr_meas = meas_storage_.pseudorange[dbidx];
      MeasFreq freq = meas_storage_.freq[dbidx];
      int prn = meas_storage_.prn[dbidx];
      int cycleslip = meas_storage_.cycle_slip[dbidx];
      double carrier_phase_meas = meas_storage_.carrier_phase[dbidx];

      vmd.idx_to_carrier_phase_map[dbidx] = carrier_phase_meas;

      // Get Error information
      double sigma_pr_noise = meas_storage_.sigma_pseudorange[dbidx];
      double sigma_carrier_noise = meas_storage_.sigma_carrier_phase[dbidx];
      double var_pr = pow(sigma_pr_noise, 2) + pow(sigma_ure, 2) + pow(sigma_iono, 2);
      double var_carrier = pow(sigma_carrier_noise, 2) + pow(sigma_ure, 2) + pow(sigma_iono, 2);
      double var_graphic
          = pow(sigma_pr_noise, 2) / 4.0 + pow(sigma_carrier_noise, 2) / 4.0 + pow(sigma_ure, 2);
      Vec3d r_ephem_mci = meas_storage_.gnss_rv_ephem_mci[dbidx].head<3>();
      double clk_ephem = meas_storage_.gnss_clock_ephem[dbidx];
      double pr_delay = meas_storage_.pr_delay[dbidx];
      double carrier_delay = meas_storage_.carrier_delay[dbidx];

      // Compute Minimum Satellite Altitude
      min_alt = ComputeMinEarthAltitude(t_tai, rx_pos_mci, r_ephem_mci);

      if (IsUseMeasurement(freq, MeasType::CODE, min_alt, tidx)) {
        double R_code = var_pr;
        vmd.z_meas_vec.push_back(pr_meas);
        vmd.R_diag_vec.push_back(R_code);
        vmd.meas_freqs_vec.push_back(freq);
        vmd.meas_types_vec.push_back(MeasType::CODE);
        vmd.prns_vec.push_back(prn);
        vmd.r_ephem_mci_vec.push_back(r_ephem_mci);
        vmd.clk_ephem_vec.push_back(clk_ephem);
        vmd.true_delays_vec.push_back(pr_delay);  // Placeholder for true delays
        vmd.min_alt_vec.push_back(min_alt);
        vmd.meas_to_signal_index_vec.push_back(i);
        vmd.true_cycle_slip_vec.push_back(0);  // No cycle slip for code measurements
      }
      if (IsUseMeasurement(freq, MeasType::CARRIER, min_alt, tidx)) {
        double R_carrier = var_carrier;
        vmd.z_meas_vec.push_back(carrier_phase_meas);
        vmd.R_diag_vec.push_back(R_carrier);
        vmd.meas_freqs_vec.push_back(freq);
        vmd.meas_types_vec.push_back(MeasType::CARRIER);
        vmd.prns_vec.push_back(prn);
        vmd.r_ephem_mci_vec.push_back(r_ephem_mci);
        vmd.clk_ephem_vec.push_back(clk_ephem);
        vmd.true_delays_vec.push_back(carrier_delay);  // Placeholder for true delays
        vmd.min_alt_vec.push_back(min_alt);
        vmd.meas_to_signal_index_vec.push_back(i);
        vmd.use_signal_carrier_vec[i] = true;
        vmd.true_cycle_slip_vec.push_back(cycleslip);
      }
      if (IsUseMeasurement(freq, MeasType::GRAPHIC, min_alt, tidx)) {
        // Implement GRAPHIC measurement extraction
        vmd.z_meas_vec.push_back((pr_meas + carrier_phase_meas) / 2.0);  // Placeholder
        vmd.meas_freqs_vec.push_back(freq);
        vmd.meas_types_vec.push_back(MeasType::GRAPHIC);
        vmd.prns_vec.push_back(prn);
        vmd.r_ephem_mci_vec.push_back(r_ephem_mci);
        vmd.clk_ephem_vec.push_back(clk_ephem);
        vmd.R_diag_vec.push_back(var_graphic);
        vmd.true_delays_vec.push_back((pr_delay + carrier_delay)
                                      / 2.0);  // Placeholder for true delays
        vmd.min_alt_vec.push_back(min_alt);
        vmd.meas_to_signal_index_vec.push_back(i);
        vmd.use_signal_carrier_vec[i] = true;
        vmd.true_cycle_slip_vec.push_back(cycleslip);
      }
      if (IsUseMeasurement(freq, MeasType::TDCP, min_alt, tidx)) {
        if (tidx > 0) {
          int prev_dbidx = meas_storage_.get_index(tidx - 1, freq, prn);
          if (prev_dbidx == -1) {
            continue;  // No previous measurement for TDCP
          }
          double carrier_phase_prev = meas_storage_.carrier_phase[prev_dbidx];
          double tdcp = carrier_phase_meas - carrier_phase_prev;

          int carrier_int_prev = meas_storage_.carrier_phase_integer[prev_dbidx];
          int carrier_int_curr = meas_storage_.carrier_phase_integer[dbidx];
          double carrier_range_prev = meas_storage_.carrier_phase_range[prev_dbidx];
          double carrier_range_curr = meas_storage_.carrier_phase_range[dbidx];

          // std::cout << "tidx: " << tidx << ", prn: " << prn << ", freq: " <<
          // static_cast<int>(freq)
          //           << ", phase_meas: " << carrier_phase_meas
          //           << ", phase_prev: " << carrier_phase_prev
          //           << ", tdcp: " << tdcp << ", int_curr: " << carrier_int_curr
          //           << ", int_prev: " << carrier_int_prev << ", range_curr: "
          //           << carrier_range_curr << ", range_prev: " << carrier_range_prev << std::endl;

          vmd.z_meas_vec.push_back(tdcp);
          vmd.meas_freqs_vec.push_back(freq);
          vmd.meas_types_vec.push_back(MeasType::TDCP);
          vmd.prns_vec.push_back(prn);
          vmd.r_ephem_mci_vec.push_back(r_ephem_mci);
          vmd.clk_ephem_vec.push_back(clk_ephem);
          // previous timesteps
          double curr_carrier_sigma = meas_storage_.sigma_carrier_phase[dbidx];
          double prev_carrier_sigma = meas_storage_.sigma_carrier_phase[prev_dbidx];
          vmd.r_prev_ephem_mci_vec.push_back(meas_storage_.gnss_rv_ephem_mci[prev_dbidx].head<3>());
          vmd.clk_prev_ephem_vec.push_back(meas_storage_.gnss_clock_ephem[prev_dbidx]);
          double var_tdcp
              = pow(curr_carrier_sigma, 2) + pow(prev_carrier_sigma, 2) + pow(sigma_tdcp, 2);
          vmd.R_diag_vec.push_back(var_tdcp);
          vmd.true_delays_vec.push_back(carrier_delay - meas_storage_.carrier_delay[prev_dbidx]);
          vmd.min_alt_vec.push_back(min_alt);
          vmd.meas_to_signal_index_vec.push_back(i);
          vmd.true_cycle_slip_vec.push_back(cycleslip);
        }
      }
    }
  }

  void StateMeasManager::AddMultiFreqMeasurement(int tidx, State& x_est_full, VecMeasData& vmd) {
    // std::cout << "Generating Iono-Free Measurements at time index " << tidx << std::endl;
    // double f_L1 = 1575.42;  // L1 frequency in MHz
    // double f_L5 = 1176.45;  // L5 frequency in MHz
    // double iono_factor_L1 = f_L1 * f_L1 / (f_L1 * f_L1 - f_L5 * f_L5);
    // double iono_factor_L5 = f_L5 * f_L5 / (f_L1 * f_L1 - f_L5 * f_L5);
    double iono_factor_L1 = 2.2606043275188257;  // Pre-computed iono factors for L1/L5
    double iono_factor_L5 = 1.2606043275188257;  // Pre-computed iono factors for L1/L5
    bool use_ionofree = filter_cfg_.use_ionofree_any;
    bool use_ionofree_pr = filter_cfg_.use_ionofree_pr;
    bool use_ionofree_carrier = filter_cfg_.use_ionofree_carrier;
    Vec3 rx_pos_mci = x_est_full.segment<3>(0);  // Receiver position from state estimate
    double t_tai = t_tais_(tidx);

    // Compute Minimum Satellite Altitude
    double sigma_ure = filter_cfg_.sigma_ure;
    double sigma_iono = filter_cfg_.sigma_iono;
    double sigma_tdcp = filter_cfg_.sigma_tdcp;

    if (!use_ionofree) {
      return;
    }

    // First, find all prns that have dual frequency measurements
    vmd.freq_db_pairs = meas_storage_.GetDualFreqPairs(tidx);
    vmd.use_signal_carrier_ionofree_vec.resize(vmd.freq_db_pairs.size(), false);

    int pair_idx = -1;

    for (const auto& dbidx_pair : vmd.freq_db_pairs) {
      pair_idx++;
      int dbidx_L1 = dbidx_pair.first;
      int dbidx_L5 = dbidx_pair.second;
      int prn_df = meas_storage_.prn[dbidx_L1];

      Vec3d r_ephem_mci_L1 = meas_storage_.gnss_rv_ephem_mci[dbidx_L1].head<3>();
      double clk_ephem_L1 = meas_storage_.gnss_clock_ephem[dbidx_L1];

      if (dbidx_L1 == -1 || dbidx_L5 == -1) {
        // L1 or L5 measurement missing, skip
        continue;
      }

      if (dbidx_L1 != -1 && dbidx_L5 != -1) {
        double pr_L1 = meas_storage_.pseudorange[dbidx_L1];
        double pr_L5 = meas_storage_.pseudorange[dbidx_L5];
        double ionofree_pr = iono_factor_L1 * pr_L1 - iono_factor_L5 * pr_L5;
        double var_pr_L1 = pow(meas_storage_.sigma_pseudorange[dbidx_L1], 2);
        double var_pr_L5 = pow(meas_storage_.sigma_pseudorange[dbidx_L5], 2);
        double var_ionofree_pr = iono_factor_L1 * iono_factor_L1 * var_pr_L1
                                 + iono_factor_L5 * iono_factor_L5 * var_pr_L5 + pow(sigma_ure, 2);

        double cp_L1 = vmd.idx_to_carrier_phase_map[dbidx_L1];
        double cp_L5 = vmd.idx_to_carrier_phase_map[dbidx_L5];
        double ionofree_carrier_phase = iono_factor_L1 * cp_L1 - iono_factor_L5 * cp_L5;
        double var_cp_L1 = pow(meas_storage_.sigma_carrier_phase[dbidx_L1], 2);
        double var_cp_L5 = pow(meas_storage_.sigma_carrier_phase[dbidx_L5], 2);
        double var_ionofree_carrier = iono_factor_L1 * iono_factor_L1 * var_cp_L1
                                      + iono_factor_L5 * iono_factor_L5 * var_cp_L5;

        // Delays
        double pr_delay_L1 = meas_storage_.pr_delay[dbidx_L1];
        double pr_delay_L5 = meas_storage_.pr_delay[dbidx_L5];
        double carrier_delay_L1 = meas_storage_.carrier_delay[dbidx_L1];
        double carrier_delay_L5 = meas_storage_.carrier_delay[dbidx_L5];
        double ionofree_pr_delay = iono_factor_L1 * pr_delay_L1 - iono_factor_L5 * pr_delay_L5;
        double ionofree_carrier_delay
            = iono_factor_L1 * carrier_delay_L1 - iono_factor_L5 * carrier_delay_L5;

        // int cycleslip
        int cycleslip_L1 = meas_storage_.cycle_slip[dbidx_L1];
        int cycleslip_L5 = meas_storage_.cycle_slip[dbidx_L5];
        int cycleslip_L1L5 = (cycleslip_L1 || cycleslip_L5) ? 1 : 0;

        // Min Altitude
        double min_alt = ComputeMinEarthAltitude(t_tai, rx_pos_mci, r_ephem_mci_L1);

        // Pseudorange
        if (IsUseMeasurement(MeasFreq::L1_L5, MeasType::IONOFREE_CODE, min_alt, tidx)) {
          vmd.z_meas_vec.push_back(ionofree_pr);
          vmd.meas_freqs_vec.push_back(MeasFreq::L1_L5);
          vmd.meas_types_vec.push_back(MeasType::IONOFREE_CODE);
          vmd.prns_vec.push_back(prn_df);
          vmd.min_alt_vec.push_back(min_alt);
          vmd.r_ephem_mci_vec.push_back(r_ephem_mci_L1);
          vmd.clk_ephem_vec.push_back(clk_ephem_L1);
          vmd.R_diag_vec.push_back(var_ionofree_pr);
          vmd.true_delays_vec.push_back(ionofree_pr_delay);
          vmd.meas_to_signal_index_vec.push_back(pair_idx);
          vmd.true_cycle_slip_vec.push_back(0);
        }

        // Carrier Phase
        if (IsUseMeasurement(MeasFreq::L1_L5, MeasType::IONOFREE_CARRIER, min_alt, tidx)) {
          // Ambiguity state indices
          vmd.z_meas_vec.push_back(ionofree_carrier_phase);
          vmd.meas_freqs_vec.push_back(MeasFreq::L1_L5);
          vmd.meas_types_vec.push_back(MeasType::IONOFREE_CARRIER);
          vmd.prns_vec.push_back(prn_df);
          vmd.r_ephem_mci_vec.push_back(r_ephem_mci_L1);
          vmd.clk_ephem_vec.push_back(clk_ephem_L1);
          vmd.min_alt_vec.push_back(min_alt);
          vmd.R_diag_vec.push_back(var_ionofree_carrier);
          vmd.true_delays_vec.push_back(ionofree_carrier_delay);
          vmd.use_signal_carrier_ionofree_vec[pair_idx] = true;
          vmd.iono_free_pr_minus_cp.push_back(ionofree_pr - ionofree_carrier_phase);
          vmd.iono_free_pr_minus_cp_sigma.push_back(
              std::sqrt(var_ionofree_pr + var_ionofree_carrier));
          vmd.meas_to_signal_index_vec.push_back(pair_idx);
          vmd.true_cycle_slip_vec.push_back(cycleslip_L1L5);
        }

        // TDCP
        if (IsUseMeasurement(MeasFreq::L1_L5, MeasType::IONOFREE_TDCP, min_alt, tidx)) {
          if (tidx > 0) {
            int prev_dbidx_L1 = meas_storage_.get_index(tidx - 1, MeasFreq::L1, prn_df);
            int prev_dbidx_L5 = meas_storage_.get_index(tidx - 1, MeasFreq::L5, prn_df);
            if (prev_dbidx_L1 == -1 || prev_dbidx_L5 == -1) {
              continue;  // No previous measurement for TDCP
            }
            double cp_L1_prev = meas_storage_.carrier_phase[prev_dbidx_L1];
            double cp_L5_prev = meas_storage_.carrier_phase[prev_dbidx_L5];
            double ionofree_cp_prev = iono_factor_L1 * cp_L1_prev - iono_factor_L5 * cp_L5_prev;
            double ionofree_tdcp = ionofree_carrier_phase - ionofree_cp_prev;
            double sigma_tdcp_ionofree
                = filter_cfg_.sigma_tdcp_ionofree;  // Could be different for iono-free

            double var_cp_L1_prev = pow(meas_storage_.sigma_carrier_phase[prev_dbidx_L1], 2);
            double var_cp_L5_prev = pow(meas_storage_.sigma_carrier_phase[prev_dbidx_L5], 2);
            double var_ionofree_cp_prev = iono_factor_L1 * iono_factor_L1 * var_cp_L1_prev
                                          + iono_factor_L5 * iono_factor_L5 * var_cp_L5_prev;
            double var_ionofree_tdcp
                = var_ionofree_carrier + var_ionofree_cp_prev + pow(sigma_tdcp_ionofree, 2);
            double ionofree_carrier_delay_prev
                = iono_factor_L1 * meas_storage_.carrier_delay[prev_dbidx_L1]
                  - iono_factor_L5 * meas_storage_.carrier_delay[prev_dbidx_L5];

            vmd.z_meas_vec.push_back(ionofree_tdcp);
            vmd.meas_freqs_vec.push_back(MeasFreq::L1_L5);
            vmd.meas_types_vec.push_back(MeasType::IONOFREE_TDCP);
            vmd.prns_vec.push_back(prn_df);
            vmd.r_ephem_mci_vec.push_back(r_ephem_mci_L1);
            vmd.clk_ephem_vec.push_back(clk_ephem_L1);
            vmd.min_alt_vec.push_back(min_alt);
            vmd.R_diag_vec.push_back(var_ionofree_tdcp);
            vmd.true_delays_vec.push_back(ionofree_carrier_delay - ionofree_carrier_delay_prev);
            vmd.r_prev_ephem_mci_vec.push_back(
                meas_storage_.gnss_rv_ephem_mci[prev_dbidx_L1].head<3>());
            vmd.clk_prev_ephem_vec.push_back(meas_storage_.gnss_clock_ephem[prev_dbidx_L1]);
            vmd.meas_to_signal_index_vec.push_back(pair_idx);
            vmd.true_cycle_slip_vec.push_back(cycleslip_L1L5);
          }  // end if tidx > 0
        }  // end if use iono-free TDCP
      }  // end if both dbidx_L1 and dbidx_L5 exist
    }
  }

  VecXd StateMeasManager::UpdateIntegerAmbiguities(int tidx, State& x_est_full, MatXd& P_est_full,
                                                   VecMeasData& vmd) {
    int n_signals = static_cast<int>(vmd.db_indices.size());

    // Count the number of ambiguities -----------------------------------
    int n_ambig = 0;
    // Single Frequency
    for (int i = 0; i < n_signals; i++) {
      if (vmd.use_signal_carrier_vec[i]) {
        n_ambig++;
      }
    }
    // Iono-free
    for (int j = 0; j < static_cast<int>(vmd.freq_db_pairs.size()); j++) {
      if (vmd.use_signal_carrier_ionofree_vec[j]) {
        n_ambig++;
      }
    }

    std::vector<double> N_integer_est_signal(n_ambig, 0.0);
    vmd.ambig_full_state_indices.clear();
    vmd.ambig_est_state_index_for_allmeas_vec.clear();
    bool use_ionofree = filter_cfg_.use_ionofree_any;
    int n_ambig_idx = 0;

    if (filter_cfg_.est_integers) {
      /**********************************************************************
       * Step 1: Update Integer Ambiguities for Single-Frequency Signals
       * *********************************************************************/
      std::vector<FreqPrnPair> curr_pairs = meas_storage_.freqs_prns_at_tidx[tidx];

      for (int i = 0; i < n_signals; i++) {
        if (!vmd.use_signal_carrier_vec[i]) {
          continue;  // Skip if this signal is not used for measurement
        }
        FreqPrnPair curr_pair = curr_pairs[i];
        bool new_acquisition = meas_storage_.acquired_this_step[vmd.db_indices[i]];
        int state_index = filter_state_.freq_prn_to_full_state_index[curr_pair];
        vmd.ambig_full_state_indices.push_back(state_index);

        if (new_acquisition || first_measurement_) {
          N_integer_est_signal[n_ambig_idx]
              = GenerateNewAmbiguityEstimate(vmd.db_indices[i]);        // New acquisition
          x_est_full[state_index] = N_integer_est_signal[n_ambig_idx];  // Update state vector
          P_est_full(state_index, state_index)
              = GenerateNewAmbiguityCovariance(vmd.db_indices[i]);  // Large initial variance
          n_ambig_idx++;
          continue;
        } else {
          N_integer_est_signal[n_ambig_idx]
              = x_est_full[state_index];  // Carry over previous estimate
          n_ambig_idx++;
          continue;
        }
      }

      // Update the number of ambiguities in the filter state
      filter_state_.N_ambiguity_est = static_cast<int>(vmd.ambig_full_state_indices.size());

      // Clear all ambiguity states not in current measurements
      for (int i = 0; i < filter_state_.N_ambiguity_full; i++) {
        int state_index = filter_state_.N_rva + filter_state_.N_clk + filter_state_.N_srp + i;
        if (std::find(vmd.ambig_full_state_indices.begin(), vmd.ambig_full_state_indices.end(),
                      state_index)
            == vmd.ambig_full_state_indices.end()) {
          x_est_full[state_index] = 0.0;
          P_est_full(state_index, state_index) = 0.0;  // Clear variance
        }
      }

      /**********************************************************************
       * Step 2: Update Integer Ambiguities for Multi-Frequency Signals
       * *********************************************************************/
      if (use_ionofree) {
        int pair_size = static_cast<int>(vmd.freq_db_pairs.size());
        int ionofree_idx = 0;

        for (int j = 0; j < pair_size; j++) {
          if (!vmd.use_signal_carrier_ionofree_vec[j]) {
            continue;  // Skip if this signal is not used for iono-free carrier phase
          }
          int dbidx_L1 = vmd.freq_db_pairs[j].first;
          int dbidx_L5 = vmd.freq_db_pairs[j].second;
          bool new_acquisition = meas_storage_.acquired_this_step[dbidx_L1]
                                 || meas_storage_.acquired_this_step[dbidx_L5];
          FreqPrnPair curr_pair = {MeasFreq::L1_L5, meas_storage_.prn[dbidx_L1]};
          int state_index = filter_state_.freq_prn_to_full_state_index[curr_pair];
          vmd.ambig_full_state_indices.push_back(state_index);

          if (new_acquisition || first_measurement_) {
            N_integer_est_signal[n_ambig_idx] = vmd.iono_free_pr_minus_cp[ionofree_idx];
            x_est_full[state_index] = N_integer_est_signal[n_ambig_idx];  // Update state vector
            P_est_full(state_index, state_index)
                = std::pow(vmd.iono_free_pr_minus_cp_sigma[ionofree_idx], 2);
            ionofree_idx++;
            n_ambig_idx++;
            continue;
          } else {
            N_integer_est_signal[n_ambig_idx]
                = x_est_full[state_index];  // Carry over previous estimate
            n_ambig_idx++;
            continue;
          }
        }  // For each frequency-PRN pair
      }  //  if (use_ionofree)

    } else {
      // Not estimating ambiguities
      filter_state_.N_ambiguity_est = 0;
      vmd.ambig_full_state_indices = {};
      N_integer_est_signal = {};
    }

    /**********************************************************************
     * Step 3: Update Number of States
     * *********************************************************************/
    filter_state_.N_ekf_sat_est = filter_state_.N_rva + filter_state_.N_clk + filter_state_.N_srp
                                  + filter_state_.N_ambiguity_est;

    // Get the Integer Ambiguities
    VecXd N_integer_est;
    if (filter_cfg_.est_integers == false) {
      N_integer_est = VecXd::Zero(0);  // Empty if not estimating integers
    } else {
      N_integer_est = VecXd::Zero(static_cast<int>(N_integer_est_signal.size()));
      for (size_t i = 0; i < N_integer_est_signal.size(); i++) {
        N_integer_est(static_cast<int>(i)) = N_integer_est_signal[i];
      }
    }

    // Store the corresponding ambiguity state index (in estimated state vector) for each
    // measurement
    for (size_t i = 0; i < vmd.z_meas_vec.size(); i++) {
      // Single Frequency Measureemnt
      MeasType meas_type = vmd.meas_types_vec[i];
      MeasFreq meas_freq = vmd.meas_freqs_vec[i];
      int prn = vmd.prns_vec[i];
      FreqPrnPair freq_prn = {meas_freq, prn};

      if (meas_type == MeasType::CARRIER || meas_type == MeasType::GRAPHIC
          || meas_type == MeasType::IONOFREE_CARRIER) {
        int signal_idx = vmd.meas_to_signal_index_vec[i];
        int full_idx = filter_state_.freq_prn_to_full_state_index[freq_prn];
        int est_idx = FullIndexToEstIndex(full_idx, vmd.ambig_full_state_indices);
        vmd.ambig_est_state_index_for_allmeas_vec.push_back(est_idx);
      } else {
        vmd.ambig_est_state_index_for_allmeas_vec.push_back(-1);  // No ambiguity associated
      }
    }

    return N_integer_est;
  }

  MeasData StateMeasManager::GetMeasurement(State& x_est_full, MatXd& P_est_full, int tidx) {
    // Placeholder for getting measurements based on the current state estimate
    // Implement the logic to extract or compute measurements from x_est at time index tidx

    // First extract the coorrect measurements from meas_storage_
    VecMeasData vmd;
    vmd.db_indices = meas_storage_.get_index_at_tidx(tidx);
    int n_signals = static_cast<int>(
        vmd.db_indices.size());  // number of measurements (including multiple frequencies)

    if (n_signals == 0) {
      // No measurements at this time index
      filter_state_.N_ambiguity_est = 0;  // No ambiguities estimated
      return CreateEmptyMeasData(tidx);
    }

    // Single frequency measurements ----------------------------------------------
    AddSingleFreqMeasurement(tidx, x_est_full, vmd);

    // Iono-free Combination --------------------------------
    AddMultiFreqMeasurement(tidx, x_est_full, vmd);

    // Extract integer ambiguities from state vector
    // -------------------------------------------------------
    VecXd N_integer_est;
    N_integer_est = UpdateIntegerAmbiguities(tidx, x_est_full, P_est_full, vmd);

    // Convert z_meas_vec to VecXd --------------------------------
    VecXd z_meas = VecXd::Zero(static_cast<int>(vmd.z_meas_vec.size()));
    VecXd R_diag = VecXd::Zero(static_cast<int>(vmd.R_diag_vec.size()));
    VecXd true_delays = VecXd::Zero(static_cast<int>(vmd.true_delays_vec.size()));
    VecXd min_altitudes = VecXd::Zero(static_cast<int>(vmd.min_alt_vec.size()));
    for (size_t i = 0; i < vmd.z_meas_vec.size(); i++) {
      z_meas(static_cast<int>(i)) = vmd.z_meas_vec[i];
      R_diag(static_cast<int>(i)) = vmd.R_diag_vec[i];
      true_delays(static_cast<int>(i)) = vmd.true_delays_vec[i];
      min_altitudes(static_cast<int>(i)) = vmd.min_alt_vec[i];
    }

    int n_meas = static_cast<int>(vmd.z_meas_vec.size());

    // Add to MeasData structure
    MeasData meas_data;
    meas_data.tidx = tidx;
    meas_data.tspan = meas_storage_.tspan_at_tidx[tidx];
    meas_data.t_tai = meas_storage_.tai_at_tidx[tidx];
    meas_data.n_meas = static_cast<int>(vmd.z_meas_vec.size());
    meas_data.z = z_meas;
    meas_data.R_diag = R_diag;
    meas_data.N_integer_est = N_integer_est;
    meas_data.ambig_est_state_index_for_allmeas = vmd.ambig_est_state_index_for_allmeas_vec;
    meas_data.ambig_full_state_indices = vmd.ambig_full_state_indices;
    meas_data.prns = vmd.prns_vec;
    meas_data.types = vmd.meas_types_vec;
    meas_data.freqs = vmd.meas_freqs_vec;
    meas_data.r_ephem_mci = vmd.r_ephem_mci_vec;
    meas_data.clock_ephem = vmd.clk_ephem_vec;
    meas_data.r_prev_ephem_mci = vmd.r_prev_ephem_mci_vec;
    meas_data.clock_prev_ephem = vmd.clk_prev_ephem_vec;
    meas_data.true_delays = true_delays;
    meas_data.est_min_alt = min_altitudes;
    meas_data.is_outlier.resize(n_meas, false);
    meas_data.cycle_slip_correction.resize(n_meas, 0.0);
    meas_data.true_cycle_slips = vmd.true_cycle_slip_vec;
    meas_data.outlier_types.resize(n_meas, OUTLIER_TYPE::NO_OUTLIER);
    meas_data.outlier_thresholds.resize(n_meas, 0.0);
    if (tidx > 0) {
      meas_data.tspan_prev = meas_storage_.tspan_at_tidx[tidx - 1];
      meas_data.t_tai_prev = meas_storage_.tai_at_tidx[tidx - 1];
    } else {
      meas_data.tspan_prev = 0.0;
      meas_data.t_tai_prev = 0.0;
    }

    first_measurement_ = false;  // Clear first measurement flag after first call

    return meas_data;
  }

  double StateMeasManager::ComputeMinEarthAltitude(double t_tai, const Vec3& rx_pos_mci,
                                                   const Vec3& tx_pos_mci) {
    // Convert Moon-centered inertial -> Earth-centered inertial
    double t_tai_rx = t_tai;
    double t_tai_tx = t_tai_rx - (rx_pos_mci - tx_pos_mci).norm().val() / C;

    Vec3 rx_pos_eci = ConvertFrame(Real(t_tai_rx), rx_pos_mci, Frame::MOON_CI, Frame::ECEF, false);
    Vec3 tx_pos_eci = ConvertFrame(Real(t_tai_tx), tx_pos_mci, Frame::MOON_CI, Frame::ECEF, false);

    // Segment from transmitter to receiver in ECI
    Vec3 d = rx_pos_eci - tx_pos_eci;  // direction along the segment
    double d2 = d.squaredNorm().val();

    // Earth center is at origin in ECI, so minimize || tx + s*d || over s in [0,1]
    double s = 0.0;
    if (d2 > 0.0) {
      // Unclamped minimizer for point-to-line distance (origin to line)
      s = -tx_pos_eci.dot(d) / d2;

      // Clamp to the segment [tx, rx]
      if (s < 0.0) s = 0.0;
      if (s > 1.0) s = 1.0;
    } else {
      // Degenerate segment: tx == rx
      s = 0.0;
    }

    Vec3 closest = tx_pos_eci + s * d;
    double min_dist_to_center_m = closest.norm().val();

    // Earth radius (meters) — replace with your project’s constant/getter
    // Examples you might have: Constants::R_EARTH, BodyConst(Body::EARTH).radius_m, etc.)
    const double min_altitude_m = min_dist_to_center_m - R_EARTH;
    return min_altitude_m;  // can be negative if the segment intersects Earth
  }

  std::string OutlierToString(OUTLIER_TYPE outlier_type) {
    switch (outlier_type) {
      case OUTLIER_TYPE::NO_OUTLIER: return "None";
      case OUTLIER_TYPE::LARGE_RESIDUAL: return "Residual";
      case OUTLIER_TYPE::UNRELIABLE_STATE: return "Unreliable";
      case OUTLIER_TYPE::CYCLE_SLIP_DISCARDED: return "CS Discarded";
      case OUTLIER_TYPE::CYCLE_SLIP_CORRECTED: return "CS Corrected";
      default: return "Unknown";
    }
  }

  void StateMeasManager::PrintMeasData(const MeasData& meas_data, const VecXd& z_pred) const {
    // Basic Info
    VecXd offset_vec = VecXd::Constant(meas_data.N_integer_est.size(), AMBIGUITY_OFFSET_N);
    VecXd N_est_vec = meas_data.N_integer_est;

    std::cout << "[Measurements]" << std::endl;
    std::cout << "tspan (s)                  : " << meas_data.tspan << std::endl;
    std::cout << "TAI Time (s)               : " << meas_data.t_tai << std::endl;
    if (filter_cfg_.est_integers) {
      std::cout << "Estimated Integer Ambiguities: " << (N_est_vec + offset_vec).transpose()
                << std::endl;
      std::cout << "State Indices of Ambiguities   : ";
      for (const auto& idx : meas_data.ambig_full_state_indices) {
        std::cout << idx << " ";
      }
      std::cout << std::endl;
      std::cout << "State Indices of Estimated Ambiguities for All Measurements: ";
      for (const auto& idx : meas_data.ambig_est_state_index_for_allmeas) {
        std::cout << idx << " ";
      }
      std::cout << std::endl;
    }

    // Variables
    VecXd dz = meas_data.z - z_pred;

    // Table
    std::cout << std::string(200, '-') << std::endl;
    std::cout << std::left << std::setw(6) << "Index" << std::setw(6) << "PRN" << std::setw(10)
              << "Freq" << std::setw(20) << "Type" << std::setw(15) << "Measurement"
              << std::setw(15) << "Predicted" << std::setw(15) << "Residual" << std::setw(15)
              << "Delays" << std::setw(15) << "Std Dev" << std::setw(15) << "Threshold"
              << std::setw(16) << "Min Altitude" << std::setw(12) << "Outlier" << std::setw(15)
              << "True Slip" << std::setw(15) << "Correction" << std::setw(18) << "Outlier Type"
              << std::endl;

    std::cout << std::string(200, '-') << std::endl;

    // Rows
    for (size_t i = 0; i < meas_data.z.size(); i++) {
      std::cout << std::left << std::setw(6) << i << std::setw(6) << meas_data.prns[i]
                << std::setw(10) << enum_name(meas_data.freqs[i]) << std::setw(20)
                << enum_name(meas_data.types[i]) << std::setw(15) << std::scientific
                << meas_data.z(i) << std::setw(15) << std::scientific << z_pred(i) << std::setw(15)
                << std::scientific << dz(i) << std::setw(15) << std::scientific
                << meas_data.true_delays(i) << std::setw(15) << std::scientific
                << std::sqrt(meas_data.R_diag(i)) << std::setw(15) << std::scientific
                << meas_data.outlier_thresholds[i] << std::setw(16) << std::fixed
                << meas_data.est_min_alt(i) / 1000 << std::setw(12)
                << (meas_data.is_outlier[i] ? "Yes" : "No") << std::setw(15)
                << meas_data.true_cycle_slips[i] << std::setw(15)
                << meas_data.cycle_slip_correction[i] << std::setw(18)
                << OutlierToString(meas_data.outlier_types[i]) << std::endl;
    }

    std::cout << std::string(200, '-') << std::endl;
    std::cout << " " << std::endl;
  }

  void StateMeasManager::PrintSingleMeasData(const MeasData& meas_data, const VecXd& z_pred,
                                             int meas_idx) const {
    // Basic Info
    if (meas_idx == 0) {
      // headers
      std::cout << "[Measurements]" << std::endl;
      std::cout << "tspan (s)                  : " << meas_data.tspan << std::endl;
      std::cout << "TAI Time (s)               : " << meas_data.t_tai << std::endl;
      if (filter_cfg_.est_integers) {
        VecXd offset_vec = VecXd::Constant(meas_data.N_integer_est.size(), AMBIGUITY_OFFSET_N);
        std::cout << "Estimated Integer Ambiguities: "
                  << meas_data.N_integer_est.transpose() + offset_vec.transpose() << std::endl;
        std::cout << "State Indices of Ambiguities   : ";
        for (const auto& idx : meas_data.ambig_full_state_indices) {
          std::cout << idx << " ";
        }
        std::cout << std::endl;
        std::cout << "State Indices of Estimated Ambiguities for All Measurements: ";
        for (const auto& idx : meas_data.ambig_est_state_index_for_allmeas) {
          std::cout << idx << " ";
        }
        std::cout << std::endl;
      }

      // Table
      std::cout << std::string(160, '-') << std::endl;
      std::cout << std::left << std::setw(6) << "Index" << std::setw(6) << "PRN" << std::setw(10)
                << "Freq" << std::setw(20) << "Type" << std::setw(18) << "Measurement"
                << std::setw(18) << "Predicted" << std::setw(18) << "Residual" << std::setw(18)
                << "Delays" << std::setw(18) << "Std Dev" << std::setw(15) << "Min Altitude"
                << std::endl;

      std::cout << std::string(160, '-') << std::endl;
    }

    double dz = meas_data.z(meas_idx) - z_pred(0);

    // Rows
    std::cout << std::left << std::setw(6) << meas_idx << std::setw(6) << meas_data.prns[meas_idx]
              << std::setw(10) << enum_name(meas_data.freqs[meas_idx]) << std::setw(20)
              << enum_name(meas_data.types[meas_idx]) << std::setw(18) << std::scientific
              << meas_data.z(meas_idx) << std::setw(18) << std::scientific << z_pred(0)
              << std::setw(18) << std::scientific << dz << std::setw(18) << std::scientific
              << meas_data.true_delays(meas_idx) << std::setw(15) << std::scientific
              << std::sqrt(meas_data.R_diag(meas_idx)) << std::setw(12) << std::fixed
              << meas_data.est_min_alt(meas_idx) / 1000 << std::endl;

    if (meas_idx == meas_data.z.size() - 1) {
      std::cout << std::string(160, '-') << std::endl;
      std::cout << " " << std::endl;
    }
  }

  VecX StateMeasManager::ComputePredictedMeasurement(const VecX& x_state, const MeasData& meas_data,
                                                     MatXd* H, int meas_idx) {
    // Placeholder for computing predicted measurement based on the current state
    // Implement the logic to compute the predicted measurement for meas_index

    int n_meas = meas_data.z.size();
    if (meas_idx >= 0 && meas_idx < n_meas) {
      // Compute predicted measurement for a specific measurement index
      n_meas = 1;
    }
    std::vector<int> meas_indices;
    if (meas_idx >= 0 && meas_idx < meas_data.z.size()) {
      meas_indices.push_back(meas_idx);
    } else {
      meas_indices.resize(n_meas);
      for (int i = 0; i < n_meas; i++) {
        meas_indices[i] = i;
      }
    }

    VecX y_pred = VecX::Zero(n_meas);

    bool compute_H = (H != nullptr);

    if (compute_H) {
      H->resize(n_meas, x_state.size());
      H->setZero();
    }

    // Extract parameters
    Vec6 rx_posvel_mci = x_state.segment<6>(0);  // Receiver position and velocity in MCI
    Vec3 rx_pos_mci = x_state.segment<3>(0);     // Receiver position in MCI
    Real rx_clk = x_state(filter_state_.N_rva);  // Receiver clock bias

    // Extract necessary state components
    int tdcp_idx = 0;
    Real relativistic_delay = 0.0;
    if (!filter_cfg_.joint_orbit_clock_dynamics) {
      relativistic_delay
          = ComputeRelativisticDelayLT(rx_posvel_mci, Frame::MOON_CI, meas_data.tspan) * C;
    }

    double lamda_L1 = get_wavelength(MeasFreq::L1);
    double lamda_L5 = get_wavelength(MeasFreq::L5);
    double f_L1 = 1575.42;  // L1 frequency in MHz
    double f_L5 = 1176.45;  // L5 frequency in MHz
    double iono_factor_L1 = f_L1 * f_L1 / (f_L1 * f_L1 - f_L5 * f_L5);
    double iono_factor_L5 = f_L5 * f_L5 / (f_L1 * f_L1 - f_L5 * f_L5);

    for (int ii = 0; ii < n_meas; ii++) {
      int i = meas_indices[ii];
      Vec3d tx_pos_mci = meas_data.r_ephem_mci[i];
      double clk_tx = meas_data.clock_ephem[i];
      Real ambiguity_range = 0.0;
      MeasType meas_type = meas_data.types[i];

      Real shapiro_delay
          = ComputeShapiroDelay(meas_data.t_tai, tx_pos_mci, rx_pos_mci, Frame::MOON_CI) * C;
      Real total_delay = shapiro_delay + relativistic_delay;

      // Debugging output
      switch (meas_type) {
        case MeasType::CODE:
        case MeasType::IONOFREE_CODE: {
          y_pred(ii) = (rx_pos_mci - tx_pos_mci).norm() + (rx_clk - clk_tx) + total_delay;
          if (compute_H) {
            H->row(ii).segment<3>(0) = (rx_pos_mci - tx_pos_mci).normalized().transpose();
            H->row(ii)(filter_state_.N_rva) = 1.0;
          }

          break;
        }
        case MeasType::CARRIER: {
          int ambig_est_idx = meas_data.ambig_est_state_index_for_allmeas[i];
          double lamda = get_wavelength(meas_data.freqs[i]);
          ambiguity_range = (x_state(ambig_est_idx) + AMBIGUITY_OFFSET_N) * lamda;
          y_pred(ii) = (rx_pos_mci - tx_pos_mci).norm() + (rx_clk - clk_tx) - ambiguity_range
                       + total_delay;
          if (compute_H) {
            H->row(ii).segment<3>(0) = (rx_pos_mci - tx_pos_mci).normalized().transpose();
            H->row(ii)(filter_state_.N_rva) = 1.0;
            H->row(ii)(ambig_est_idx) = -lamda;
          }
          break;
        }
        case MeasType::GRAPHIC: {
          int ambig_est_idx = meas_data.ambig_est_state_index_for_allmeas[i];
          double lamda = get_wavelength(meas_data.freqs[i]);
          ambiguity_range = (x_state(ambig_est_idx) + AMBIGUITY_OFFSET_N) * lamda;
          Real pr_pred = (rx_pos_mci - tx_pos_mci).norm() + (rx_clk - clk_tx);
          Real carrier_phase_m_pred = pr_pred - ambiguity_range;
          y_pred(ii) = (pr_pred + carrier_phase_m_pred) / 2.0 + total_delay;
          if (compute_H) {
            H->row(ii).segment<3>(0) = (rx_pos_mci - tx_pos_mci).normalized().transpose();
            H->row(ii)(filter_state_.N_rva) = 1.0;
            H->row(ii)(ambig_est_idx) = -0.5 * lamda;
          }
          break;
        }
        case MeasType::TDCP:
        case MeasType::IONOFREE_TDCP: {
          VecXd tx_pos_mci_prev = meas_data.r_prev_ephem_mci[tdcp_idx];
          double clk_tx_prev = meas_data.clock_prev_ephem[tdcp_idx];
          tdcp_idx++;
          Real carrier_curr_m = (rx_pos_mci - tx_pos_mci).norm() + (rx_clk - clk_tx) + total_delay;

          VecX rx_pos_prev
              = x_state.segment<3>(filter_state_.N_ekf_sat_est);  // Assuming static for simplicity
          VecX rx_posvel_prev
              = x_state.segment<6>(filter_state_.N_ekf_sat_est);  // Assuming static for simplicity
          Real rx_clk_prev = x_state(filter_state_.N_ekf_sat_est
                                     + filter_state_.N_rva);  // Assuming static for simplicity

          // Retrieve previous receiver position and clock from augmented state
          Real shapiro_delay_prev = ComputeShapiroDelay(meas_data.t_tai_prev, tx_pos_mci_prev,
                                                        rx_pos_prev, Frame::MOON_CI)
                                    * C;

          Real relativistic_delay_prev = 0.0;

          if (!filter_cfg_.joint_orbit_clock_dynamics) {
            relativistic_delay_prev
                = ComputeRelativisticDelayLT(rx_posvel_prev, Frame::MOON_CI, meas_data.tspan_prev)
                  * C;
          }

          Real total_delay_prev = shapiro_delay_prev + relativistic_delay_prev;
          Real carrier_prev_m = (rx_pos_prev - tx_pos_mci_prev).norm() + (rx_clk_prev - clk_tx_prev)
                                + total_delay_prev;

          if (compute_H) {
            H->row(ii).segment<3>(0) = (rx_pos_mci - tx_pos_mci).normalized().transpose();
            H->row(ii)(filter_state_.N_rva) = 1.0;
            H->row(ii).segment<3>(filter_state_.N_ekf_sat_est)
                = -(rx_pos_prev - tx_pos_mci_prev).normalized().transpose();
            H->row(ii)(filter_state_.N_ekf_sat_est + filter_state_.N_rva) = -1.0;
          }

          y_pred(ii) = carrier_curr_m - carrier_prev_m;

          // std::cout << "N_ekf_sat_est: " << filter_state_.N_ekf_sat_est << std::endl;
          // std::cout << "TDCP idx: " << tdcp_idx << std::endl;
          // std::cout << "[Current]" << std::endl;
          // std::cout << "Carrier Current: " << carrier_curr << std::endl;
          // std::cout << " Rx Position Current: " << rx_pos_mci.transpose() << std::endl;
          // std::cout << " Rx Clock Current: " << rx_clk << std::endl;
          // std::cout << " Tx Position Current: " << tx_pos_mci.transpose() << std::endl;
          // std::cout << " Tx Clock Current: " << clk_tx << std::endl;
          // std::cout << " " << std::endl;
          // std::cout << "[Previous]" << std::endl;
          // std::cout << "Carrier Previous: " << carrier_prev << std::endl;
          // std::cout<< "  Rx Position Previous: " << rx_pos_prev.transpose() << std::endl;
          // std::cout<< "  Rx Position Diff: " << (rx_pos_mci - rx_pos_prev).transpose() <<
          // std::endl; std::cout<< "  Rx Clock Previous: " << rx_clk_prev << std::endl; std::cout<<
          // "  Tx Position Previous: " << tx_pos_mci_prev.transpose() << std::endl; std::cout<< "
          // Tx Clock Previous: " << clk_tx_prev << std::endl; std::cout<< "Carrier Difference: " <<
          // y_pred(ii) << std::endl; std::cout << " " << std::endl; std::cout << " --------  " <<
          // std::endl;

          break;
        }
        case MeasType::IONOFREE_CARRIER: {
          Real pr = (rx_pos_mci - tx_pos_mci).norm() + (rx_clk - clk_tx) + total_delay;
          int ambig_est_idx_L1L5 = meas_data.ambig_est_state_index_for_allmeas[i];
          Real ambiguity_range_L1L5 = (x_state(ambig_est_idx_L1L5) + AMBIGUITY_OFFSET_N);
          y_pred(ii) = pr - ambiguity_range_L1L5;

          if (compute_H) {
            H->row(ii).segment<3>(0) = (rx_pos_mci - tx_pos_mci).normalized().transpose();
            H->row(ii)(filter_state_.N_rva) = 1.0;
            H->row(ii)(ambig_est_idx_L1L5) = -1.0;
          }
          break;
        }
        default:
          throw std::runtime_error("Unknown measurement type in ComputePredictedMeasurement.");
      }
    }

    return y_pred;
  }

  void StateMeasManager::DetectCycleSlip(MeasData& meas_data, const VecX& dz, const MatXd* S,
                                         int meas_idx) {
    int n_meas = meas_data.z.size();
    std::vector<int> meas_indices;
    if (meas_idx >= 0 && meas_idx < n_meas) {
      n_meas = 1;
      meas_indices.push_back(meas_idx);
    } else {
      meas_indices.resize(n_meas);
      for (int i = 0; i < n_meas; i++) {
        meas_indices[i] = i;
      }
    }

    for (int i = 0; i < n_meas; i++) {
      int idx = meas_indices[i];
      MeasType meas_type = meas_data.types[idx];

      double residual = dz(i);
      double innovation_var = (*S)(i, i);
      double lambda_signal = get_wavelength(meas_data.freqs[idx]);

      bool is_carrier = true;
      double detect_sigma = 0.0;

      if (meas_type != MeasType::CARRIER && meas_type != MeasType::GRAPHIC
          && meas_type != MeasType::IONOFREE_CARRIER && meas_type != MeasType::TDCP
          && meas_type != MeasType::IONOFREE_TDCP) {
        // Carrier related
        is_carrier = false;
        detect_sigma = 3.0;  // Use a fixed sigma threshold for non-carrier measurements
      } else {
        is_carrier = true;
        detect_sigma = filter_cfg_.cycle_slip_detect_sigma;
      }

      double threshold = detect_sigma * std::sqrt(innovation_var);  // e.g., 4-sigma threshold
      meas_data.outlier_thresholds[idx] = threshold;

      if (!is_carrier) {
        // Just check if residual is within reasonable bounds for non-carrier measurements
        if (std::abs(residual) >= threshold) {
          meas_data.is_outlier[idx] = true;
          meas_data.outlier_types[idx] = OUTLIER_TYPE::LARGE_RESIDUAL;
        } else {
          meas_data.is_outlier[idx] = false;
          meas_data.outlier_types[idx] = OUTLIER_TYPE::NO_OUTLIER;
        }
        meas_data.cycle_slip_correction[idx] = 0.0;
        continue;  // Only check cycle slips for carrier-related measurements
      }

      // Cycle Slip Detection Logic for Carrier Measurements
      if (threshold >= lambda_signal) {
        // Lambda < Threshold  -> unreliable state estimate for residual check
        meas_data.is_outlier[idx] = true;  // Disable cycle slip detection
        meas_data.outlier_types[idx] = OUTLIER_TYPE::UNRELIABLE_STATE;
        meas_data.cycle_slip_correction[idx] = 0.0;
      } else if ((std::abs(residual) >= lambda_signal) && (threshold < lambda_signal)) {
        // Threshold < Lambda < Residual
        meas_data.is_outlier[idx] = true;

        // Check if we can coorect for the signal wavelength
        int sign = (residual > 0) ? 1 : -1;
        double z_correct = 0.0;

        if (std::abs(residual) > lambda_signal) {
          // Large residual, correct by multiple wavelengths
          int n_cycles = static_cast<int>(std::round(std::abs(residual) / lambda_signal));
          z_correct = -sign * n_cycles * lambda_signal;
          residual = residual + z_correct;
        } else {
          // Small residual, correct by one wavelength
          z_correct = -sign * lambda_signal;
          residual = residual + z_correct;
        }

        // Check if correction is successful
        if (std::abs(residual) < lambda_signal && std::abs(residual) < threshold
            && filter_cfg_.fix_cycleslip) {
          // Successful correction
          // Update the measurement to reflect the correction
          meas_data.cycle_slip_correction[idx] = z_correct;
          meas_data.is_outlier[idx] = true;
          meas_data.outlier_types[idx] = OUTLIER_TYPE::CYCLE_SLIP_CORRECTED;
        } else {
          // Correction unsuccessful, likely due to large state uncertainty -> Do not use
          // measurement
          meas_data.cycle_slip_correction[idx] = 0.0;
          meas_data.is_outlier[idx] = true;
          meas_data.outlier_types[idx] = OUTLIER_TYPE::CYCLE_SLIP_DISCARDED;
        }
      } else if (std::abs(residual) > threshold) {
        // Threshold < Residual < Lambda -> simply mark as outlier, no correction
        meas_data.is_outlier[idx] = true;
        meas_data.cycle_slip_correction[idx] = 0.0;
        meas_data.outlier_types[idx] = OUTLIER_TYPE::LARGE_RESIDUAL;
      } else {
        // residual < threshold < Lambda -> safe measurement
        meas_data.is_outlier[idx] = false;
        meas_data.cycle_slip_correction[idx] = 0.0;
        meas_data.outlier_types[idx] = OUTLIER_TYPE::NO_OUTLIER;
      }
    }  // For each measurement
  }

  void StateMeasManager::RemoveCorrectCycleSlip(MeasData& meas_data, KalmanFilter* filter,
                                                int meas_idx) {
    int n_meas = meas_data.z.size();
    std::vector<int> meas_indices;
    if (meas_idx >= 0 && meas_idx < n_meas) {
      n_meas = 1;
      meas_indices.push_back(meas_idx);
    } else {
      meas_indices.resize(n_meas);
      for (int i = 0; i < n_meas; i++) {
        meas_indices[i] = i;
      }
    }

    std::vector<bool> use_meas(n_meas, true);
    int n_use_meas = 0;

    VecXd dz = filter->GetMeasurementResidual();

    for (int i = 0; i < n_meas; i++) {
      int idx = meas_indices[i];
      if (meas_data.is_outlier[idx]) {
        if (meas_data.cycle_slip_correction[idx] != 0.0) {
          // Correction applied -> just change meas_data.z and leave H, R unchanged
          dz(i) += meas_data.cycle_slip_correction[idx];
          use_meas[i] = true;
          n_use_meas = n_use_meas + 1;  // Keep this measurement
        } else {
          use_meas[i] = false;  // Remove this measurement
        }
      } else {
        // No cycle slip detected -> keep measurement
        use_meas[i] = true;
        n_use_meas = n_use_meas + 1;
      }
    }

    // Remove dz
    MatXd H = filter->GetMeasurementJacobian();
    MatXd R = filter->GetMeasurementNoiseCov();
    MatXd S = filter->GetInnovationCov();

    VecXd dz_new = VecXd::Zero(n_use_meas);
    MatXd H_new = MatXd::Zero(n_use_meas, H.cols());
    MatXd R_new = MatXd::Zero(n_use_meas, n_use_meas);
    MatXd S_new = MatXd::Zero(n_use_meas, n_use_meas);

    int n = 0;
    for (int i = 0; i < n_meas; i++) {
      if (!use_meas[i]) continue;
      dz_new(n) = dz(i);
      H_new.row(n) = H.row(i);

      int m = 0;
      for (int j = 0; j < n_meas; j++) {
        if (!use_meas[j]) continue;
        R_new(n, m) = R(i, j);
        S_new(n, m) = S(i, j);
        m++;
      }
      n++;
    }

    // Update
    filter->SetMeasurementResidual(dz_new);
    filter->SetMeasurementJacobian(H_new);
    filter->SetMeasurementNoiseCov(R_new);
    filter->SetInnovationCov(S_new);
  }

  FaultDetectionFunction StateMeasManager::CreateFaultDetectionFunction(MeasData& meas_data,
                                                                        int meas_idx,
                                                                        bool verbose) {
    // Placeholder for creating fault detection function
    // Implement the logic to create fault detection based on meas_data and meas_idx
    FaultDetectionFunction fault_detection_func
        = [this, &meas_data, meas_idx, verbose](KalmanFilter* filter) {
            MatXd P = filter->GetCovariance();
            VecXd dz = filter->GetMeasurementResidual();
            MatXd S = filter->GetInnovationCov();
            MatXd H = filter->GetMeasurementJacobian();
            MatXd R = filter->GetMeasurementNoiseCov();

            // First detect cycle-slips and record that to meas_data
            this->DetectCycleSlip(meas_data, dz, &S, meas_idx);

            // Print Measurement Data
            VecXd z_pred = filter->GetPredictedMeasurement();

            // If verbose, print measurement data
            if (verbose) {
              if (meas_idx != -1) {
                this->PrintSingleMeasData(meas_data, z_pred, meas_idx);
              } else {
                this->PrintMeasData(meas_data, z_pred);
              }
            }

            // Then use that information to remove or correct cycle-slipped measurements
            // Inside it fixes dz, H, R accordingly
            this->RemoveCorrectCycleSlip(meas_data, filter, meas_idx);
          };

    return fault_detection_func;
  }

  FilterMeasurementFunction StateMeasManager::CreateMeasurementFunction(MeasData& meas_data,
                                                                        bool verbose) {
    bool use_ad = filter_cfg_.use_meas_autodiff;

    // Construct Measurement Function
    FilterMeasurementFunction meas_func
        = [this, verbose, use_ad, meas_data](const State& x_state, MatXd* H, MatXd* R) -> VecXd {
      VecX z_pred = VecX::Zero(meas_data.z.size());

      if (use_ad) {
        // This lambda captures the current StateMeasManager and meas_data
        // First construct a function that returns predicted measurements based on x_state
        std::function<VecX(const VecX& x)> meas_func = [this, &meas_data](const VecX& x) -> VecX {
          return this->ComputePredictedMeasurement(x, meas_data, nullptr);
        };

        VecX x0_tmp = x_state.cast<double>();  // Refresh all the jacobian chains

        // Use AutoDiff to compute H matrix
        jacobian(meas_func, wrt(x0_tmp), at(x0_tmp), z_pred, *H);

      } else {
        z_pred = ComputePredictedMeasurement(x_state, meas_data, H);
      }

      *R = meas_data.R_diag.asDiagonal();

      return z_pred.cast<double>();
    };

    return meas_func;
  }

  FilterMeasurementFunction StateMeasManager::CreateSingleMeasurementFunction(MeasData& meas_data,
                                                                              int meas_idx,
                                                                              bool verbose) {
    bool use_ad = filter_cfg_.use_meas_autodiff;

    // Construct Measurement Function
    FilterMeasurementFunction meas_func = [this, verbose, &meas_data, use_ad, meas_idx](
                                              const State& x_state, MatXd* H, MatXd* R) -> VecXd {
      VecX z_pred = VecX::Zero(1);

      if (use_ad) {
        // This lambda captures the current StateMeasManager and meas_data
        // First construct a function that returns predicted measurements based on x_state
        std::function<VecX(const VecX& x)> meas_func
            = [this, &meas_data, meas_idx](const VecX& x) -> VecX {
          return this->ComputePredictedMeasurement(x, meas_data, nullptr, meas_idx);
        };

        VecX x0_tmp = x_state.cast<double>();  // Refresh all the jacobian chains

        // Use AutoDiff to compute H matrix
        jacobian(meas_func, wrt(x0_tmp), at(x0_tmp), z_pred, *H);

      } else {
        z_pred = ComputePredictedMeasurement(x_state, meas_data, H, meas_idx);
      }

      *R = MatXd::Zero(1, 1);
      (*R)(0, 0) = meas_data.R_diag(meas_idx);

      return z_pred.cast<double>();
    };

    return meas_func;
  }

  FilterDynamicsFunction StateMeasManager::CreateDynamicsFunction(int N_ambiguity, bool verbose) {
    FilterDynamicsFunction dyn_func
        = [this, N_ambiguity, verbose](const State& x, Real t0, Real tf, const State* u,
                                       MatXd* F) -> State {
      // This lambda captures the current StateMeasManager
      // Implement the dynamics function logic here

      int n_state = x.size() - N_ambiguity;
      MatXd F_short = MatXd::Zero(n_state, n_state);

      // First, simply call the joint state dynamics function for rva
      FilterDynamicsFunction joint_dyn_func = this->joint_state_.GetDynamicsFunction();
      State x_next_short = joint_dyn_func(x, t0, tf, u, &F_short);

      // Next State
      State x_next = State::Zero(x.size());
      x_next.head(n_state) = x_next_short;

      // Extend the state transition matrix F to include ambiguities
      F->resize(x.size(), x.size());
      F->setZero(x.size(), x.size());
      F->topLeftCorner(n_state, n_state) = F_short;

      if (N_ambiguity > 0) {
        Real dt = tf - t0;
        double tau_ambiguity = filter_cfg_.tau_ambiguity;
        double q_ambiguity = filter_cfg_.q_ambiguity;
        double exp_factor = exp(-dt / tau_ambiguity);

        x_next.tail(N_ambiguity) = exp_factor * x.tail(N_ambiguity);  // Ambiguities remain constant
        F->bottomRightCorner(N_ambiguity, N_ambiguity)
            = exp_factor * MatXd::Identity(N_ambiguity, N_ambiguity);
      }

      if (verbose) {
        auto old_flags = std::cout.flags();
        auto old_precision = std::cout.precision();
        std::cout << "[Dynamics Model] " << std::endl;
        std::cout << std::scientific << std::setprecision(3);
        std::cout << "Dynamics State Transition x_next: \n" << x_next.transpose() << std::endl;
        std::cout << "Dynamics State Transition F: \n" << *F << std::endl;
        std::cout << " " << std::endl;
        std::cout.flags(old_flags);
        std::cout.precision(old_precision);
      }

      return x_next;
    };

    return dyn_func;
  }

  void StateMeasManager::UpdateProcessNoise(Ptr<KalmanFilter> filter, const VecXd& z) {
    int n_meas = z.size();

    // Adaptive Process Noise
    if (filter_cfg_.use_adaptive_process_noise) {
      if (n_meas > 0) {
        VecXd dx = filter->GetStateCorrection().head(6);
        MatXd Sigma_dx = filter->GetStateCorrectionCov().topLeftCorner(6, 6);
        MatXd P_bar = filter->GetCovarianceBar().topLeftCorner(6, 6);
        MatXd P_post = filter->GetCovariancePost().topLeftCorner(6, 6);
        filter_cfg_.process_noise.Update(dx, Sigma_dx, P_bar, P_post);
      } else {
        // Reset Process Noise stocks
        filter_cfg_.process_noise.Reset();
      }
    }
  }

  ProcessNoiseFunction StateMeasManager::CreateProcessNoiseFunction(Ptr<KalmanFilter> filter,
                                                                    int N_ambiguity, bool verbose) {
    ProcessNoiseFunction proc_noise_func
        = [this, filter, N_ambiguity, verbose](const State& x, Real t0, Real tf) -> MatXd {
      int n_state = x.size() - N_ambiguity;
      int proc_noise_size = 0;
      int N_rva = filter_state_.N_rva;
      int N_clk = filter_state_.N_clk;
      int N_srp = filter_state_.N_srp;

      MatXd Q = MatXd::Zero(x.size(), x.size());
      MatXd Q_rva = MatXd::Zero(N_rva, N_rva);
      Real dt = tf - t0;

      // Position, Velocity, Acceleration
      if (filter_cfg_.use_adaptive_process_noise) {
        Q_rva = filter_cfg_.process_noise.ComputeProcessNoise();
        Q.topLeftCorner(N_rva, N_rva) = Q_rva;
        proc_noise_size += N_rva;
      } else {
        Q_rva = ProcessNoisePosVel(filter_cfg_.Q_a, dt);
        Q.topLeftCorner(N_rva, N_rva) = Q_rva;
        proc_noise_size += N_rva;
      }

      // Clock
      MatXd Q_clk(3, 3);
      if (N_clk == 3) {
        // Q_clk = filter_cfg_.clock_dynamics_->ThreeStateNoise(filter_cfg_.clk_model_sat, tf - t0);
        Q_clk = C * C
                * ThreeStateClockNoise(filter_cfg_.clk_model_sat, tf - t0, filter_cfg_.inflate_q2);
      } else if (N_clk == 2) {
        Q_clk.resize(2, 2);
        // Q_clk = filter_cfg_.clock_dynamics_->TwoStateNoise(filter_cfg_.clk_model_sat, tf - t0);
        Q_clk = C * C
                * TwoStateClockNoise(filter_cfg_.clk_model_sat, tf - t0, filter_cfg_.inflate_q2);
      } else {
        // Error
        throw std::runtime_error("Unsupported number of clock states: " + std::to_string(N_clk));
      }
      Q.block(N_rva, N_rva, N_clk, N_clk) = Q_clk;
      proc_noise_size += N_clk;

      // SRP
      if (N_srp > 0) {
        double Q_srp = 1e-14;  // Very small process noise for SRP
        Q(N_rva + N_clk, N_rva + N_clk) = Q_srp;
        proc_noise_size += N_srp;
      }

      // Integer Ambiguity
      if (N_ambiguity > 0) {
        double tau_ambiguity = filter_cfg_.tau_ambiguity;
        double q_ambiguity = filter_cfg_.q_ambiguity;
        MatXd Q_ambiguity = MatXd::Zero(N_ambiguity, N_ambiguity);
        double exp_factor = exp(-2 * dt / tau_ambiguity);
        Q_ambiguity.diagonal().array() = q_ambiguity * tau_ambiguity / 2 * (1 - exp_factor);
        Q.block(N_rva + N_clk + N_srp, N_rva + N_clk + N_srp, N_ambiguity, N_ambiguity)
            = Q_ambiguity;
        proc_noise_size += N_ambiguity;
      }

      if (filter_cfg_.use_udu) {
        int start_idx = 0;
        MatXd Q_new = MatXd::Zero(proc_noise_size, proc_noise_size);
        MatXd G = MatXd::Zero(x.size(), proc_noise_size);
        // For Position, Velocity, Acc
        VecMatPair DU_rva = UDUDecomposition(Q_rva);
        Q_new.topLeftCorner(N_rva, N_rva) = DU_rva.first.asDiagonal();
        G.topLeftCorner(N_rva, N_rva) = DU_rva.second;
        // For Clock
        start_idx += N_rva;
        VecMatPair DU_clk = UDUDecomposition(Q_clk);
        Q_new.block(start_idx, start_idx, N_clk, N_clk) = DU_clk.first.asDiagonal();
        G.block(start_idx, start_idx, N_clk, N_clk) = DU_clk.second;
        // For SRP
        start_idx += N_clk;
        if (N_srp > 0) {
          Q_new.block(start_idx, start_idx, N_srp, N_srp)
              = Q.block(N_rva + N_clk, N_rva + N_clk, N_srp, N_srp);
          G.block(start_idx, start_idx, N_srp, N_srp) = MatXd::Identity(N_srp, N_srp);
        }
        // For Integer Ambiguity
        start_idx += N_srp;
        if (N_ambiguity > 0) {
          Q_new.block(start_idx, start_idx, N_ambiguity, N_ambiguity)
              = Q.block(N_rva + N_clk + N_srp, N_rva + N_clk + N_srp, N_ambiguity, N_ambiguity);
          G.block(start_idx, start_idx, N_ambiguity, N_ambiguity)
              = MatXd::Identity(N_ambiguity, N_ambiguity);
        }
        filter->SetProcessNoiseMappingMatrix(G);
        Q = Q_new;
      }

      if (verbose) {
        std::cout << "[Process Noise Function]" << std::endl;
        if (filter_cfg_.use_udu) {
          std::cout << "Process Noise Mapping Matrix G (size: "
                    << filter->GetProcessNoiseMappingMatrix().rows() << "x"
                    << filter->GetProcessNoiseMappingMatrix().cols() << ")\n"
                    << filter->GetProcessNoiseMappingMatrix() << std::endl;
        }
        if (N_ambiguity > 0) {
          std::cout << "q_ambiguity: " << filter_cfg_.q_ambiguity << std::endl;
          std::cout << "tau_ambiguity: " << filter_cfg_.tau_ambiguity << std::endl;
          std::cout << "exp_factor: " << exp(-2 * (tf - t0) / filter_cfg_.tau_ambiguity)
                    << std::endl;
        }
        std::cout << "q2 inflate factor: " << filter_cfg_.inflate_q2 << std::endl;
        std::cout << std::scientific << std::setprecision(3);
        std::cout << "Process Noise (size: " << Q.rows() << "x" << Q.cols() << ")" << std::endl;

        auto old_flags = std::cout.flags();
        auto old_precision = std::cout.precision();
        std::cout << "rva: \n" << Q.topLeftCorner(N_rva, N_rva) << std::endl;
        std::cout << "clk: \n" << Q.block(N_rva, N_rva, N_clk, N_clk) << std::endl;
        if (Q.rows() > N_rva + N_clk) {
          std::cout << "srp: \n"
                    << Q.block(N_rva + N_clk, N_rva + N_clk, N_srp, N_srp) << std::endl;
        }
        if (N_ambiguity > 0) {
          std::cout << "ambiguity: \n"
                    << Q.block(N_rva + N_clk + N_srp, N_rva + N_clk + N_srp, N_ambiguity,
                               N_ambiguity)
                    << std::endl;
        }
        std::cout << " " << std::endl;
        std::cout.flags(old_flags);
        std::cout.precision(old_precision);
      }

      return Q;
    };

    return proc_noise_func;
  }

}  // namespace filtering_sim

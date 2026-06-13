#pragma once

#include <lupnt/lupnt.h>

#include "filter_setup.h"
#include "gnssmeas_loader.h"
#include "measurement_storage.h"
#include "simulation_config.h"
#include "udu_filter.h"

namespace filtering_sim {
  using namespace lupnt;

  using StateCovPair = std::pair<VecXd, MatXd>;

  constexpr double AMBIGUITY_OFFSET_M = 384400000.0;  // Mean offset for ambiguity states
  constexpr double AMBIGUITY_OFFSET_N
      = AMBIGUITY_OFFSET_M / 0.19;  // Mean offset for ambiguity states in N units

  enum OUTLIER_TYPE {
    NO_OUTLIER = 0,
    LARGE_RESIDUAL = 1,
    UNRELIABLE_STATE = 2,
    CYCLE_SLIP_DISCARDED = 3,
    CYCLE_SLIP_CORRECTED = 4
  };

  struct MeasData {
    int tidx;                                   // Time index
    double tspan;                               // Measurement time
    double t_tai;                               // TAI time of measurement
    int n_meas;                                 // Number of measurements
    VecXd z;                                    // Measurement vector
    VecXd R_diag;                               // Measurement noise diagonal
    VecXd N_integer_est;                        // Integer ambiguity estimates
    VecXd ambig_range_est;                      // Estimated integer ranges
    VecXd true_delays;                          // True delays for measurements
    VecXd est_min_alt;                          // Estimated minimum altitudes for measurements
    std::vector<bool> is_outlier;               // Cycle slip flags for each measurement (n_meas, )
    std::vector<double> cycle_slip_correction;  // Cycle slip magnitudes for each measurement

    std::vector<int>
        ambig_est_state_index_for_allmeas;      // Indices of ambiguity states in the filter state
    std::vector<int> ambig_full_state_indices;  // Indices of ambiguity states in the filter
                                                // state (flattened)
    std::vector<int> prns;                      // PRNs corresponding to measurements
    std::vector<MeasFreq> freqs;                // Measurement frequencies
    std::vector<MeasType> types;                // Measurement types
    std::vector<FreqPrnPair> freq_prn_pairs;    // Frequency-PRN pairs
    std::vector<Vec3d> r_ephem_mci;             // Ephemeris positions in MCI frame
    std::vector<double> clock_ephem;            // GNSS clock ephemeris values
    std::vector<Vec3d>
        r_prev_ephem_mci;  // Previous ephemeris positions in MCI frame (for TDCP, different size)
    std::vector<double>
        clock_prev_ephem;  // Previous GNSS clock ephemeris values (for TDCP, different size)
    std::vector<int> true_cycle_slips;        // True cycle slip indicators for each measurement
    std::vector<OUTLIER_TYPE> outlier_types;  // Outlier types for each measurement
    std::vector<double> outlier_thresholds;   // Outlier detection thresholds for each measurement
    VecXd S_diag;                             // Innovation covariance diagonal
    double tspan_prev;
    double t_tai_prev;
  };

  struct VecMeasData {
    std::vector<int> db_indices;  // (n_signal, ) Indices of measurements in the database
    std::vector<std::pair<int, int>> freq_db_pairs = {};  // (n_pair, ) Frequency-PRN pairs
    std::vector<double> z_meas_vec;                       // (n_meas, ) Measurement vector
    std::vector<double> R_diag_vec;                       // (n_meas, ) Measurement noise diagonal
    std::vector<Vec3d> r_ephem_mci_vec;       // (n_meas, ) Ephemeris positions in MCI frame
    std::vector<double> clk_ephem_vec;        // (n_meas, ) GNSS clock ephemeris values
    std::vector<Vec3d> r_prev_ephem_mci_vec;  // (n_meas, ) Previous ephemeris positions in MCI
                                              // frame (for TDCP, different size)
    std::vector<double> clk_prev_ephem_vec;  // (n_meas, ) Previous GNSS clock ephemeris values (for
                                             // TDCP, different size)
    std::vector<MeasType> meas_types_vec;    // (n_meas, ) Measurement types
    std::vector<MeasFreq> meas_freqs_vec;    // (n_meas, ) Measurement frequencies
    std::vector<int> prns_vec;               // (n_meas, ) PRNs corresponding to measurements
    std::vector<double> true_delays_vec;     // (n_meas, ) True delays for each measurement
    std::vector<int>
        meas_to_signal_index_vec;  // (n_meas, ) Mapping from measurement index to signal index
    std::vector<bool>
        use_signal_carrier_vec;  // (n_signal, )  Whether to use carrier phase for each signal
    std::vector<bool> use_signal_carrier_ionofree_vec;  // (n_pair, )  Whether to use iono-free
                                                        // carrier phase for each signal
    std::vector<double> iono_free_pr_minus_cp;        // (n_ionofree, ) Iono-free pseudorange minus
                                                      // carrier phase for each frequency-PRN pair
    std::vector<double> iono_free_pr_minus_cp_sigma;  // (n_ionofree, ) Iono-free carrier phase for
                                                      // each frequency-PRN pair

    std::vector<int> ambig_full_state_indices;  // (n_ambig, ) Indices of ambiguity states in the
                                                // full filter state
    std::vector<int>
        ambig_est_state_index_for_allmeas_vec;  // (n_meas, ) Indices of ambiguity states in the est
                                                // filter state for each measurement
    std::vector<int> true_cycle_slip_vec;       // Indices of measurements with true cycle slips

    std::vector<double> min_alt_vec;                 // Placeholder for min altitudes if needed
    std::map<int, double> idx_to_carrier_phase_map;  // Map to store carrier ranges for TDCP
  };

  class StateMeasManager {
  private:
    FilterConfig filter_cfg_;
    FilterState filter_state_;
    SimulationConfig sim_config_;
    MeasStorage meas_storage_;
    JointState joint_state_;
    bool first_measurement_ = true;

    // Time span of measurements
    VecXd tspan_;
    VecXd t_tais_;

    // Temporary variables
    Vec3 rx_pos_prev_;
    Real rx_clk_prev_;

    std::vector<std::filesystem::path> filtering_results_path_;
    std::vector<std::filesystem::path> smoothing_results_path_;

    bool run_filter_ = true;

  public:
    StateMeasManager(const SimulationConfig& sim_config, const GnssMeasLoader& gnss_meas_loader,
                     int seed, const std::filesystem::path& output_dir, bool verbose = false);
    ~StateMeasManager();

    /********************************************************
     * Setup functions
     *******************************************************/

    /**
     * Setup GNSS measurements from the measurement loader
     * @param gnss_meas_loader GNSS measurement loader
     * @param seed Random seed for any stochastic processes
     */
    void SetupGnssMeasurement(const GnssMeasLoader& gnss_meas_loader, int seed);

    /**
     * Initialize ambiguity states in the filter state
     * @param verbose Flag to enable verbose output
     * @return Updated FilterState with initialized ambiguity states
     */
    void InitializeAmbiguityStates(bool verbose);

    /********************************************************
     * Getters for internal states and configurations
     *******************************************************/
    const SimulationConfig& GetSimulationConfig() const { return sim_config_; }
    const MeasStorage& GetMeasStorage() const { return meas_storage_; }
    const FilterState& GetFilterState() const { return filter_state_; }
    FilterState& GetNonconstFilterState() { return filter_state_; }
    void SetFilterState(const FilterState& filter_state) { filter_state_ = filter_state; }
    const FilterConfig& GetFilterConfig() const { return filter_cfg_; }
    const JointState& GetJointState() const { return joint_state_; }
    JointState& GetJointState() { return joint_state_; }
    const std::filesystem::path& GetFilteringResultsPath(int iter) const {
      return filtering_results_path_[iter];
    }
    const std::filesystem::path& GetSmoothingResultsPath(int iter) const {
      return smoothing_results_path_[iter];
    }
    bool IsRunFilter() const { return run_filter_; }

    /********************************************************
     * PRN Utilities
     *******************************************************/
    int GetUniquePrnCount() const { return meas_storage_.GetUniquePrnCount(); }

    std::vector<int> GetAllUniquePrns() const;

    std::vector<int> GetTrackedPrnsAtTimeIndex(int tidx) const;

    int GetTrueAmbiguity(int tidx, MeasFreq freq, int prn) const;

    /********************************************************
     * State Conversion Utilities
     *******************************************************/
    StateCovPair FullStateToEstState(const VecXd& x_full, const MatXd& P_full,
                                     std::vector<int>& ambig_state_indices, bool is_tdcp);
    StateCovPair EstStateToFullState(const VecXd& x_est, const MatXd& P_est,
                                     std::vector<int>& ambig_state_indices, bool is_tdcp);
    int FullIndexToEstIndex(int full_index, const std::vector<int>& ambig_state_indices) const;
    int EstIndexToFullIndex(int est_index, const std::vector<int>& ambig_state_indices) const;

    /********************************************************
     * Measurement Manipulation at Each Time Step
     *******************************************************/
    /**
     * Get measurements based on the current state estimate
     * @param x_est Current state estimate vector
     * @param P_est Current state covariance matrix
     * @param tidx Current time index
     * @return MeasData structure containing measurement information
     */

    MeasData CreateEmptyMeasData(int tidx);
    void AddSingleFreqMeasurement(int tidx, State& x_est_full, VecMeasData& vmd);
    void AddMultiFreqMeasurement(int tidx, State& x_est_full, VecMeasData& vmd);
    VecXd UpdateIntegerAmbiguities(int tidx, State& x_est_full, MatXd& P_est_full,
                                   VecMeasData& vmd);
    bool IsUseMeasurement(MeasFreq freq, MeasType type, double min_alt, int tidx) const;

    /**
     * Get measurements based on the current state estimate
     * @param x_est_full Current full state estimate vector
     * @param P_est_full Current full state covariance matrix
     * @param tidx Current time index
     * @return MeasData structure containing measurement information
     */
    MeasData GetMeasurement(State& x_est_full, MatXd& P_est_full, int tidx);

    VecX ComputePredictedMeasurement(const VecX& x_state, const MeasData& meas_data,
                                     MatXd* H = nullptr, int meas_idx = -1);

    double ComputeMinEarthAltitude(double t_tai, const Vec3& rx_pos_mci, const Vec3& tx_pos_mci);

    void DetectCycleSlip(MeasData& meas_data, const VecX& dz, const MatXd* S, int meas_idx = -1);
    void RemoveCorrectCycleSlip(MeasData& meas_data, KalmanFilter* filter, int meas_idx = -1);
    FaultDetectionFunction CreateFaultDetectionFunction(MeasData& meas_data, int meas_idx,
                                                        bool verbose);

    /**
     * Update the filter measurement function and state based on the current time index
     * @param filter Pointer to the filter to be updated
     * @param x_est Current state estimate vector
     * @param P_est Current state covariance matrix
     * @param tidx Current time index
     */
    FilterMeasurementFunction CreateMeasurementFunction(MeasData& meas_data, bool verbose);
    /**
     * Create a measurement function for a single measurement (for UDU filters that process
     * measurement one by one)
     * @param filter Pointer to the filter to be updated
     * @param meas_data Measurement data structure
     * @param meas_idx Index of the measurement to process
     * @param verbose Whether to enable verbose logging
     * @return FilterMeasurementFunction object for the specified measurement
     */
    FilterMeasurementFunction CreateSingleMeasurementFunction(MeasData& meas_data, int meas_idx,
                                                              bool verbose);

    /**
     * Create a dynamics function for the filter
     * @param N_ambiguity Number of ambiguity states
     * @param verbose Whether to enable verbose logging
     * @return FilterDynamicsFunction object
     */
    FilterDynamicsFunction CreateDynamicsFunction(int N_ambiguity, bool verbose);

    /**
     * Create a process noise function for the filter
     * @param N_ambiguity Number of ambiguity states
     * @param verbose Whether to enable verbose logging
     * @return ProcessNoiseFunction object
     */
    void UpdateProcessNoise(Ptr<KalmanFilter> filter, const VecXd& z);
    ProcessNoiseFunction CreateProcessNoiseFunction(Ptr<KalmanFilter> filter, int N_ambiguity,
                                                    bool verbose);

    double GenerateNewAmbiguityEstimate(int meas_idx);
    double GenerateNewAmbiguityCovariance(int meas_index);

    void PrintMeasData(const MeasData& meas_data, const VecXd& z_pred) const;
    void PrintSingleMeasData(const MeasData& meas_data, const VecXd& z_pred, int meas_idx) const;

    int GetAmbiguityNum() const { return filter_state_.N_ambiguity_est; }
    int GetFullAmbiguityNum() const { return filter_state_.N_ambiguity_full; }
  };
}  // namespace filtering_sim

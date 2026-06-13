#include "src/filter_execution.h"

#include "lupnt/core/definitions.h"
#include "src/augmented_state_ekf.h"
#include "src/data_logger.h"
#include "src/filter_setup.h"
#include "src/udu_filter.h"

namespace filtering_sim {
  using namespace lupnt;

  Ptr<KalmanFilter> SetupFilter(StateMeasManager& smm, State& x_init, MatXd& P_init,
                                bool& is_tdcp) {
    Ptr<KalmanFilter> filter;
    FilterConfig filter_cfg = smm.GetFilterConfig();
    JointState joint_state = smm.GetJointState();

    is_tdcp = false;

    if (filter_cfg.use_tdcp || filter_cfg.use_ionofree_tdcp) {
      std::cout << "Using Delayed EKF for TDCP or ionofree TDCP measurements." << std::endl;
      if (filter_cfg.use_udu) {
        filter = MakePtr<UDUDelayedEKF>();
      } else {
        filter = MakePtr<DelayedEKF>();
      }
      is_tdcp = true;
    } else {
      if (filter_cfg.use_udu) {
        std::cout << "Using EKF with UDU decomposition for covariance." << std::endl;
        filter = MakePtr<UDUEKF>();
      } else {
        filter = MakePtr<EKF>();
      }
    }
    filter->SetDynamicsFunction(joint_state.GetDynamicsFunction());
    std::cout << "Set process noise function." << std::endl;
    filter->SetProcessNoiseFunction(joint_state.GetProcessNoiseFunction());
    std::cout << "Set initial state and covariance." << std::endl;

    InitializeFilterStates(filter, x_init, P_init, is_tdcp, filter_cfg.use_udu);

    std::cout << "Initialized EKF with state dimension: " << filter->GetState().size() << std::endl;

    return filter;
  }

  void InitializeFilterStates(Ptr<KalmanFilter> filter, State& x_init, MatXd& P_init, bool is_tdcp,
                              bool is_udu) {
    if (is_tdcp) {
      int n_state = x_init.size();
      // Expand to delayed state
      State delayed_state(n_state * 2);
      delayed_state.head(n_state) = x_init;
      delayed_state.tail(n_state) = x_init;
      filter->SetState(delayed_state);

      MatXd aug_cov = MatXd::Zero(2 * n_state, 2 * n_state);
      aug_cov.topLeftCorner(n_state, n_state) = P_init;
      aug_cov.topRightCorner(n_state, n_state) = P_init;
      aug_cov.bottomLeftCorner(n_state, n_state) = P_init;
      aug_cov.bottomRightCorner(n_state, n_state) = P_init;
      filter->SetCovariance(aug_cov);
    } else {
      filter->SetState(x_init);
      filter->SetCovariance(P_init);
    }
  }

  void RunFilterTimeStep(Ptr<KalmanFilter> filter, StateMeasManager& smm, EstData& est_data,
                         std::vector<MeasData>& meas_data_log, int tidx, bool is_tdcp,
                         bool use_meas_log, double& dz_sum, bool verbose) {
    // Variable
    StateCovPair pair_full_state_cov;
    FilterConfig filter_cfg = smm.GetFilterConfig();

    bool est_integers = filter_cfg.est_integers;

    // Store time index
    est_data.tidx = tidx;

    // Turn on/off verbose for each step
    std::tuple<bool, bool, bool, bool> verbose_settings = GetVerboseSettings(tidx, filter_cfg);
    bool verbose_dynamics = std::get<0>(verbose_settings);
    bool verbose_predict = std::get<1>(verbose_settings);
    bool verbose_measurement = std::get<2>(verbose_settings);
    bool verbose_update = std::get<3>(verbose_settings);

    // 1. Update Dynamics Function
    // ----------------------------------------------------------------
    FilterDynamicsFunction dyn_func_ambig
        = smm.CreateDynamicsFunction(est_data.ambig_full_state_indices.size(), verbose_dynamics);
    filter->SetDynamicsFunction(dyn_func_ambig);

    // 2. Update process noise function ------------------------------------------------
    ProcessNoiseFunction proc_noise_func = smm.CreateProcessNoiseFunction(
        filter, est_data.ambig_full_state_indices.size(), verbose_dynamics);
    filter->SetProcessNoiseFunction(proc_noise_func);

    // 3. Prediction step ------------------------------------------------
    double t_start = filter->GetTime().val();
    double t_end = smm.GetFilterState().ts_filter(tidx);
    filter->Predict(t_end);
    est_data.x_est = filter->GetState();
    est_data.P_est = filter->GetCovariance();
    PrintPredictionStep(tidx, est_data.x_est, est_data.P_est, smm, verbose_predict);

    // 4. Construct Full State for Measurement Generation ------------------------
    if (est_integers) {
      pair_full_state_cov = smm.EstStateToFullState(
          est_data.x_est, est_data.P_est, est_data.ambig_full_state_indices, est_data.is_tdcp);
      est_data.x_est_full = pair_full_state_cov.first;
      est_data.P_est_full = pair_full_state_cov.second;
    } else {
      est_data.x_est_full = est_data.x_est;
      est_data.P_est_full = est_data.P_est;
    }

    // 5. Generate measurements for this time step
    // ------------------------------------------------
    MeasData meas_data;
    // Load measureements from storage for second iteration onwards
    if (use_meas_log) {
      meas_data = meas_data_log[tidx];
      // Check that tidx matches tidx
      if (meas_data.tidx != tidx) {
        throw std::runtime_error("Measurement data time index does not match current time index.");
      }
    } else {
      // If not available, generate from state
      meas_data = smm.GetMeasurement(est_data.x_est_full, est_data.P_est_full, tidx);
      est_data.ambig_full_state_indices = meas_data.ambig_full_state_indices;
      meas_data_log[tidx] = meas_data;
    }

    if (est_integers) {
      pair_full_state_cov
          = smm.FullStateToEstState(est_data.x_est_full, est_data.P_est_full,
                                    est_data.ambig_full_state_indices, est_data.is_tdcp);
      est_data.x_est = pair_full_state_cov.first;
      est_data.P_est = pair_full_state_cov.second;

      // Update filter state and covariance with estimated integers
      // We need to do this since number of integers will change over time
      filter->SetState(est_data.x_est);
      filter->SetCovariance(est_data.P_est);
    }

    // 6. Measurement Update ----------------------------------------------------------------
    if (filter_cfg.use_udu && filter_cfg.udu_predict_single) {
      std::vector<double> dz_vector;
      for (int i = 0; i < meas_data.z.size(); i++) {
        auto meas_func_single
            = smm.CreateSingleMeasurementFunction(meas_data, i, verbose_measurement);
        // Update each measurement individually using UDU decomposition
        filter->SetFaultDetectionFunction(
            smm.CreateFaultDetectionFunction(meas_data, i, verbose_update));
        filter->SetMeasurementFunction(meas_func_single);
        filter->Update(meas_data.z.segment(i, 1));
      }
    } else {
      FilterMeasurementFunction meas_func
          = smm.CreateMeasurementFunction(meas_data, verbose_measurement);
      filter->SetFaultDetectionFunction(
          smm.CreateFaultDetectionFunction(meas_data, -1, verbose_update));
      filter->SetMeasurementFunction(meas_func);
      filter->Update(meas_data.z);
    }

    est_data.x_est = filter->GetState();
    est_data.P_est = filter->GetCovariance();

    if (filter_cfg.use_adaptive_process_noise) {
      smm.UpdateProcessNoise(filter, meas_data.z);
    }

    // 7. Log Results ----------------------------------------------------------------
    if (tidx >= filter_cfg.smooth_start_tidx) {
      // Only sum measurement residuals during smoothing period to directly compare with smoother
      VecXd z_pred = smm.ComputePredictedMeasurement(filter->GetState(), meas_data, nullptr, -1);
      VecXd dz = meas_data.z - z_pred;
      dz_sum += ComputeMeasResidualSum(dz, meas_data.R_diag, tidx, verbose_update);
    }
    pair_full_state_cov = smm.EstStateToFullState(est_data.x_est, est_data.P_est,
                                                  est_data.ambig_full_state_indices, is_tdcp);
    est_data.x_est_full = pair_full_state_cov.first;
    est_data.P_est_full = pair_full_state_cov.second;
    LogFilteringResults(filter, smm.GetNonconstFilterState(), pair_full_state_cov,
                        est_data.ambig_full_state_indices, meas_data.z.size(), tidx,
                        verbose_update);
  }

  double ComputeMeasResidualSum(VecXd& dz, const VecXd& R_diag, int tidx, bool verbose) {
    double dz_sum = 0.0;
    for (int i = 0; i < dz.size(); i++) {
      if (isnan(dz(i)) || isnan(R_diag(i)) || R_diag(i) <= 0) {
        dz(i) = 0.0;
      }
      double ratio = (dz(i) * dz(i)) / R_diag(i);
      if (ratio > 25.0) {
        continue;  // Skip large residuals
      }
      dz_sum += ratio;
      if (verbose) {
        std::cout << "Measurement Residual at tidx " << tidx << ", meas " << i << " / " << dz.size()
                  << ": dz = " << dz(i) << ", R_diag = " << R_diag(i) << " ratio = " << ratio
                  << " dz_sum = " << dz_sum << std::endl;
      }
    }
    if (verbose) {
      std::cout << " " << std::endl;
    }
    return dz_sum;
  }

  void RunSmoothingIteration(Ptr<KalmanFilter> filter, StateMeasManager& smm,
                             std::vector<MeasData>& meas_data_log, double& dz_sum, int tidx,
                             bool is_tdcp) {
    FilterConfig filter_cfg = smm.GetFilterConfig();

    // Update Smoother
    filter->UpdateSmoother(tidx);

    // Compute Measurement Residuals for Convergence Check
    bool verbose_sm = GetSmootherVerboseSetting(tidx, filter_cfg);
    VecXd z_true = meas_data_log[tidx + 1].z;
    if (z_true.size() > 0) {
      VecXd R_diag = meas_data_log[tidx + 1].R_diag;
      MatXd R_mat = R_diag.asDiagonal();
      MatXd H;
      VecXd x_sm = filter->GetSmoothedState(tidx + 1);       // N
      MatXd P_sm = filter->GetSmoothedCovariance(tidx + 1);  // N x N
      VecXd x_pred;
      if (is_tdcp) {
        // For delayed state, extract current state part
        x_pred.resize(x_sm.size() * 2);
        x_pred.head(x_sm.size()) = x_sm;
        x_pred.tail(x_sm.size()) = filter->GetSmoothedState(tidx);
      } else {
        x_pred = x_sm;
      }
      VecXd z_pred = smm.ComputePredictedMeasurement(x_pred, meas_data_log[tidx + 1], &H, -1);
      VecXd dz = z_true - z_pred;

      // Accumulate residuals
      dz_sum += ComputeMeasResidualSum(dz, R_diag, tidx, verbose_sm);
    }
    if (verbose_sm) {
      std::cout << " " << std::endl;
    }

    // Log Smoother Results
    LogSmoothingResults(filter, smm.GetNonconstFilterState(), dz_sum, tidx, verbose_sm);
  }

  void PrintMeasResidualSum(VecXd& dz_sum_log, int iter) {
    std::cout << "[measurement residual sum] " << std::endl;
    auto old_flags = std::cout.flags();
    auto old_precision = std::cout.precision();
    std::cout << std::scientific << std::setprecision(6);
    for (int i = 0; i <= iter; i++) {
      double dz_ratio = 0.0;
      if (i > 0) {
        dz_ratio = std::abs(dz_sum_log(i) - dz_sum_log(i - 1)) / dz_sum_log(i - 1);
      }
      std::cout << "  Iteration " << i + 1 << "|  |dz|/|R| = " << dz_sum_log(i)
                << "  |dz_(k-1) - dz(k)|/|dz_(k-1)| = " << dz_ratio << std::endl;
    }
    std::cout << " " << std::endl;
    std::cout.flags(old_flags);
    std::cout.precision(old_precision);
  }

  void ExecuteFilter(StateMeasManager& smm, int mci, bool verbose) {
    // Placeholder for filter execution logic
    // Implement the filter execution steps here, such as prediction, update, and logging

    if (!smm.IsRunFilter()) {
      std::cout << "Filter execution is disabled. Exiting ExecuteFilter." << std::endl;
      return;
    }

    // Print Settings
    auto old_flags = std::cout.flags();
    auto old_precision = std::cout.precision();

    // Setup EKF and Smoother ----------------------------------------------------------------
    FilterConfig filter_cfg = smm.GetFilterConfig();
    bool is_tdcp = false;

    State x_init = smm.GetJointState().GetStmState();
    MatXd P_init = smm.GetFilterState().P_est;
    Ptr<KalmanFilter> filter = SetupFilter(smm, x_init, P_init, is_tdcp);

    int iter = 0;  // For now, only single iteration
    int filter_start_tidx = 1;
    int smooth_start_tidx = filter_cfg.smooth_start_tidx;

    // Run Iteration -------------------------------------------------------------------------
    std::vector<MeasData> meas_data_log(filter_cfg.N_t_filter);
    std::vector<MatXd> P_est_log(filter_cfg.N_t_filter);
    std::vector<std::vector<int>> ambig_indices_log(filter_cfg.N_t_filter);

    VecXd dz_sum_filter_log(filter_cfg.N_smooth_iter);
    VecXd dz_sum_smooth_log(filter_cfg.N_smooth_iter);

    meas_data_log[0] = smm.CreateEmptyMeasData(0);

    for (iter = 0; iter < filter_cfg.N_smooth_iter; iter++) {
      std::cout << " " << std::endl;
      std::cout << "---------------------------------------------------------------------------"
                << std::endl;
      std::cout << "Starting filtering and smoothing iteration " << iter + 1 << " of "
                << filter_cfg.N_smooth_iter << std::endl;
      std::cout << "---------------------------------------------------------------------------"
                << std::endl;
      std::cout << " " << std::endl;

      EstData est_data;
      est_data.tidx = 0;
      est_data.x_est = filter->GetState();
      est_data.P_est = filter->GetCovariance();
      est_data.x_est_prior = filter->GetStatePrior();
      est_data.P_est_prior = filter->GetCovariancePrior();
      est_data.x_est_full = est_data.x_est;
      est_data.P_est_full = est_data.P_est;
      est_data.ambig_full_state_indices = {};  // Indices of ambiguity states in the filter state
      est_data.is_tdcp = is_tdcp;
      est_data.is_udu = (filter_cfg.use_udu);

      if (filter_cfg.run_smoother) {
        filter->InitializeLogger(filter_cfg.N_t_filter);
      }

      P_est_log[0] = est_data.P_est;
      ambig_indices_log[0] = est_data.ambig_full_state_indices;

      bool use_meas_log = false;
      if (iter > 0) {
        use_meas_log = true;  // Reuse stored measurement for the second iteration onwards
        filter_start_tidx = smooth_start_tidx + 1;
        est_data.ambig_full_state_indices = ambig_indices_log[smooth_start_tidx];
      }

      // Run filter -------------------------------------------------------------------------
      int total_filter_iters = filter_cfg.N_t_filter - filter_start_tidx;
      double dz_sum = 0.0;
      auto pbar_filter = Logger::GetProgressBar(
          total_filter_iters, "Running filter: iter " + std::to_string(iter), "Filter");
      for (int tidx = filter_start_tidx; tidx < filter_cfg.N_t_filter; tidx++) {
        RunFilterTimeStep(filter, smm, est_data, meas_data_log, tidx, is_tdcp, use_meas_log, dz_sum,
                          verbose);

        // Log results to smoother
        if (filter_cfg.run_smoother) filter->LogFilterEstimate(tidx);

        // Log Pest
        P_est_log[tidx] = est_data.P_est;
        ambig_indices_log[tidx] = est_data.ambig_full_state_indices;

        pbar_filter->Update();
      }
      pbar_filter->Finish();

      // Log filter residual sum
      dz_sum_filter_log(iter) = dz_sum;
      PrintMeasResidualSum(dz_sum_filter_log, iter);

      // Save Filtering Results
      SaveFilteringResults(smm.GetFilteringResultsPath(iter), smm.GetFilterState(),
                           dz_sum_filter_log);

      if (filter_cfg.run_smoother == false) {
        std::cout << "Smoother execution is disabled. Exiting smoothing step." << std::endl;
        return;  // Skip smoothing if not enabled
      }

      // Run Smoother ---------------------------------------------------------------------------
      auto pbar_smooth
          = Logger::GetProgressBar(filter_cfg.N_t_filter - smooth_start_tidx,
                                   "Running smoother: iter " + std::to_string(iter), "Smoother");

      filter->InitializeSmootherState();  // Initialize smoother state

      int counter = 0;
      dz_sum = 0.0;
      int num_meas = 0;
      for (int tidx = filter_cfg.N_t_filter - 2; tidx >= smooth_start_tidx; tidx--) {
        RunSmoothingIteration(filter, smm, meas_data_log, dz_sum, tidx, is_tdcp);
        pbar_smooth->Update();
      }
      pbar_smooth->Finish();

      // Log smoother residual sum
      dz_sum_smooth_log(iter) = dz_sum;
      PrintMeasResidualSum(dz_sum_smooth_log, iter);

      // Save Smoothing Results
      SaveSmoothingResults(smm.GetSmoothingResultsPath(iter), smm.GetFilterState(),
                           smooth_start_tidx, dz_sum_smooth_log);

      // Re-initialize filter for next iteration ------------------------------------------------
      // Set state and covariance from smoother results at the last time step
      State x_init = filter->GetSmoothedState(smooth_start_tidx);
      MatXd P_init = P_est_log[smooth_start_tidx];

      InitializeFilterStates(filter, x_init, P_init, is_tdcp, filter_cfg.use_udu);
    }
  }
}  // namespace filtering_sim

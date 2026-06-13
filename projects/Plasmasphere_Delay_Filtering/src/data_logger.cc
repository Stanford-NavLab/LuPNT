#include "src/data_logger.h"

#include <lupnt/lupnt.h>

#include <iomanip>

#include "src/filter_setup.h"
#include "src/state_meas_manager.h"

namespace filtering_sim {
  using namespace lupnt;

  MatXd RotCartToRtn(const Vec3& r, const Vec3& v) {
    Vec3 r_hat = r.normalized();
    Vec3 v_hat = v.normalized();
    Vec3 h_hat = r_hat.cross(v_hat).normalized();
    Vec3 theta_hat = h_hat.cross(r_hat).normalized();

    MatXd R(3, 3);
    R.col(0) = r_hat;
    R.col(1) = theta_hat;
    R.col(2) = h_hat;
    return R;
  }

  void LogFilteringResults(Ptr<KalmanFilter> filter, FilterState& state,
                           StateCovPair pair_full_state_cov,
                           std::vector<int> ambig_full_state_indices, int num_meas, int tidx,
                           bool verbose) {
    // Extract state
    State x_est = filter->GetState();
    MatXd P_est = filter->GetCovariance();
    VecXd x_est_full = pair_full_state_cov.first;
    MatXd P_est_full = pair_full_state_cov.second;

    VecXd ambig_est_full(state.N_ambiguity_full), ambig_sigma_full(state.N_ambiguity_full);
    if (state.N_ambiguity_est > 0) {
      ambig_est_full
          = x_est_full.segment(state.N_rva + state.N_clk + state.N_srp, state.N_ambiguity_full);
      ambig_sigma_full = P_est_full
                             .block(state.N_rva + state.N_clk + state.N_srp,
                                    state.N_rva + state.N_clk + state.N_srp, state.N_ambiguity_full,
                                    state.N_ambiguity_full)
                             .diagonal()
                             .cwiseSqrt();
    }

    // Log to FilterState
    state.rva_est_sat.row(tidx) = x_est.head(state.N_rva).transpose();
    state.rva_sigma_sat.row(tidx)
        = filter->GetCovariance().block(0, 0, state.N_rva, state.N_rva).diagonal().cwiseSqrt();
    state.clk_est_sat.row(tidx) = x_est.segment(state.N_rva, state.N_clk).transpose();
    state.clk_sigma_sat.row(tidx) = filter->GetCovariance()
                                        .block(state.N_rva, state.N_rva, state.N_clk, state.N_clk)
                                        .diagonal()
                                        .cwiseSqrt();
    state.srp_est_sat.row(tidx) = x_est.segment(state.N_rva + state.N_clk, state.N_srp).transpose();
    state.srp_sigma_sat.row(tidx)
        = filter->GetCovariance()
              .block(state.N_rva + state.N_clk, state.N_rva + state.N_clk, state.N_srp, state.N_srp)
              .diagonal()
              .cwiseSqrt();
    if (state.N_ambiguity_est > 0) {
      state.ambiguity_est_sat.row(tidx) = ambig_est_full.transpose();
      state.ambiguity_sigma_sat.row(tidx) = ambig_sigma_full.transpose();
    }

    // Measurement
    state.num_meas(tidx) = double(num_meas);

    // Convert to RTN frame
    MatXd R_rtn = RotCartToRtn(state.rva_true_sat.row(tidx).head<3>(),
                               state.rva_true_sat.row(tidx).segment<3>(3));
    Vec3 pos_err_xyz
        = state.rva_est_sat.row(tidx).head<3>() - state.rva_true_sat.row(tidx).head<3>();
    Vec3 vel_err_xyz
        = state.rva_est_sat.row(tidx).segment<3>(3) - state.rva_true_sat.row(tidx).segment<3>(3);
    Vec3 pos_err_rtn = R_rtn.transpose() * pos_err_xyz;
    Vec3 vel_err_rtn = R_rtn.transpose() * vel_err_xyz;
    Mat3 pos_err_cov_rtn = R_rtn.transpose() * filter->GetCovariance().block(0, 0, 3, 3) * R_rtn;
    Mat3 vel_err_cov_rtn = R_rtn.transpose() * filter->GetCovariance().block(3, 3, 3, 3) * R_rtn;
    Vec3 pos_sigma_rtn = pos_err_cov_rtn.diagonal().cwiseSqrt();
    Vec3 vel_sigma_rtn = vel_err_cov_rtn.diagonal().cwiseSqrt();

    state.rtn_err_pos.row(tidx) = pos_err_rtn.transpose();
    state.rtn_err_vel.row(tidx) = vel_err_rtn.transpose();
    state.rtn_err_pos_sigma.row(tidx) = pos_sigma_rtn.transpose();
    state.rtn_err_vel_sigma.row(tidx) = vel_sigma_rtn.transpose();

    if (verbose) {
      Eigen::IOFormat matFmt(Eigen::StreamPrecision, Eigen::Aligned,
                             "  ",  // coeff separator
                             "\n",  // row separator
                             "", "", "", "");

      // Print SettingsfE
      auto old_flags = std::cout.flags();
      auto old_precision = std::cout.precision();

      double pos_err = pos_err_xyz.norm();
      double vel_err = vel_err_xyz.norm();
      double clk_err
          = (state.clk_est_sat.row(tidx).col(0) - state.clk_true_sat.row(tidx).col(0)).norm();
      double srp_err = (state.srp_est_sat.row(tidx) - state.srp_true_sat.row(tidx)).norm();
      double pos_sigma = sqrt(state.rva_sigma_sat.row(tidx).head<3>().array().square().sum());
      double vel_sigma = sqrt(state.rva_sigma_sat.row(tidx).segment<3>(3).array().square().sum());
      double clk_sigma = sqrt(state.clk_sigma_sat.row(tidx).col(0).array().square().sum());
      double srp_sigma = sqrt(state.srp_sigma_sat.row(tidx).array().square().sum());

      double ambig_err = 0.0;
      double ambig_var = 0.0;
      double L1_lambda = get_wavelength(MeasFreq::L1);
      double ambig_sigma = 0.0;

      std::cout << "[UPDATE STEP]" << std::endl;
      std::cout << "Time index: " << tidx << std::endl;
      std::cout << "Number of Ambiguities Estimated: " << state.N_ambiguity_est << std::endl;
      std::cout << "Error after update (x_post): " << std::endl;
      std::cout << std::fixed << std::setprecision(3);
      std::cout << "  Pos: " << pos_err << " m  (3-sigma: " << 3 * pos_sigma << " m)" << std::endl;
      std::cout << "  Vel: " << vel_err * 1000 << " mm/s  (3-sigma: " << 3 * vel_sigma * 1000
                << " mm/s)" << std::endl;
      std::cout << "  Clk: " << clk_err << " m  (3-sigma: " << 3 * clk_sigma << " m)" << std::endl;
      std::cout << std::scientific << std::setprecision(3);
      std::cout << "  SRP: " << srp_err << " m^2/kg  (3-sigma: " << 3 * srp_sigma << " m^2/kg)"
                << std::endl;
      std::cout << "  Ambiguities: " << std::endl;
      if (state.N_ambiguity_est > 0) {
        for (int i = 0; i < state.N_ambiguity_est; i++) {
          int ambig_idx = ambig_full_state_indices[i] - (state.N_rva + state.N_clk + state.N_srp);
          double true_ambig = state.ambiguity_true_sat(tidx, ambig_idx);
          double est_ambig = state.ambiguity_est_sat(tidx, ambig_idx) + AMBIGUITY_OFFSET_N;
          double ambig_err_tmp = abs(est_ambig - true_ambig) * L1_lambda;
          ambig_err += ambig_err_tmp;
          double ambig_var_tmp = pow(state.ambiguity_sigma_sat(tidx, ambig_idx) * L1_lambda, 2);
          ambig_var += ambig_var_tmp;
          std::cout << "     Meas: " << i << " Full Idx: " << ambig_idx << " True: " << true_ambig
                    << " Est: " << est_ambig << " Error: " << ambig_err_tmp << " m"
                    << " (3-sigma: " << 3 * sqrt(ambig_var_tmp) << " m)" << std::endl;
        }
        ambig_err = ambig_err / state.N_ambiguity_est;
        ambig_var = ambig_var / state.N_ambiguity_est;
        ambig_sigma = sqrt(ambig_var);
      }
      VecXd dz = filter->GetMeasurementResidual();
      VecXd dx = filter->GetStateCorrection();
      MatXd S = filter->GetInnovationCov();

      if (state.N_ambiguity_est > 0) {
        std::cout << "  Ambiguities (Total): " << ambig_err << "  m (3-sigma: " << 3 * ambig_sigma
                  << " m)" << std::endl;
      }
      std::cout << "Number of measurements: " << num_meas << std::endl;
      std::cout << "Measurement residual (dz): " << dz.transpose() << std::endl;
      std::cout << "State correction (dx): " << dx.transpose() << std::endl;
      std::cout << "Innovation covariance (S): " << S.diagonal().transpose().array().cwiseSqrt()
                << std::endl;
      std::cout << std::scientific << std::setprecision(2);
      // std::cout << "Updated covariance (P_post) (size: " << P_est.rows() << "x" << P_est.cols()
      //           << ")" << std::endl
      //           << P_est.format(matFmt) << std::endl;
      std::cout << " " << std::endl;
      std::cout << " " << std::endl;

      // Restore
      std::cout.flags(old_flags);
      std::cout.precision(old_precision);
    }
  }  // namespace LogResults

  void PrintPredictionStep(int tidx, const VecXd& x_est, const MatXd& P_est,
                           const StateMeasManager& smm, bool verbose) {
    Eigen::IOFormat matFmt(Eigen::StreamPrecision, Eigen::Aligned,
                           "  ",  // coeff separator
                           "\n",  // row separator
                           "", "", "", "");

    if (verbose) {
      // Print Settings
      auto old_flags = std::cout.flags();
      auto old_precision = std::cout.precision();
      std::cout << "[PREDICTION STEP]" << std::endl;
      std::cout << "Time index: " << tidx << std::endl;
      std::cout << "Predicted state (x_prior): " << std::endl;
      std::cout << "  Pos: " << x_est.segment<3>(0).transpose() << std::endl;
      std::cout << "  Vel: " << x_est.segment<3>(3).transpose() << std::endl;
      std::cout << "  Clk: " << x_est.segment<2>(6).transpose() << std::endl;
      std::cout << "  SRP: "
                << x_est
                       .segment(smm.GetFilterState().N_rva + smm.GetFilterState().N_clk,
                                smm.GetFilterState().N_srp)
                       .transpose()
                << std::endl;
      std::cout << std::scientific << std::setprecision(2);
      // std::cout << "Predicted covariance (P_prior) (size: " << P_est.rows() << "x" <<
      // P_est.cols()
      //           << ")" << std::endl
      //           << P_est.format(matFmt) << std::endl;
      std::cout << " " << std::endl;

      // Restore
      std::cout.flags(old_flags);
      std::cout.precision(old_precision);
    }
  }

  std::tuple<bool, bool, bool, bool> GetVerboseSettings(int tidx, const FilterConfig& filter_cfg) {
    // Verbose settings interpretation:
    //   -1: disable
    //   0: first step only
    //   1: every step
    //   2-: every N steps
    auto check_verbose = [tidx](int setting) {
      if (setting < 0) {
        return false;
      } else if (setting == 0) {
        return tidx == 1;
      } else if (setting == 1) {
        return true;
      } else {
        return (tidx % setting) == 0;
      }
    };

    bool verbose_dynamics = check_verbose(filter_cfg.verbose_filter_dynamics);
    bool verbose_predict = check_verbose(filter_cfg.verbose_filter_predict);
    bool verbose_measurement = check_verbose(filter_cfg.verbose_filter_measurement);
    bool verbose_update = check_verbose(filter_cfg.verbose_filter_update);

    bool verbose = verbose_dynamics || verbose_predict || verbose_measurement || verbose_update;

    if (verbose) {
      std::cout << std::string(200, '=') << std::endl;
      std::cout << "Filter Time Step: " << tidx << " / " << filter_cfg.N_t_filter - 1 << std::endl;
    }

    return {verbose_dynamics, verbose_predict, verbose_measurement, verbose_update};
  }

  bool GetSmootherVerboseSetting(int tidx, const FilterConfig& filter_cfg) {
    // Verbose settings interpretation:
    //   -1: disable
    //   0: first step only
    //   1: every step
    //   2-: every N steps
    if (filter_cfg.verbose_filter_smoothing < 0) {
      return false;
    } else if (filter_cfg.verbose_filter_smoothing == 0) {
      return tidx == filter_cfg.N_t_filter - 2;
    } else if (filter_cfg.verbose_filter_smoothing == 1) {
      return true;
    } else {
      return (tidx % filter_cfg.verbose_filter_smoothing) == 0;
    }
  }

  void SaveFilteringResults(const std::filesystem::path filepath, const FilterState& state,
                            VecXd dz_sum_iter) {
    // Create results file
    H5Easy::File results_file = GetH5File(filepath, true);
    std::cout << "Executing filter and saving results to " << filepath << std::endl;

    // True States
    Dump(results_file, "/ts_filter", state.ts_filter);
    Dump(results_file, "/rva_true_sat", state.rva_true_sat);
    Dump(results_file, "/srp_true_sat", state.srp_true_sat);
    Dump(results_file, "/clk_true_sat", state.clk_true_sat);
    Dump(results_file, "/ambiguity_true_sat", state.ambiguity_true_sat);
    // Estimated States
    Dump(results_file, "/rva_est_sat", state.rva_est_sat);
    Dump(results_file, "/clk_est_sat", state.clk_est_sat);
    Dump(results_file, "/clk_true_sat", state.clk_true_sat);
    Dump(results_file, "/srp_est_sat", state.srp_est_sat);
    Dump(results_file, "/ambiguity_est_sat", state.ambiguity_est_sat);
    // Sigmas
    Dump(results_file, "/rva_sigma_sat", state.rva_sigma_sat);
    Dump(results_file, "/srp_sigma_sat", state.srp_sigma_sat);
    Dump(results_file, "/clk_sigma_sat", state.clk_sigma_sat);
    Dump(results_file, "/ambiguity_sigma_sat", state.ambiguity_sigma_sat);
    Dump(results_file, "/num_meas", state.num_meas);
    // RTN Errors
    Dump(results_file, "/pos_rtn_err", state.rtn_err_pos);
    Dump(results_file, "/vel_rtn_err", state.rtn_err_vel);
    Dump(results_file, "/pos_rtn_sigma", state.rtn_err_pos_sigma);
    Dump(results_file, "/vel_rtn_sigma", state.rtn_err_vel_sigma);
    // Measurement residual sum
    Dump(results_file, "/dz_sum_iter", dz_sum_iter);
    std::cout << "Saved filter results to " << filepath << std::endl;
  }

  void LogSmoothingResults(Ptr<KalmanFilter> filter, FilterState& state, double dz_sum, int tidx,
                           bool verbose) {
    // Extract smoothed state
    State x_smooth = filter->GetSmoothedState(tidx);
    MatXd P_smooth = filter->GetSmoothedCovariance(tidx);

    // Log to FilterState
    state.rva_est_sat.row(tidx) = x_smooth.head(state.N_rva).transpose();
    state.clk_est_sat.row(tidx) = x_smooth.segment(state.N_rva, state.N_clk).transpose();
    state.srp_est_sat.row(tidx)
        = x_smooth.segment(state.N_rva + state.N_clk, state.N_srp).transpose();

    // Log sigmas
    state.rva_sigma_sat.row(tidx)
        = P_smooth.block(0, 0, state.N_rva, state.N_rva).diagonal().cwiseSqrt();
    state.clk_sigma_sat.row(tidx)
        = P_smooth.block(state.N_rva, state.N_rva, state.N_clk, state.N_clk).diagonal().cwiseSqrt();
    state.srp_sigma_sat.row(tidx)
        = P_smooth
              .block(state.N_rva + state.N_clk, state.N_rva + state.N_clk, state.N_srp, state.N_srp)
              .diagonal()
              .cwiseSqrt();

    // Convert to RTN frame
    MatXd R_rtn = RotCartToRtn(state.rva_true_sat.row(tidx).head<3>(),
                               state.rva_true_sat.row(tidx).segment<3>(3));
    Vec3 pos_err_xyz
        = state.rva_est_sat.row(tidx).head<3>() - state.rva_true_sat.row(tidx).head<3>();
    Vec3 vel_err_xyz
        = state.rva_est_sat.row(tidx).segment<3>(3) - state.rva_true_sat.row(tidx).segment<3>(3);
    Vec3 pos_err_rtn = R_rtn.transpose() * pos_err_xyz;
    Vec3 vel_err_rtn = R_rtn.transpose() * vel_err_xyz;
    Mat3 pos_err_cov_rtn = R_rtn.transpose() * P_smooth.block(0, 0, 3, 3) * R_rtn;
    Mat3 vel_err_cov_rtn = R_rtn.transpose() * P_smooth.block(3, 3, 3, 3) * R_rtn;
    Vec3 pos_sigma_rtn = pos_err_cov_rtn.diagonal().cwiseSqrt();
    Vec3 vel_sigma_rtn = vel_err_cov_rtn.diagonal().cwiseSqrt();

    state.rtn_err_pos.row(tidx) = pos_err_rtn.transpose();
    state.rtn_err_vel.row(tidx) = vel_err_rtn.transpose();
    state.rtn_err_pos_sigma.row(tidx) = pos_sigma_rtn.transpose();
    state.rtn_err_vel_sigma.row(tidx) = vel_sigma_rtn.transpose();

    if (verbose) {
      double pos_err = pos_err_xyz.norm();
      double vel_err = vel_err_xyz.norm();
      double clk_err
          = (state.clk_est_sat.row(tidx).col(0) - state.clk_true_sat.row(tidx).col(0)).norm();
      double srp_err = (state.srp_est_sat.row(tidx) - state.srp_true_sat.row(tidx)).norm();
      double pos_sigma = sqrt(state.rva_sigma_sat.row(tidx).head<3>().array().square().sum());
      double vel_sigma = sqrt(state.rva_sigma_sat.row(tidx).segment<3>(3).array().square().sum());
      double clk_sigma = sqrt(state.clk_sigma_sat.row(tidx).col(0).array().square().sum());
      double srp_sigma = sqrt(state.srp_sigma_sat.row(tidx).array().square().sum());

      std::cout << "[SMOOTHING STEP]" << std::endl;
      std::cout << "Time index: " << tidx << std::endl;
      std::cout << "Number of Ambiguities Estimated: " << state.N_ambiguity_est << std::endl;
      std::cout << "Error after update (x_post): " << std::endl;
      std::cout << std::fixed << std::setprecision(3);
      std::cout << "  Pos: " << pos_err << " m  (3-sigma: " << 3 * pos_sigma << " m)" << std::endl;
      std::cout << "  Vel: " << vel_err * 1000 << " mm/s  (3-sigma: " << 3 * vel_sigma * 1000
                << " mm/s)" << std::endl;
      std::cout << "  Clk: " << clk_err << " m  (3-sigma: " << 3 * clk_sigma << " m)" << std::endl;
      std::cout << std::scientific << std::setprecision(3);
      std::cout << "  SRP: " << srp_err << " m^2/kg  (3-sigma: " << 3 * srp_sigma << " m^2/kg)"
                << std::endl;
      std::cout << "Residual RMS after smoothing: " << dz_sum << std::endl;
      std::cout << " " << std::endl;
    }

  }  // namespace LogSmoothingResults

  void SaveSmoothingResults(const std::filesystem::path filepath, const FilterState& state,
                            int smooth_start_tidx, VecXd dz_sum_iter) {
    // Create results file
    H5Easy::File results_file = GetH5File(filepath, true);
    std::cout << "Saving smoothing results to " << filepath << std::endl;

    int len_sm = state.N_t - smooth_start_tidx;

    // True States
    Dump(results_file, "/ts_filter", state.ts_filter.segment(smooth_start_tidx, len_sm));
    Dump(results_file, "/rva_true_sat",
         state.rva_true_sat.block(smooth_start_tidx, 0, len_sm, state.N_rva));
    Dump(results_file, "/clk_true_sat",
         state.clk_true_sat.block(smooth_start_tidx, 0, len_sm, state.N_clk));
    Dump(results_file, "/srp_true_sat", state.srp_true_sat.segment(smooth_start_tidx, len_sm));
    // Smoothed States
    Dump(results_file, "/ts_filter", state.ts_filter.segment(smooth_start_tidx, len_sm));
    Dump(results_file, "/rva_smooth_sat",
         state.rva_est_sat.block(smooth_start_tidx, 0, len_sm, state.N_rva));
    Dump(results_file, "/clk_smooth_sat",
         state.clk_est_sat.block(smooth_start_tidx, 0, len_sm, state.N_clk));
    Dump(results_file, "/srp_smooth_sat",
         state.srp_est_sat.block(smooth_start_tidx, 0, len_sm, state.N_srp));
    // Sigmas
    Dump(results_file, "/rva_sigma_sat",
         state.rva_sigma_sat.block(smooth_start_tidx, 0, len_sm, state.N_rva));
    Dump(results_file, "/clk_sigma_sat",
         state.clk_sigma_sat.block(smooth_start_tidx, 0, len_sm, state.N_clk));
    Dump(results_file, "/srp_sigma_sat",
         state.srp_sigma_sat.block(smooth_start_tidx, 0, len_sm, state.N_srp));
    Dump(results_file, "/clk_sigma_sat",
         state.clk_sigma_sat.block(smooth_start_tidx, 0, len_sm, state.N_clk));
    // RTN Errors
    Dump(results_file, "/pos_rtn_err", state.rtn_err_pos.block(smooth_start_tidx, 0, len_sm, 3));
    Dump(results_file, "/vel_rtn_err", state.rtn_err_vel.block(smooth_start_tidx, 0, len_sm, 3));
    Dump(results_file, "/pos_rtn_sigma",
         state.rtn_err_pos_sigma.block(smooth_start_tidx, 0, len_sm, 3));
    Dump(results_file, "/vel_rtn_sigma",
         state.rtn_err_vel_sigma.block(smooth_start_tidx, 0, len_sm, 3));
    Dump(results_file, "/dz_sum_iter", dz_sum_iter);
    std::cout << "Saved smoothing results to " << filepath << std::endl;
  }

}  // namespace filtering_sim

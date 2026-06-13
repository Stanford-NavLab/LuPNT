#pragma once

#include <lupnt/lupnt.h>

#include <iomanip>

#include "src/data_logger.h"
#include "src/state_meas_manager.h"

namespace filtering_sim {
  using namespace lupnt;

  MatXd RotCartToRtn(const Vec3& r, const Vec3& v);

  void LogFilteringResults(Ptr<KalmanFilter> filter, FilterState& state,
                           StateCovPair pair_full_state_cov,
                           std::vector<int> ambig_full_state_indices, int num_meas, int tidx,
                           bool verbose);

  void LogSmoothingResults(Ptr<KalmanFilter> filter, FilterState& state, double dz_sum, int tidx,
                           bool verbose);

  void PrintPredictionStep(int tidx, const VecXd& x_est, const MatXd& P_est,
                           const StateMeasManager& smm, bool verbose);

  std::tuple<bool, bool, bool, bool> GetVerboseSettings(int tidx, const FilterConfig& filter_cfg);

  bool GetSmootherVerboseSetting(int tidx, const FilterConfig& filter_cfg);

  void SaveFilteringResults(const std::filesystem::path filepath, const FilterState& state,
                            VecXd dz_sum_iter);

  void SaveSmoothingResults(const std::filesystem::path filepath, const FilterState& state,
                            int smooth_start_tidx, VecXd dz_sum_iter);

}  // namespace filtering_sim

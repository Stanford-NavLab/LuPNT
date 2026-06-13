#pragma once

#include <lupnt/lupnt.h>
#include <yaml-cpp/yaml.h>

#include <cassert>
#include <filesystem>
#include <set>
#include <vector>

#include "src/gnssmeas_loader.h"
#include "src/measurement_storage.h"
#include "src/simulation_config.h"
#include "src/state_meas_manager.h"

namespace filtering_sim {
  using namespace lupnt;
  template <typename T> using vec = std::vector<T>;

  struct EstData {
    int tidx = 0;
    State x_est;
    MatXd P_est;
    State x_est_full;
    MatXd P_est_full;
    // Used for Smoother logging
    State x_est_prior;
    MatXd P_est_prior;
    std::vector<int> ambig_full_state_indices;
    bool is_tdcp;
    bool is_udu;
    MeasData meas_data;
  };

  Ptr<KalmanFilter> SetupFilter(StateMeasManager& smm, State& x_init, MatXd& P_init, bool& is_tdcp);

  void InitializeFilterStates(Ptr<KalmanFilter> filter, State& x_init, MatXd& P_init, bool is_tdcp,
                              bool is_udu);

  void RunFilterTimeStep(Ptr<KalmanFilter> filter, StateMeasManager& smm, EstData& est_data,
                         std::vector<MeasData>& meas_data_log, int tidx, bool is_tdcp,
                         bool use_meas_log, double& dz_sum, bool verbose);

  void RunSmoothingIteration(Ptr<KalmanFilter> filter, StateMeasManager& smm,
                             std::vector<MeasData>& meas_data_log, double& dz_sum, int tidx,
                             bool is_tdcp);

  double ComputeMeasResidualSum(VecXd& dz, const VecXd& R_diag, int tidx, bool verbose);

  void PrintMeasResidualSum(VecXd& dz_sum_log, int iter);

  void ExecuteFilter(StateMeasManager& meas_manager, int mci, bool verbose = false);

}  // namespace filtering_sim

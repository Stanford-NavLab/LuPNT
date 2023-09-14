/**
 * @file StateEstimationApp.cpp
 * @author Stanford University NAV Lab
 * @brief Base class for state estimation applications
 * @version
 * @date
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "StateEstimationApp.h"

namespace LPT {

void StateEstimationApp::SetDynamicsFunction() {
  auto dynfunc = [&](ad::VectorXreal x, ad::real t_curr, ad::real t_end,
                     Eigen::MatrixXd &Phi) {
    std::vector<IState *> state_vec = state_vec_.GetJointState();

    // Iterate for each dynamics and corresponding state (e.g. orbit and
    // dynamics)
    int start_idx = 0;
    for (int i = 0; i < dynamics_vec_.size(); i++) {
      int state_size = state_vec[i]->GetStateSize();
      Eigen::MatrixXd Phi_tmp(state_size, state_size);
      ad::VectorXreal x_seg(state_size);
      for (int j = 0; j < state_size; j++) {
        x_seg(j) = x(start_idx + j);
      }
      dynamics_vec_[i]->PropagateWithStm(x_seg, t_curr, t_end, Phi_tmp);
      Phi.block(start_idx, start_idx, state_size, state_size) = Phi_tmp;
      for (int j = 0; j < state_size; j++) {
        x(start_idx + j) = x_seg(j);
      }
      // Add states
      start_idx += state_size;
    }

    return x;
  };

  dynamics_func_ = dynfunc;
};

}  // namespace LPT
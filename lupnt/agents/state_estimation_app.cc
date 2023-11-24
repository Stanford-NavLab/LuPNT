/**
 * @file StateEstimationApp.cpp
 * @author Stanford NAV LAB
 * @brief Base class for State Estimation Application
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#include "state_estimation_app.h"

namespace lupnt {

void StateEstimationApp::SetDynamicsFunction() {
  auto dynfunc = [&](VectorXreal x, real t_curr, real t_end, MatrixXd &Phi) {
    std::vector<IState *> state_vec = state_vec_.GetJointState();

    // Iterate for each dynamics and corresponding state (e.g. orbit and
    // dynamics)
    int start_idx = 0;
    for (int i = 0; i < dynamics_vec_.size(); i++) {
      int state_size = state_vec[i]->GetSize();
      MatrixXd Phi_tmp(state_size, state_size);
      VectorXreal x_seg(state_size);
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

}  // namespace lupnt
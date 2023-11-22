/**
 * @file filters.cpp
 * @author Stanford NAVLAB
 * @brief Implementation of Filters
 * @version 0.1
 * @date 2023-09-09
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "filters.h"

namespace lupnt {

// Extended Kalman Filter
void EKF::Predict(real t_end) {
  MatrixXd Phi(x.size(), x.size());
  MatrixXd Q(x.size(), x.size());
  x = this->dynamics(x, t_curr, t_end, Phi);
  Q = this->process_noise(x, t_curr, t_end);
  P = Phi * P * Phi.transpose() + Q;
  xbar = x;
  Pbar = P;
}

void EKF::Update(VectorXreal z_obs) {
  int n = x.size();
  int m = z_obs.size();
  MatrixXd H(m, n);
  MatrixXd R(m, m);
  VectorXreal z_predict = this->measurement(xbar, H, R);
  S = R + H * P * H.transpose();
  K = P * H.transpose() * S.inverse();
  dx = K * (z_obs - z_predict);
  x = x + dx;
  I = MatrixXd::Identity(x.size(), x.size());
  P = (I - K * H) * P * (I - K * H).transpose() +
      K * R * K.transpose();  // Joseph form
}

}  // namespace lupnt
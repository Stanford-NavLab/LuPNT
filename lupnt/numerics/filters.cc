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

/*****************************************************
 *   Extended Kalman Filter
 *****************************************************/

/**
 * @brief Predict step
 *
 * @param t_end
 */
void EKF::Predict(real t_end) {
  int n = x.size();
  MatrixXd Phi(n, n);
  MatrixXd Q(n, n);
  P.resize(n, n);
  x = this->dynamics(x, t_curr, t_end, Phi);
  Q = this->process_noise(x, t_curr, t_end);
  P = Phi * P * Phi.transpose() + Q;
  xbar = x;
  Pbar = P;
  t_curr = t_end;
}

/**
 * @brief Update step
 *
 * @param z_obs observed measurement
 */
void EKF::Update(VectorX z_obs) {
  int n = x.size();
  int m = z_obs.size();

  // allocate memory (without this, MatrixXd will cause segfault)
  MatrixXd H(m, n);
  MatrixXd R(m, m);
  S.resize(m, m);
  K.resize(n, m);
  dx.resize(n);

  VectorX z_predict = this->measurement(xbar, H, R);

  // Update step
  S = R + H * P * H.transpose();        // Measurement information
  K = P * H.transpose() * S.inverse();  // Kalman gain
  dx = K * (z_obs - z_predict);
  x = x + dx;
  I = MatrixXd::Identity(n, n);
  MatrixXd G = MatrixXd(n, n);
  G = I - K * H;
  P = G * P * G.transpose() + K * R * K.transpose();  // Joseph form
}

/*****************************************************
 *   EKF Smoother  (Todo)
 *****************************************************/

/*****************************************************
 *   Information Filter  (Todo)
 *****************************************************/

/*****************************************************
 *   Square Root Information Filter  (Todo)
 *****************************************************/

/*****************************************************
 *   Batch Filter  (Todo)
 *****************************************************/

/*****************************************************
 *   Unscented Kalman Filter  (Todo)
 *****************************************************/

/*****************************************************
 *   Particle Filter (Todo)
 *****************************************************/

}  // namespace lupnt
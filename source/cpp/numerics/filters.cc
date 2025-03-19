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

#include "lupnt/numerics/filters.h"

namespace lupnt {

  /*****************************************************
   *   Extended Kalman Filter
   *****************************************************/

  void EKF::Predict(Real t_end) {
    int n = x_.size();
    Phi_.resize(n, n);
    P_.resize(n, n);
    Q_.resize(n, n);

    Q_ = this->process_noise_(x_, t_curr_, t_end);
    x_ = this->dynamics_(x_, t_curr_, t_end, Phi_);

    P_ = Phi_ * P_ * Phi_.transpose() + Q_;
    xbar_ = x_;
    Pbar_ = P_;
    t_curr_ = t_end;
  }

  int EKF::RemoveOutliers(int m_orig, bool debug) {
    std::vector<int> is_outlier(m_orig);
    VecXd ratio(m_orig);
    int n_valid = 0;
    int n = x_.size();

    for (int i = 0; i < m_orig; i++) {
      ratio(i) = abs(dy_(i).val() / sqrt(S_(i, i)));
      if (ratio(i) > outlier_threshold_) {
        is_outlier[i] = 1;
      } else {
        is_outlier[i] = 0;
        n_valid++;
      }
    }

    // Remove outliers
    int m = n_valid;

    if (m == m_orig) {
      return m;  // all measurement valid, nothing to change
    } else {
      if (debug) {
        std::cout << "Removing " << m_orig - m << "/" << m_orig
                  << " outliers - ratio: " << ratio.transpose() << std::endl;
      }
    }

    MatXd H_new(m, n);
    MatXd R_new(m, m);
    MatXd S_new(m, m);
    VecX dy_new(m);
    int j = 0;
    int l = 0;

    for (int i = 0; i < m_orig; i++) {
      if (is_outlier[i] == 0) {
        H_new.row(j) = H_.row(i);
        dy_new(j) = dy_(i);

        l = 0;
        for (int k = 0; k < m_orig; k++) {
          if (is_outlier[k] == 0) {
            R_new(j, l) = R_(i, k);
            S_new(j, l) = S_(i, k);
            l++;
          }
        }

        j++;
      }
    }

    // set new value
    H_ = H_new;
    R_ = R_new;
    S_ = S_new;
    dy_ = dy_new;

    return m;
  }

  /**
   * @brief Update step
   *
   * @param z_true observed measurement
   */
  void EKF::Update(VecX z_true_in, bool debug) {
    z_true_ = z_true_in;

    int n = x_.size();
    int m = z_true_.size();

    if (m == 0) {
      return;  // no measurement, nothing to update
    }

    // allocate memory (without this, VecXd will cause segfault)
    H_.resize(m, n);
    S_.resize(m, m);
    dy_.resize(m);
    dx_.resize(n);
    R_.resize(m, m);

    z_pred_ = measurement_(x_, H_, R_);

    S_ = R_ + H_ * P_ * H_.transpose();  // Measurement information
    dy_ = z_true_ - z_pred_;

    // Remove outliers
    m = RemoveOutliers(m, debug);

    // re-allocate memory
    K_.resize(n, m);

    // Update step
    K_ = P_ * H_.transpose() * S_.inverse();  // Kalman gain
    dx_ = K_ * dy_;
    x_ = x_ + dx_;
    MatXd I = MatXd::Identity(n, n);
    MatXd G(n, n);
    G = I - K_ * H_;
    P_ = G * P_ * G.transpose() + K_ * R_ * K_.transpose();  // Joseph form
  }

  void EKF::Step(Real t_end, VecX z_obs, bool debug) {
    Predict(t_end);
    Update(z_obs, debug);
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

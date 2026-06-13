#include "src/augmented_state_ekf.h"

namespace filtering_sim {
  using namespace lupnt;

  void DelayedEKF::Predict(Real t, const State* u) {
    // Checks
    LUPNT_CHECK(x_.size() > 0, "(Predict) State vector is not set", "ASEKF");
    LUPNT_CHECK(P_.size() > 0, "(Predict) Covariance matrix is not set", "ASEKF");
    LUPNT_CHECK(f_dyn_, "(Predict) Dynamics function not set", "ASEKF");
    LUPNT_CHECK(f_proc_, "(Predict) Process noise function not set", "ASEKF");

    int N_state = int(x_.size() / 2);
    VecXd x_prev = x_.head(N_state);
    MatXd P_prev = P_.topLeftCorner(N_state, N_state);

    // Ensure the last half of x_ is the same
    x_.tail(N_state) = x_prev;
    P_.bottomRightCorner(N_state, N_state) = P_prev;
    P_.topRightCorner(N_state, N_state) = P_prev;
    P_.bottomLeftCorner(N_state, N_state) = P_prev;

    // Process noise
    Q_ = f_proc_(x_prev, t_, t);
    LUPNT_CHECK(!Q_.hasNaN(), "(Predict) Process noise covariance has NaN", "ASEKF");

    // Predict state
    x_.head(N_state) = f_dyn_(x_prev, t_, t, u, &F_);
    LUPNT_CHECK(!x_.hasNaN(), "(Predict) Predicted state has NaN", "ASEKF");
    // Update State
    MatXd P_cross = F_ * P_prev;
    MatXd P_next = P_cross * F_.transpose() + Q_;

    // Covariance Update
    P_.topLeftCorner(N_state, N_state) = P_next;
    P_.topRightCorner(N_state, N_state) = P_cross;
    P_.bottomLeftCorner(N_state, N_state) = P_cross.transpose();
    P_.bottomRightCorner(N_state, N_state) = P_prev;

    // Update time
    t_ = t;

    // Store
    x_prior_ = x_;
    P_prior_ = P_;
  }

  void DelayedEKF::Update(const VecX& z_true) {
    LUPNT_CHECK(!z_true.hasNaN(), "(Update) True measurement has NaN", "ASEKF");
    LUPNT_CHECK(f_meas_, "(Update) Measurement function not set", "ASEKF");

    // Store
    z_true_ = z_true;
    x_post_ = x_;
    P_post_ = P_;

    // Size
    int N_state = int(x_.size() / 2);
    int N_state_aug = int(x_.size());
    int N_meas = z_true_.size();
    if (N_meas == 0) return;  // no measurement, nothing to update

    z_prior_ = f_meas_(x_, &H_, &R_);
    LUPNT_CHECK(!z_prior_.hasNaN(), "(Update) Predicted measurement has NaN", "ASEKF");

    S_ = R_ + H_ * P_ * H_.transpose();  // Measurement information
    LUPNT_CHECK(!S_.hasNaN(), "(Update) Innovation covariance has NaN", "ASEKF");
    dz_ = z_true_ - z_prior_;
    LUPNT_CHECK(!dz_.hasNaN(), "(Update) Measurement residual has NaN", "ASEKF");

    // Outlier rejection
    if (use_custom_fault_detection_ && f_fault_det_) {
      f_fault_det_(this);
    } else {
      RemoveOutliers();
    }

    if (dz_.size() == 0) return;  // All measurements are outliers

    // Update step
    K_ = P_.block(0, 0, N_state, N_state_aug) * H_.transpose() * S_.inverse();  // (N_state, N_meas)
    LUPNT_CHECK(!K_.hasNaN(), "(Update) Kalman gain has NaN", "ASEKF");

    dx_ = K_ * dz_;  // State update (N_state,)
    Sigma_dx_ = K_ * S_ * K_.transpose();
    x_.head(N_state) = x_.head(N_state) + dx_;
    LUPNT_CHECK(!x_.hasNaN(), "(Update) State has NaN", "ASEKF");

    MatXd I = MatXd::Zero(N_state, N_state_aug);
    I.leftCols(N_state) = MatXd::Identity(N_state, N_state);

    MatXd G = I - K_ * H_;  // (N, 2N) - (N, M) * (M, 2N) = (N, 2N)
    P_.block(0, 0, N_state, N_state)
        = G * P_ * G.transpose()
          + K_ * R_ * K_.transpose();  // Joseph form (Only Update the first half)
    LUPNT_CHECK(!P_.hasNaN(), "(Update) Covariance has NaN", "ASEKF");
    LUPNT_CHECK((P_.diagonal().array() >= 0).all(), "(Update) Covariance has negative diagonal",
                "ASEKF");

    // Covariance inflation
    // double lambda = 0.0;
    // P_ = (1 + lambda) * P_;
    // P_ = G * P_;

    // Store
    x_post_ = x_;
  }
};  // namespace filtering_sim

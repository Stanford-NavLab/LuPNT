#include "src/udu_filter.h"

#include <matplot/util/colors.h>

#include "src/udu_utils.h"

namespace filtering_sim {
  using namespace lupnt;

  void UDUEKF::SetCovariance(const MatXd& P) {
    P_ = P;
    SetUd();  // Update U and D matrices based on the new covariance
  }

  void UDUEKF::SetUd() {
    int n = int(x_.size());
    std::tie(D_diag_, U_) = UDUDecomposition(P_);
  }

  void UDUEKF::Predict(Real t, const State* u) {
    // Checks
    LUPNT_CHECK(x_.size() > 0, "(Predict) State vector is not set", "UDU");
    LUPNT_CHECK(P_.size() > 0, "(Predict) Covariance matrix is not set", "UDU");
    LUPNT_CHECK(f_dyn_, "(Predict) Dynamics function not set", "UDU");
    LUPNT_CHECK(f_proc_, "(Predict) Process noise function not set", "UDU");

    int n = int(x_.size());
    int N_state = int(n / 2);  // For augmented state

    // Predict state
    x_ = f_dyn_(x_, t_, t, u, &F_);

    // Process noise
    Q_ = f_proc_(x_, t_, t);

    // Check if Q is diagonal
    LUPNT_CHECK(Q_.isDiagonal(), "(Predict) Process noise covariance is not diagonal", "UDU");

    if (!use_process_noise_mapping_) {
      G_ = MatXd::Identity(n, n);
    }
    LUPNT_CHECK(!Q_.hasNaN(), "(Predict) Process noise covariance has NaN", "UDU");

    LUPNT_CHECK(!x_.hasNaN(), "(Predict) Predicted state has NaN", "UDU");
    MatXd F_U = F_ * U_;
    int m = G_.cols();  // Process Noise dimension

    // Normal UDU Filter
    MatXd Y = MatXd::Zero(n, n + m);
    Y.block(0, 0, n, n) = F_U;
    Y.block(0, n, n, m) = G_;

    MatXd D_tilde = MatXd::Zero(n + m, n + m);
    D_tilde.block(0, 0, n, n) = D_diag_.asDiagonal();
    D_tilde.block(n, n, m, m) = Q_;

    // Update State
    std::tie(D_minus_diag_, U_minus_) = ModifiedGramSchmidt(D_tilde, Y);

    // Update time
    t_ = t;

    // UDU
    P_ = U_minus_ * D_minus_diag_.asDiagonal() * U_minus_.transpose();
    U_ = U_minus_;
    D_diag_ = D_minus_diag_;

    // Store
    x_prior_ = x_;
    P_prior_ = P_;
    // P_bar_ = F_U * D_diag_.asDiagonal() * F_U.transpose();  // Needed for adaptive filtering

    Sigma_dx_ = MatXd::Zero(n, n);  // Reset
    dx_ = VecXd::Zero(n);
  }

  void UDUEKF::CarlsonUpdate() {
    int n = int(x_.size());
    int N_meas = dz_.size();

    K_.resize(n, N_meas);
    dx_ = VecXd::Zero(n);
    D_plus_diag_ = VecXd::Zero(n);
    U_plus_ = MatXd::Identity(n, n);
    VecXd alpha(n + 1);
    VecXd beta(n);

    for (int m = 0; m < N_meas; m++) {
      VecXd Ht = H_.row(m);
      double R = R_(m, m);

      VecXd f = U_.transpose() * Ht;
      VecXd v = D_diag_.asDiagonal() * f;  // v*v = f^T * D * f = H^T * P * H

      alpha(0) = R;

      for (int k = 0; k < n; k++) {
        alpha(k + 1) = alpha(k) + v(k) * f(k);
        D_plus_diag_(k) = alpha(k) / alpha(k + 1) * D_diag_(k);
        beta(k) = v(k);
        if (k == 0) continue;
        double p_k = -f(k) / alpha(k);
        for (int j = 0; j < k; j++) {
          U_plus_(j, k) = U_(j, k) + beta(j) * p_k;
          beta(j) = beta(j) + U_(j, k) * v(k);
        }
      }

      // Set Kalman gain
      VecXd K_m = beta / alpha(n);  // N x 1
      K_.col(m) = K_m;

      if (N_meas > 1) {
        dz_(m) = dz_(m) - H_.row(m).dot(dx_);  // Modify residual for next measurement
      }

      // State Update
      VecXd dx_m = K_m * dz_(m);
      dx_ = dx_ + dx_m;

      // Update D_plus_diag_ and U_plus_ for the next measurement
      D_diag_ = D_plus_diag_;
      U_ = U_plus_;
    }

    // Update state
    x_ = x_ + dx_;

    // std::cout << "Kalman gain (K): \n" << K_.format(matFmt) << std::endl;

    // Construct Covariance
    P_ = U_ * D_diag_.asDiagonal() * U_.transpose();
  }

  void UDUEKF::Update(const VecX& z_true) {
    LUPNT_CHECK(!z_true.hasNaN(), "(Update) True measurement has NaN", "UDU");
    LUPNT_CHECK(f_meas_, "(Update)Measurement function not set", "UDU");

    // Store
    z_true_ = z_true;
    x_post_ = x_;
    P_post_ = P_;

    // Size
    int N_meas = z_true_.size();
    if (N_meas == 0) return;  // no measurement, nothing to update

    z_prior_ = f_meas_(x_, &H_, &R_);
    LUPNT_CHECK(!z_prior_.hasNaN(), "(Update) Predicted measurement has NaN", "UDU");

    dz_ = z_true_ - z_prior_;
    LUPNT_CHECK(!dz_.hasNaN(), "(Update) Measurement residual has NaN", "UDU");

    // Outlier rejection: Todo - Implement in UDU form
    S_ = R_ + H_ * P_ * H_.transpose();  // Measurement information

    if (use_custom_fault_detection_ && f_fault_det_) {
      f_fault_det_(this);
    } else {
      RemoveOutliers();
    }

    if (dz_.size() == 0) return;  // All measurements are outliers

    // Carlson rank 1 update
    CarlsonUpdate();

    LUPNT_CHECK(!x_.hasNaN(), "(Update) Updated state has NaN", "UDU");
    LUPNT_CHECK(!P_.hasNaN(), "(Update) Updated covariance has NaN", "UDU");

    // Store
    x_post_ = x_;
    P_post_ = P_;
    U_ = U_plus_;
    D_diag_ = D_plus_diag_;

    dx_ = x_post_ - x_prior_;
    Sigma_dx_ = P_post_ - P_prior_;
  }

  /*******************************************************
   * Smoother Functions
   *  For UDU, simply use the RTS smoother from EKF
   *  (the rank 2 updates requires inversion of the STM)
   ********************************************************/
  void UDUEKF::InitializeLogger(int max_tidx) {
    max_tidx_ = max_tidx;
    x_prior_log_.resize(max_tidx_);
    P_prior_log_.resize(max_tidx_);
    x_pos_log_.resize(max_tidx_);
    P_pos_log_.resize(max_tidx_);
    stm_log_.resize(max_tidx_);
    sm_tidx_ = max_tidx_ - 1;
    time_log_.resize(max_tidx_);
    x_sm_.resize(max_tidx_ + 1);
    P_sm_.resize(max_tidx_ + 1);
  }

  void UDUEKF::InitializeSmootherState() {
    if (x_sm_.size() != max_tidx_ + 1) {
      InitializeLogger(max_tidx_);
    }
    x_sm_[sm_tidx_] = x_post_;
    P_sm_[sm_tidx_] = P_post_;
    time_log_[sm_tidx_] = t_;
  }

  void UDUEKF::LogFilterEstimate(int tidx) {
    x_prior_log_[tidx] = x_prior_;
    P_prior_log_[tidx] = P_prior_;
    x_pos_log_[tidx] = x_post_;
    P_pos_log_[tidx] = P_post_;
    stm_log_[tidx] = F_;  // (tidx-1 -> tidx)
    time_log_[tidx] = t_;
  }

  // Smoother step Obtain x_tidx|N. and P_tidx|N
  void UDUEKF::UpdateSmoother(int tidx) {
    // Compute gain
    int k = tidx;
    MatXd P_kk = P_pos_log_[k];
    MatXd P_k1k = P_prior_log_[k + 1];
    MatXd stm_k = stm_log_[k + 1];
    VecXd x_kk = x_pos_log_[k];
    VecXd x_k1k = x_prior_log_[k + 1];

    MatXd A_k = P_kk * stm_k.transpose() * P_k1k.inverse();

    // Update smoothed state
    x_sm_[k] = x_kk + A_k * (x_sm_[k + 1] - x_k1k);
    // Update smoothed covariance
    P_sm_[k] = P_kk + A_k * (P_sm_[k + 1] - P_k1k) * A_k.transpose();
    sm_tidx_ = tidx;

    // Log time
    t_ = time_log_[k];
  }

  /*******************************************************************************************************
   * Delayed State UDU Smoother Update
   ******************************************************************************************************/
  void UDUDelayedEKF::StoreSubmatrices(int N_state) {
    U11_ = U_.topLeftCorner(N_state, N_state);
    U12_ = U_.topRightCorner(N_state, N_state);
    U22_ = U_.bottomRightCorner(N_state, N_state);
    D11_diag_ = D_diag_.head(N_state);
    D22_diag_ = D_diag_.tail(N_state);
    P11_ = P_.topLeftCorner(N_state, N_state);
    P12_ = P_.topRightCorner(N_state, N_state);
    P22_ = P_.bottomRightCorner(N_state, N_state);

    // std::cout << "Stored Submatrices:" << std::endl;
    // // In scientific notation for better readability
    // auto matFmt = Eigen::IOFormat(4, 0, ", ", "\n", "[", "]");
    // std::cout << "U11_: \n" << U11_.format(matFmt) << std::endl;
    // std::cout << "U12_: \n" << U12_.format(matFmt) << std::endl;
    // std::cout << "U22_: \n" << U22_.format(matFmt) << std::endl;
    // std::cout << "D11_diag_: \n" << D11_diag_.transpose().format(matFmt) << std::endl;
    // std::cout << "D22_diag_: \n" << D22_diag_.transpose().format(matFmt) << std::endl;
  }

  void UDUDelayedEKF::UpdateCurrStepUD() {
    // Compute U_ D_ U_^T = U11_ D11 U11_^T + U12 D22 U12_^T
    int n = int(U11_.rows());

    // Method 1: Rank-1 updates for each column of U12_ and corresponding D22_diag_ element
    // U_ = U11_;
    // D_diag_ = D11_diag_;
    // for (int i = 0; i < n; ++i) {
    //   std::tie(D_diag_, U_) = AgeeTurnerRankOneUpdate(U_, D_diag_, U12_.col(i), D22_diag_(i));
    //   std::cout << "D_diag_ after update " << i << ": \n" << D_diag_.transpose() << std::endl;
    // }

    // Method2: Direct reconstruction
    // MatXd P_11 = U11_ * D11_diag_.asDiagonal() * U11_.transpose()
    //             + U12_ * D22_diag_.asDiagonal() * U12_.transpose();
    MatXd P_11 = P_.topLeftCorner(n, n);
    std::tie(D_diag_, U_) = UDUDecomposition(P_11);

    // For check, compare against full covariance reconstruction
    // auto matFmt = Eigen::IOFormat(4, 0, ", ", "\n", "[", "]");
    // MatXd P_check = U_ * D_diag_.asDiagonal() * U_.transpose();
    // MatXd P_portion = P_.topLeftCorner(n, n);
    // std::cout << "P_check: \n" << P_check.format(matFmt) << std::endl;
    // std::cout << "P_portion: \n" << P_portion.format(matFmt) << std::endl;
    // std::cout << "P_check - P_portion: \n" << (P_check - P_portion).format(matFmt) << std::endl;
  }

  void UDUDelayedEKF::SetCovariance(const MatXd& P) {
    P_ = P;
    SetUd();  // Update U and D matrices based on the new covariance
  }

  void UDUDelayedEKF::SetUd() {
    int n = int(x_.size());
    int N_state = int(n / 2);  // For augmented state
    MatXd P_prev = P_.topLeftCorner(N_state, N_state);
    std::tie(D_diag_, U_) = UDUDecomposition(P_prev);
    // std::cout << "Initial U and D" << std::endl;
    // std::cout << "U: \n"
    //           << U_.format(Eigen::IOFormat(4, 0, ", ", "\n", "[", "]")) << std::endl;
    // std::cout << "D diag: \n"
    //           << D_diag_.format(Eigen::IOFormat(4, 0, ", ", "\n", "[", "]")) << std::endl;
  }

  void UDUDelayedEKF::Predict(Real t, const State* u) {
    // Checks
    LUPNT_CHECK(x_.size() > 0, "(Predict) State vector is not set", "UDU");
    LUPNT_CHECK(P_.size() > 0, "(Predict) Covariance matrix is not set", "UDU");
    LUPNT_CHECK(f_dyn_, "(Predict) Dynamics function not set", "UDU");
    LUPNT_CHECK(f_proc_, "(Predict) Process noise function not set", "UDU");

    int n = int(x_.size());
    int N_state = int(n / 2);  // For augmented state

    VecXd x_prev = x_.head(N_state);

    // Ensure the last half of x_ is the same
    x_.head(N_state) = f_dyn_(x_prev, t_, t, u, &F_);
    x_.tail(N_state) = x_prev;

    // Process Noise
    Q_ = f_proc_(x_prev, t_, t);

    // Check if Q is diagonal
    LUPNT_CHECK(Q_.isDiagonal(), "(Predict) Process noise covariance is not diagonal", "UDU");

    if (!use_process_noise_mapping_) {
      G_ = MatXd::Identity(n, n);
    }
    LUPNT_CHECK(!Q_.hasNaN(), "(Predict) Process noise covariance has NaN", "UDU");
    LUPNT_CHECK(!x_.hasNaN(), "(Predict) Predicted state has NaN", "UDU");

    MatXd F_U = F_ * U_;
    int m = G_.cols();  // Process Noise dimension

    // Construct Y and D matrices
    // Augmented State UDU Filter (For processing TDCP)
    MatXd Y = MatXd::Zero(N_state + N_state, m + N_state);
    VecXd D_tilde = VecXd::Zero(N_state + m);

    Y.topLeftCorner(N_state, m) = G_;
    Y.topRightCorner(N_state, N_state) = F_U;
    Y.bottomRightCorner(N_state, N_state) = U_;

    D_tilde.head(m) = Q_.diagonal();
    D_tilde.tail(N_state) = D_diag_;

    // Update State
    D_minus_diag_ = VecXd::Zero(N_state + N_state);
    D_minus_diag_.head(m + N_state) = D_tilde;

    // U_Minus
    //         N + m  N-m
    //   N+m   [       0]
    //   N-m   [   Y   I]
    U_minus_ = MatXd::Zero(N_state + N_state, N_state + N_state);
    U_minus_.topLeftCorner(N_state + N_state, m + N_state) = Y;
    if (m < N_state) {
      U_minus_.bottomRightCorner(N_state - m, N_state - m)
          = MatXd::Identity(N_state - m, N_state - m);
    }

    // Update time
    t_ = t;

    // UDU
    P_ = U_minus_ * D_minus_diag_.asDiagonal() * U_minus_.transpose();
    U_ = U_minus_;
    D_diag_ = D_minus_diag_;

    // Store the submatrices for smoother
    StoreSubmatrices(N_state);

    // Store
    x_prior_ = x_;
    P_prior_ = P_;
    // P_bar_ = F_U * D_diag_.asDiagonal() * F_U.transpose();  // Needed for adaptive filtering

    Sigma_dx_ = MatXd::Zero(n, n);  // Reset
    dx_ = VecXd::Zero(n);
  }

  void UDUDelayedEKF::DelayedCarlsonUpdate() {
    int n = int(x_.size());
    int N_meas = dz_.size();

    K_.resize(n, N_meas);
    dx_ = VecXd::Zero(n);
    D_plus_diag_ = VecXd::Zero(n);
    U_plus_ = MatXd::Identity(n, n);
    VecXd alpha(n + 1);
    VecXd beta(n);

    for (int m = 0; m < N_meas; m++) {
      VecXd Ht = H_.row(m);
      double R = R_(m, m);

      VecXd f = U_.transpose() * Ht;
      VecXd v = D_diag_.asDiagonal() * f;  // v*v = f^T * D * f = H^T * P * H

      alpha(0) = R;

      for (int k = 0; k < n; k++) {
        alpha(k + 1) = alpha(k) + v(k) * f(k);
        D_plus_diag_(k) = alpha(k) / alpha(k + 1) * D_diag_(k);
        beta(k) = v(k);
        if (k == 0) continue;
        double p_k = -f(k) / alpha(k);
        for (int j = 0; j < k; j++) {
          U_plus_(j, k) = U_(j, k) + beta(j) * p_k;
          beta(j) = beta(j) + U_(j, k) * v(k);
        }
      }

      // Set Kalman gain
      VecXd K_m = beta / alpha(n);  // N x 1
      K_.col(m) = K_m;

      if (N_meas > 1) {
        dz_(m) = dz_(m) - H_.row(m).dot(dx_);  // Modify residual for next measurement
      }

      // State Update
      VecXd dx_m = K_m * dz_(m);
      dx_ = dx_ + dx_m;

      // Update D_plus_diag_ and U_plus_ for the next measurement
      D_diag_ = D_plus_diag_;
      U_ = U_plus_;
    }

    // Update state
    x_ = x_ + dx_;

    // std::cout << "Kalman gain (K): \n" << K_.format(matFmt) << std::endl;

    // Construct Covariance
    P_ = U_ * D_diag_.asDiagonal() * U_.transpose();
  }

  void UDUDelayedEKF::Update(const VecX& z_true) {
    LUPNT_CHECK(!z_true.hasNaN(), "(Update) True measurement has NaN", "UDU");
    LUPNT_CHECK(f_meas_, "(Update)Measurement function not set", "UDU");

    // Store
    z_true_ = z_true;
    x_post_ = x_;
    P_post_ = P_;

    // Size
    int N_meas = z_true_.size();
    if (N_meas > 0) {
      z_prior_ = f_meas_(x_, &H_, &R_);
      LUPNT_CHECK(!z_prior_.hasNaN(), "(Update) Predicted measurement has NaN", "UDU");

      dz_ = z_true_ - z_prior_;
      LUPNT_CHECK(!dz_.hasNaN(), "(Update) Measurement residual has NaN", "UDU");

      // Outlier rejection: Todo - Implement in UDU form
      S_ = R_ + H_ * P_ * H_.transpose();  // Measurement information

      if (use_custom_fault_detection_ && f_fault_det_) {
        f_fault_det_(this);
      } else {
        RemoveOutliers();
      }
      if (dz_.size() > 0) {  // All measurements are outliers
        // Carlson rank 1 update
        DelayedCarlsonUpdate();
      }
    } else {
      U_plus_ = U_;
      D_plus_diag_ = D_diag_;
    }

    LUPNT_CHECK(!x_.hasNaN(), "(Update) Updated state has NaN", "UDU");
    LUPNT_CHECK(!P_.hasNaN(), "(Update) Updated covariance has NaN", "UDU");

    // Store
    x_post_ = x_;  // 2N
    P_post_ = P_;  // 2N x 2N

    int Nstate = int(x_.size() / 2);
    StoreSubmatrices(Nstate);  // Update submatrices after measurement update
    dx_ = x_post_.head(Nstate)
          - x_prior_.head(Nstate);  // Since there are multiple measurements, recompute dx_
    Sigma_dx_ = P_prior_.topLeftCorner(Nstate, Nstate)
                - P_post_.topLeftCorner(Nstate, Nstate);  // Compute change in covariance
    // Sigma_dx = K * S * K = PH^T * (H * P * H^T + R)^-1 * H * P
    // P+ = P - KHP- = (I - KH)P = P - PH^T (H * P * H^T + R)^-1 HP = P - Sigma_dx
    UpdateCurrStepUD();  // Update U and D for the current step
  }

  /*******************************************************
   * Smoother Functions (Delayed EKF)
   ********************************************************/
  void UDUDelayedEKF::InitializeLogger(int lent) {
    max_tidx_ = lent;
    x_pos_log_.resize(max_tidx_);
    U12_plus_log_.resize(max_tidx_);
    D22_plus_log_.resize(max_tidx_);
    U22_plus_log_.resize(max_tidx_);
    Ukk_log_.resize(max_tidx_);
    Dkk_log_.resize(max_tidx_);
    time_log_.resize(max_tidx_);

    P12_log_.resize(max_tidx_);
    P11_log_.resize(max_tidx_);
    P22_log_.resize(max_tidx_);

    // Smoother logs
    x_sm_.resize(max_tidx_);
    P_sm_.resize(max_tidx_);
  }

  void UDUDelayedEKF::InitializeSmootherState() {
    sm_tidx_ = max_tidx_ - 1;

    int n = int(x_.size());
    int N_state = int(n / 2);  // For augmented state

    x_sm_[sm_tidx_] = x_.head(N_state);
    P_sm_[sm_tidx_] = P_.topLeftCorner(N_state, N_state);
    time_log_[sm_tidx_] = t_;

    std::cout << std::scientific << std::setprecision(4) << std::endl;
    auto matFmt = Eigen::IOFormat(4, 0, ", ", "\n", "[", "]");

    // std::cout << "Initialized Delayed UDU Smoother State at tidx: " << sm_tidx_ << std::endl;
    // std::cout << "x_sm_: \n"
    //           << x_sm_[sm_tidx_].format(matFmt) << std::endl;
    // std::cout << "P_sm_: \n"
    //           << P_sm_[sm_tidx_].format(matFmt) << std::endl;
  }

  void UDUDelayedEKF::LogFilterEstimate(int tidx) {
    x_pos_log_[tidx] = x_post_;
    U12_plus_log_[tidx] = U12_;
    D22_plus_log_[tidx] = D22_diag_;
    U22_plus_log_[tidx] = U22_;
    P12_log_[tidx] = P12_;
    P11_log_[tidx] = P11_;
    P22_log_[tidx] = P22_;
    Ukk_log_[tidx] = U_;
    Dkk_log_[tidx] = D_diag_;
    time_log_[tidx] = t_;
  }

  void UDUDelayedEKF::UpdateSmoother(int tidx) {
    int n = int(x_.size());
    int N_state = int(n / 2);  // For augmented state

    int use_state = N_state;  // Do not update last two states (clock drift rate and srp)

    // Extract prior and posterior states
    VecXd x_kp1kp1 = x_pos_log_[tidx + 1].head(N_state).head(use_state);
    VecXd x_kkp1 = x_pos_log_[tidx + 1].tail(N_state).head(use_state);

    MatXd U12_plus = U12_plus_log_[tidx + 1].block(0, 0, use_state, use_state);
    VecXd D22_plus_diag = D22_plus_log_[tidx + 1].head(use_state);
    MatXd U22_plus = U22_plus_log_[tidx + 1].block(0, 0, use_state, use_state);
    MatXd Ukp1kp1 = Ukk_log_[tidx + 1].block(0, 0, use_state, use_state);
    VecXd Dkp1kp1 = Dkk_log_[tidx + 1].head(use_state);
    MatXd P_kkp1 = P22_log_[tidx + 1].block(0, 0, use_state, use_state);

    // Smoother gain with forward/backward substitution
    // Currently disabled since it generates NaNs in some cases
    // MatXd Y = BackwardSubstitution(Ukp1kp1, U12_plus);
    // MatXd Z = Y;
    // for (int i = 0; i < Z.rows(); ++i) {
    //   for (int j = 0; j < Z.cols(); ++j) {
    //     Z(i, j) *= D22_plus_diag(j) / Dkp1kp1(i);
    //   }
    // }
    // MatXd J_kp1 = ForwardSubstitution(Ukp1kp1, Z * U22_plus.transpose()).transpose();

    MatXd P_kp1kp1 = P11_log_[tidx + 1].block(0, 0, use_state, use_state);
    MatXd P_cross = P12_log_[tidx + 1].block(0, 0, use_state, use_state);
    MatXd P_kp1kp1_inv = P_kp1kp1.inverse();
    MatXd J_kp1 = P_cross.transpose() * P_kp1kp1_inv;
    VecXd diff_x = x_sm_[tidx + 1].head(use_state) - x_kp1kp1;

    // Do not update the last two states (clock drift rate and srp)
    VecXd dx = J_kp1 * (x_sm_[tidx + 1].head(use_state) - x_kp1kp1);

    x_sm_[tidx] = x_pos_log_[tidx + 1].tail(N_state);  // x_(k|k+1)
    x_sm_[tidx].segment(0, use_state) += dx;

    // Update smoothed covariance
    P_sm_[tidx] = P22_log_[tidx + 1];
    MatXd diff_P = (P_sm_[tidx + 1].block(0, 0, use_state, use_state) - P_kp1kp1);
    P_sm_[tidx] += J_kp1 * diff_P * J_kp1.transpose();

    // Update smoothed state
    // std::cout << std::scientific << std::setprecision(4) << std::endl;
    // std::cout << "Smoother Update at tidx: " << tidx << std::endl;
    // auto fmt = Eigen::IOFormat(4, 0, ", ", "\n", "[", "]");
    // std::cout << "P_cross: \n" << P_cross.format(fmt) << std::endl;
    // std::cout << "P_kp1kp1_inv: \n" << P_kp1kp1_inv.format(fmt) << std::endl;
    // std::cout << "J_kp1_direct: \n" << J_kp1_direct.format(fmt) << std::endl;
    // std::cout << "J_kp1: \n" << J_kp1.format(fmt) << std::endl;
    // std::cout << "diff_J: \n" << (J_kp1 - J_kp1_direct).format(fmt) << std::endl;
    // std::cout << "diff_x: \n" << diff_x.transpose().format(fmt) << std::endl;
    // std::cout << "x_kp1kp1: \n" << x_kp1kp1.transpose().format(fmt) << std::endl;
    // std::cout << "x_kkp1: \n" << x_kkp1.transpose().format(fmt) << std::endl;
    // std::cout << "x_sm_[tidx + 1]: \n" << x_sm_[tidx + 1].transpose().format(fmt) << std::endl;
    // std::cout << "dx: \n" << dx.transpose().format(fmt) << std::endl;
    // std::cout << "diff_P: \n" << diff_P.format(fmt) << std::endl;
    // std::cout << "P_sm_[tidx]: \n" << P_sm_[tidx].format(fmt) << std::endl;

    LUPNT_CHECK(!x_sm_[tidx].hasNaN(), "(Smoother Update) Smoothed state has NaN", "UDU");
    LUPNT_CHECK(!P_sm_[tidx].hasNaN(), "(Smoother Update) Smoothed covariance has NaN", "UDU");
    LUPNT_CHECK((P_sm_[tidx].diagonal().array() >= 0).all(),
                "(Update) Covariance has negative diagonal", "UDU");
    // Log time
    sm_tidx_ = tidx;
    t_ = time_log_[tidx];
  }

};  // namespace filtering_sim

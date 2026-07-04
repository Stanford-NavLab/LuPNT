#include "lupnt/filters/udu_filter.h"

#include "lupnt/core/asset_factory.h"
#include "lupnt/filters/udu_utils.h"

namespace lupnt {

  void UDUEKF::SetUd() { std::tie(D_diag_, U_) = UDUDecomposition(P_); }

  void UDUEKF::SetCovariance(const MatXd& P) {
    Filter::SetCovariance(P);
    SetUd();
  }

  void UDUEKF::Predict(Real t, const State* u) {
    LUPNT_CHECK(x_.size() > 0, "(Predict) State vector is not set", "UDU");
    LUPNT_CHECK(P_.size() > 0, "(Predict) Covariance matrix is not set", "UDU");
    LUPNT_CHECK(f_dyn_, "(Predict) Dynamics function not set", "UDU");
    LUPNT_CHECK(f_proc_, "(Predict) Process noise function not set", "UDU");

    Q_ = f_proc_(x_, t_, t);
    if (use_process_noise_mapping_) {
      Q_ = G_ * Q_ * G_.transpose();
    }
    LUPNT_CHECK(!Q_.hasNaN(), "(Predict) Process noise covariance has NaN", "UDU");

    x_ = f_dyn_(x_, t_, t, u, &F_);
    LUPNT_CHECK(!x_.hasNaN(), "(Predict) State has NaN", "UDU");
    LUPNT_CHECK(F_.cols() == x_.rows() && F_.rows() == x_.rows(),
                "(Predict) State and STM have different sizes", "UDU");

    P_bar_ = F_ * UDUReconstruct(U_, D_diag_) * F_.transpose();
    P_ = P_bar_ + Q_;
    P_ = 0.5 * (P_ + P_.transpose());
    LUPNT_CHECK(!P_.hasNaN(), "(Predict) Covariance has NaN", "UDU");
    SetUd();

    t_ = t;
    x_prior_ = x_;
    P_prior_ = P_;
    dx_ = VecXd::Zero(x_.size());
    Sigma_dx_ = MatXd::Zero(x_.size(), x_.size());
  }

  void UDUEKF::CarlsonUpdate() {
    const int n = static_cast<int>(x_.size());
    const int n_meas = static_cast<int>(dz_.size());

    K_.setZero(n, n_meas);
    dx_ = VecXd::Zero(n);

    for (int meas = 0; meas < n_meas; ++meas) {
      VecXd h = H_.row(meas).transpose();
      double r = R_(meas, meas);

      MatXd P_curr = UDUReconstruct(U_, D_diag_);
      Real innovation_var = r + (h.transpose() * P_curr * h)(0, 0);
      LUPNT_CHECK(innovation_var > 0.0, "(Update) Innovation variance must be positive", "UDU");

      VecXd k = (P_curr * h / innovation_var.val()).cast<double>();
      K_.col(meas) = k;

      Real residual = dz_(meas);
      if (meas > 0) residual -= H_.row(meas).dot(dx_);
      dx_ += k * residual.val();

      MatXd I = MatXd::Identity(n, n);
      MatXd G = I - k * h.transpose();
      MatXd P_next = G * P_curr * G.transpose() + (k * r) * k.transpose();
      P_next = 0.5 * (P_next + P_next.transpose());
      std::tie(D_diag_, U_) = UDUDecomposition(P_next);
      P_ = P_next;
    }

    x_ = x_ + dx_;
  }

  void UDUEKF::Update(const VecX& z_true) {
    LUPNT_CHECK(!z_true.hasNaN(), "(Update) True measurement has NaN", "UDU");
    LUPNT_CHECK(f_meas_, "(Update) Measurement function not set", "UDU");

    z_true_ = z_true;
    x_post_ = x_;
    P_post_ = P_;

    const int n_meas = static_cast<int>(z_true_.size());
    if (n_meas == 0) return;

    z_prior_ = f_meas_(x_, &H_, &R_);
    LUPNT_CHECK(!z_prior_.hasNaN(), "(Update) Predicted measurement has NaN", "UDU");

    S_ = R_ + H_ * P_ * H_.transpose();
    LUPNT_CHECK(!S_.hasNaN(), "(Update) Innovation covariance has NaN", "UDU");

    dz_ = z_true_ - z_prior_;
    LUPNT_CHECK(!dz_.hasNaN(), "(Update) Measurement residual has NaN", "UDU");

    if (use_custom_fault_detection_ && f_fault_det_) {
      f_fault_det_(this);
    } else {
      RemoveOutliers();
    }
    if (dz_.size() == 0) return;

    CarlsonUpdate();

    LUPNT_CHECK(!x_.hasNaN(), "(Update) Updated state has NaN", "UDU");
    LUPNT_CHECK(!P_.hasNaN(), "(Update) Updated covariance has NaN", "UDU");
    LUPNT_CHECK((P_.diagonal().array() >= 0).all(), "(Update) Covariance has negative diagonal",
                "UDU");

    P_ = UDUReconstruct(U_, D_diag_);
    x_post_ = x_;
    P_post_ = P_;
    Sigma_dx_ = P_post_ - P_prior_;
  }

  REGISTER_FACTORY_CLASS(Filter, UDUEKF)

  void UDUStochasticCloningEKF::InferBaseStateSize() {
    if (base_state_size_ == 0 && x_.size() > 0) {
      LUPNT_CHECK(x_.size() % 2 == 0, "Cloned UDU state must have even size", "UDU");
      base_state_size_ = static_cast<int>(x_.size() / 2);
    }
  }

  void UDUStochasticCloningEKF::SetBaseStateSize(int base_state_size) {
    LUPNT_CHECK(base_state_size > 0, "Base state size must be positive", "UDU");
    base_state_size_ = base_state_size;
  }

  void UDUStochasticCloningEKF::SetState(const State& x) {
    Filter::SetState(x);
    InferBaseStateSize();
    LUPNT_CHECK(x_.size() == 2 * base_state_size_,
                "Cloned UDU state size must be twice the base state size", "UDU");
  }

  void UDUStochasticCloningEKF::SetCovariance(const MatXd& P) {
    Filter::SetCovariance(P);
    InferBaseStateSize();
    LUPNT_CHECK(P_.rows() == 2 * base_state_size_ && P_.cols() == 2 * base_state_size_,
                "Cloned UDU covariance size must be twice the base state size", "UDU");
    SetUd();
  }

  void UDUStochasticCloningEKF::Predict(Real t, const State* u) {
    LUPNT_CHECK(x_.size() > 0, "(Predict) State vector is not set", "UDU");
    LUPNT_CHECK(P_.size() > 0, "(Predict) Covariance matrix is not set", "UDU");
    LUPNT_CHECK(f_dyn_, "(Predict) Dynamics function not set", "UDU");
    LUPNT_CHECK(f_proc_, "(Predict) Process noise function not set", "UDU");
    InferBaseStateSize();

    const int n = base_state_size_;
    LUPNT_CHECK(x_.size() == 2 * n, "Cloned UDU state has incorrect size", "UDU");

    State x_prev_current = x_.head(n);
    x_prev_current.SetFrame(x_.GetFrame());  // State's Eigen-slice ctor leaves frame_ unset
    MatXd F_base;
    State x_current = f_dyn_(x_prev_current, t_, t, u, &F_base);
    LUPNT_CHECK(F_base.rows() == n && F_base.cols() == n, "(Predict) Base STM has incorrect size",
                "UDU");

    State x_pred(2 * n);
    x_pred.head(n) = x_current;
    x_pred.tail(n) = x_prev_current;
    x_pred.SetFrame(x_.GetFrame());
    x_ = x_pred;

    Q_ = f_proc_(x_prev_current, t_, t);
    if (use_process_noise_mapping_) Q_ = G_ * Q_ * G_.transpose();
    LUPNT_CHECK(Q_.rows() == n && Q_.cols() == n, "(Predict) Base process noise has incorrect size",
                "UDU");
    LUPNT_CHECK(!Q_.hasNaN(), "(Predict) Process noise covariance has NaN", "UDU");

    F_.setZero(2 * n, 2 * n);
    F_.topLeftCorner(n, n) = F_base;
    F_.bottomLeftCorner(n, n) = MatXd::Identity(n, n);

    MatXd Q_aug = MatXd::Zero(2 * n, 2 * n);
    Q_aug.topLeftCorner(n, n) = Q_;

    P_bar_ = F_ * UDUReconstruct(U_, D_diag_) * F_.transpose();
    P_ = P_bar_ + Q_aug;
    P_ = 0.5 * (P_ + P_.transpose());
    LUPNT_CHECK(!x_.hasNaN(), "(Predict) Predicted cloned state has NaN", "UDU");
    LUPNT_CHECK(!P_.hasNaN(), "(Predict) Cloned covariance has NaN", "UDU");
    SetUd();

    t_ = t;
    x_prior_ = x_;
    P_prior_ = P_;
    dx_ = VecXd::Zero(x_.size());
    Sigma_dx_ = MatXd::Zero(x_.size(), x_.size());
  }

  REGISTER_FACTORY_CLASS(Filter, UDUStochasticCloningEKF)

}  // namespace lupnt

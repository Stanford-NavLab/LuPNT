#include "lupnt/numerics/filters/srif.h"

#include "lupnt/core/asset_factory.h"

namespace lupnt {

  namespace {

    /// Upper-triangular R factor of a Householder QR of `A` (p x q, p >= q), size q x q.
    /// Any right-hand-side columns appended to `A` are triangularized along with it, so
    /// the top rows of those columns come out as `Q^T rhs` -- the transformed data.
    MatXd QrRFactor(const MatXd& A) {
      const int q = static_cast<int>(A.cols());
      Eigen::HouseholderQR<MatXd> qr(A);
      return qr.matrixQR().topRows(q).template triangularView<Eigen::Upper>();
    }

  }  // namespace

  void SRIF::SetInfoSqrt() {
    Eigen::LLT<MatXd> llt(P_.inverse());
    LUPNT_CHECK(llt.info() == Eigen::Success,
                "(SetInfoSqrt) Inverse covariance is not positive definite", "SRIF");
    // R_info_^T R_info_ = P^{-1}, upper triangular: R = L^T with P^{-1} = L L^T.
    R_info_ = MatXd(llt.matrixL().transpose());
  }

  void SRIF::SetCovariance(const MatXd& P) {
    Filter::SetCovariance(P);
    SetInfoSqrt();
  }

  void SRIF::Predict(Real t, const State* u) {
    LUPNT_CHECK(x_.size() > 0, "(Predict) State vector is not set", "SRIF");
    LUPNT_CHECK(P_.size() > 0, "(Predict) Covariance matrix is not set", "SRIF");
    LUPNT_CHECK(f_dyn_, "(Predict) Dynamics function not set", "SRIF");
    LUPNT_CHECK(f_proc_, "(Predict) Process noise function not set", "SRIF");

    // Process noise (evaluated at the prior state, same convention as EKF::Predict).
    Q_ = f_proc_(x_, t_, t);
    LUPNT_CHECK(!Q_.hasNaN(), "(Predict) Process noise covariance has NaN", "SRIF");
    if (use_process_noise_mapping_) {
      Q_ = G_ * Q_ * G_.transpose();
      LUPNT_CHECK(!Q_.hasNaN(), "(Predict) Mapped process noise covariance has NaN", "SRIF");
    }

    // Dynamics (state + state-transition matrix F).
    x_ = f_dyn_(x_, t_, t, u, &F_);
    LUPNT_CHECK(!x_.hasNaN(), "(Predict) State has NaN", "SRIF");
    LUPNT_CHECK(F_.cols() == x_.rows() && F_.rows() == x_.rows(),
                "(Predict) State and STM have different sizes", "SRIF");

    const int n = static_cast<int>(x_.size());
    MatXd Phi_inv = F_.inverse();
    P_bar_ = F_ * P_ * F_.transpose();  // F P F^T (before process noise)

    // Process-noise square root S (Q = S S^T); zero when there is no process noise.
    MatXd S = MatXd::Zero(n, n);
    if (!Q_.isZero()) {
      Eigen::LLT<MatXd> llt(Q_);
      LUPNT_CHECK(llt.info() == Eigen::Success,
                  "(Predict) Process noise covariance is not positive definite", "SRIF");
      S = llt.matrixL();
    }

    // Dyer-McReynolds time-update array (the information estimate is re-centered on x_
    // each step, so the data column b is zero):
    //   [ I              0          ]
    //   [ -R Phi^{-1} S  R Phi^{-1} ]
    // A full QR triangularization zeros the bottom-left block; the bottom-right block is
    // the a-priori information square root at time t.
    MatXd RPhiInv = R_info_ * Phi_inv;
    MatXd A = MatXd::Zero(2 * n, 2 * n);
    A.topLeftCorner(n, n) = MatXd::Identity(n, n);
    A.bottomLeftCorner(n, n) = -RPhiInv * S;
    A.bottomRightCorner(n, n) = RPhiInv;
    R_info_ = QrRFactor(A).bottomRightCorner(n, n);

    MatXd R_inv = R_info_.inverse();
    P_ = R_inv * R_inv.transpose();
    LUPNT_CHECK(!P_.hasNaN(), "(Predict) Covariance has NaN", "SRIF");

    t_ = t;
    x_prior_ = x_;
    P_prior_ = P_;
    dx_ = VecXd::Zero(n);
    Sigma_dx_ = MatXd::Zero(n, n);
  }

  void SRIF::InformationUpdate() {
    const int n = static_cast<int>(x_.size());
    const int m = static_cast<int>(dz_.size());

    // Whiten the measurements by R^{-1/2}: with R = L L^T, the data equation
    // L^{-1} dz = (L^{-1} H) dx + unit-variance noise.
    Eigen::LLT<MatXd> llt(R_);
    LUPNT_CHECK(llt.info() == Eigen::Success,
                "(Update) Measurement noise covariance is not positive definite", "SRIF");
    MatXd L = llt.matrixL();
    MatXd H_w = L.template triangularView<Eigen::Lower>().solve(H_);    // m x n
    VecXd dz_w = L.template triangularView<Eigen::Lower>().solve(dz_);  // m

    // QR of the measurement-update array (the data column b is zero because the estimate
    // is re-centered on x_):
    //   [ R_info_  0    ]
    //   [ H_w      dz_w ]
    // The triangularized result gives the updated information square root R_hat (top-left)
    // and the correction data b_hat (top of the last column).
    MatXd A(n + m, n + 1);
    A.topLeftCorner(n, n) = R_info_;
    A.block(0, n, n, 1).setZero();
    A.bottomLeftCorner(m, n) = H_w;
    A.block(n, n, m, 1) = dz_w;

    MatXd Rtri = QrRFactor(A);
    R_info_ = Rtri.topLeftCorner(n, n);
    VecXd b_hat = Rtri.block(0, n, n, 1);

    dx_ = R_info_.template triangularView<Eigen::Upper>().solve(b_hat);
    x_ = x_ + dx_;

    MatXd R_inv = R_info_.inverse();
    P_ = R_inv * R_inv.transpose();
  }

  void SRIF::Update(const VecX& z_true) {
    LUPNT_CHECK(!z_true.hasNaN(), "(Update) True measurement has NaN", "SRIF");
    LUPNT_CHECK(f_meas_, "(Update) Measurement function not set", "SRIF");

    z_true_ = z_true;
    x_post_ = x_;
    P_post_ = P_;

    const int n_meas = static_cast<int>(z_true_.size());
    if (n_meas == 0) return;  // no measurement, nothing to update

    z_prior_ = f_meas_(x_, &H_, &R_);
    LUPNT_CHECK(!z_prior_.hasNaN(), "(Update) Predicted measurement has NaN", "SRIF");

    S_ = R_ + H_ * P_ * H_.transpose();
    LUPNT_CHECK(!S_.hasNaN(), "(Update) Innovation covariance has NaN", "SRIF");

    dz_ = z_true_ - z_prior_;
    LUPNT_CHECK(!dz_.hasNaN(), "(Update) Measurement residual has NaN", "SRIF");

    // Outlier rejection (reuses EKF's residual-ratio test / custom fault detection).
    if (use_custom_fault_detection_ && f_fault_det_) {
      f_fault_det_(this);
    } else {
      RemoveOutliers();
    }
    if (dz_.size() == 0) return;  // all measurements are outliers

    InformationUpdate();

    LUPNT_CHECK(!x_.hasNaN(), "(Update) Updated state has NaN", "SRIF");
    LUPNT_CHECK(!P_.hasNaN(), "(Update) Updated covariance has NaN", "SRIF");
    LUPNT_CHECK((P_.diagonal().array() >= 0).all(), "(Update) Covariance has negative diagonal",
                "SRIF");

    x_post_ = x_;
    P_post_ = P_;
    Sigma_dx_ = P_post_ - P_prior_;
  }

  MatXd CwnaProcessNoise(double dt, double accel_psd) {
    MatXd Q = MatXd::Zero(6, 6);
    const double q3 = accel_psd * dt * dt * dt / 3.0;
    const double q2 = accel_psd * dt * dt / 2.0;
    const double q1 = accel_psd * dt;
    for (int k = 0; k < 3; ++k) {
      Q(k, k) = q3;
      Q(k, k + 3) = q2;
      Q(k + 3, k) = q2;
      Q(k + 3, k + 3) = q1;
    }
    return Q;
  }

  REGISTER_FACTORY_CLASS(Filter, SRIF)

}  // namespace lupnt

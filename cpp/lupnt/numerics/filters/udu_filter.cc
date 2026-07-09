#include "lupnt/numerics/filters/udu_filter.h"

#include "lupnt/core/asset_factory.h"
#include "lupnt/numerics/filters/udu_utils.h"

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

    // Process noise, evaluated at the prior state (same convention as EKF::Predict). The UDU
    // time update requires a *diagonal* Q: correlated process noise (e.g. state-noise-
    // compensation position/velocity, or the 2-state clock block) must be supplied through
    // the mapping matrix G (Q = G diag(Q) G^T) via SetProcessNoiseMappingMatrix, exactly as
    // the reference filter does. With no mapping, G defaults to identity and Q must already be
    // diagonal.
    Q_ = f_proc_(x_, t_, t);
    LUPNT_CHECK(!Q_.hasNaN(), "(Predict) Process noise covariance has NaN", "UDU");

    x_ = f_dyn_(x_, t_, t, u, &F_);
    LUPNT_CHECK(!x_.hasNaN(), "(Predict) State has NaN", "UDU");
    LUPNT_CHECK(F_.cols() == x_.rows() && F_.rows() == x_.rows(),
                "(Predict) State and STM have different sizes", "UDU");

    const int n = static_cast<int>(x_.size());
    if (!use_process_noise_mapping_) G_ = MatXd::Identity(n, n);
    LUPNT_CHECK(Q_.isDiagonal(),
                "(Predict) UDU process noise must be diagonal; supply correlated noise through "
                "the mapping matrix G (SetProcessNoiseMappingMatrix)",
                "UDU");
    const int m = static_cast<int>(G_.cols());  // reduced process-noise dimension
    LUPNT_CHECK(G_.rows() == n && Q_.rows() == m && Q_.cols() == m,
                "(Predict) Process-noise mapping G must be n x m with a matching m x m Q", "UDU");

    P_bar_ = F_ * UDUReconstruct(U_, D_diag_) * F_.transpose();  // F P F^T (for adaptive filtering)

    // Numerically stable UDU time update (Thornton MWGS): with Y = [F U | G] and
    // D_tilde = blkdiag(D, Q), Y D_tilde Y^T = F P F^T + G Q G^T is factored *without ever
    // forming or re-decomposing P*, keeping the (diagonal) process-noise weight separate from
    // the state factor. Reconstructing P, adding Q, and re-factoring (the previous
    // implementation) rounds variances many orders of magnitude below the largest one down to
    // zero -- e.g. the clock-drift term collapses after a few steps.
    MatXd Y = MatXd::Zero(n, n + m);
    Y.leftCols(n) = F_ * U_;
    Y.rightCols(m) = G_;
    MatXd D_tilde = MatXd::Zero(n + m, n + m);
    D_tilde.topLeftCorner(n, n) = D_diag_.asDiagonal();
    D_tilde.bottomRightCorner(m, m) = Q_;
    std::tie(D_diag_, U_) = ModifiedGramSchmidt(D_tilde, Y);

    P_ = UDUReconstruct(U_, D_diag_);
    LUPNT_CHECK(!P_.hasNaN(), "(Predict) Covariance has NaN", "UDU");

    t_ = t;
    x_prior_ = x_;
    P_prior_ = P_;
    dx_ = VecXd::Zero(n);
    Sigma_dx_ = MatXd::Zero(n, n);
  }

  void UDUEKF::CarlsonUpdate() {
    const int n = static_cast<int>(x_.size());
    const int n_meas = static_cast<int>(dz_.size());

    K_.setZero(n, n_meas);
    dx_ = VecXd::Zero(n);

    // Carlson rank-1 (Bierman-Thornton) sequential update: each scalar measurement updates the
    // U/D factors in place via the alpha/beta recursion, so the covariance is never
    // reconstructed or re-decomposed. This preserves variances spanning many orders of
    // magnitude, unlike a reconstruct/Joseph/re-factor update (which rounds the smallest ones
    // to zero). Requires a diagonal R (measurements processed one component at a time).
    VecXd alpha(n + 1);
    VecXd beta(n);
    MatXd U_plus = MatXd::Identity(n, n);
    VecXd D_plus = VecXd::Zero(n);

    for (int meas = 0; meas < n_meas; ++meas) {
      VecXd h = H_.row(meas).transpose();
      double r = R_(meas, meas);

      VecXd f = U_.transpose() * h;        // f = U^T h
      VecXd v = D_diag_.asDiagonal() * f;  // v = D f, so f.dot(v) = h^T P h

      alpha(0) = r;
      for (int k = 0; k < n; ++k) {
        alpha(k + 1) = alpha(k) + v(k) * f(k);  // alpha grows from r > 0, stays positive
        D_plus(k) = alpha(k) / alpha(k + 1) * D_diag_(k);
        beta(k) = v(k);
        if (k == 0) continue;
        double p_k = -f(k) / alpha(k);
        for (int j = 0; j < k; ++j) {
          U_plus(j, k) = U_(j, k) + beta(j) * p_k;
          beta(j) = beta(j) + U_(j, k) * v(k);
        }
      }

      VecXd k_gain = beta / alpha(n);  // Kalman gain for this scalar measurement
      K_.col(meas) = k_gain;

      Real residual = dz_(meas);
      if (meas > 0) residual -= H_.row(meas).dot(dx_);  // correct residual for prior updates
      dx_ += k_gain * residual.val();

      D_diag_ = D_plus;
      U_ = U_plus;
    }

    x_ = x_ + dx_;
    P_ = UDUReconstruct(U_, D_diag_);
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

    // Base-state process noise (diagonal, correlated noise mapped through the base G -- see
    // UDUEKF::Predict). G is the n x m base mapping; with no mapping it defaults to identity.
    Q_ = f_proc_(x_prev_current, t_, t);
    LUPNT_CHECK(!Q_.hasNaN(), "(Predict) Process noise covariance has NaN", "UDU");
    if (!use_process_noise_mapping_) G_ = MatXd::Identity(n, n);
    LUPNT_CHECK(Q_.isDiagonal(),
                "(Predict) UDU process noise must be diagonal; supply correlated noise through "
                "the mapping matrix G (SetProcessNoiseMappingMatrix)",
                "UDU");
    const int m = static_cast<int>(G_.cols());  // reduced base process-noise dimension
    LUPNT_CHECK(G_.rows() == n && Q_.rows() == m && Q_.cols() == m,
                "(Predict) Base process-noise mapping G must be n x m with a matching m x m Q",
                "UDU");

    const int n2 = 2 * n;
    F_.setZero(n2, n2);
    F_.topLeftCorner(n, n) = F_base;
    F_.bottomLeftCorner(n, n) = MatXd::Identity(n, n);

    P_bar_ = F_ * UDUReconstruct(U_, D_diag_) * F_.transpose();  // (for adaptive filtering)

    // Augmented Thornton MWGS time update: the process noise acts on the current (top) block
    // only, mapped through the base G. With Y = [F U | E], E = [G; 0], and
    // D_tilde = blkdiag(D, Q), we get Y D_tilde Y^T = F P F^T + blkdiag(G Q G^T, 0), factored
    // without forming/re-decomposing the full 2n covariance (which would round the smallest
    // variances to zero).
    MatXd E = MatXd::Zero(n2, m);
    E.topRows(n) = G_;
    MatXd Y = MatXd::Zero(n2, n2 + m);
    Y.leftCols(n2) = F_ * U_;
    Y.rightCols(m) = E;
    MatXd D_tilde = MatXd::Zero(n2 + m, n2 + m);
    D_tilde.topLeftCorner(n2, n2) = D_diag_.asDiagonal();
    D_tilde.bottomRightCorner(m, m) = Q_;
    std::tie(D_diag_, U_) = ModifiedGramSchmidt(D_tilde, Y);

    P_ = UDUReconstruct(U_, D_diag_);
    LUPNT_CHECK(!x_.hasNaN(), "(Predict) Predicted cloned state has NaN", "UDU");
    LUPNT_CHECK(!P_.hasNaN(), "(Predict) Cloned covariance has NaN", "UDU");

    t_ = t;
    x_prior_ = x_;
    P_prior_ = P_;
    dx_ = VecXd::Zero(x_.size());
    Sigma_dx_ = MatXd::Zero(x_.size(), x_.size());
  }

  REGISTER_FACTORY_CLASS(Filter, UDUStochasticCloningEKF)

}  // namespace lupnt

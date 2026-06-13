#include "lupnt/filters/adaptive_process_noise.h"

#include "lupnt/core/definitions.h"
#include "lupnt/filters/filter_utils.h"

namespace lupnt {

  void ProcessNoise::Reset() {
    dx_list_.clear();
    P_post_list_.clear();
    P_bar_list_.clear();
    Sigma_dx_list_.clear();
    step_ = 0;
  }

  void ProcessNoise::Update(const VecXd& dx, const MatXd& Sigma_dx, const MatXd& P_bar,
                            const MatXd& P_post) {
    // LUPNT_CHECK((Sigma_dx.diagonal().array() >= 0).all(), "Sigma_dx has negative diagonal",
    //             "AdaptiveModel");
    // LUPNT_CHECK((P_bar.diagonal().array() >= 0).all(), "P_bar has negative diagonal",
    //             "AdaptiveModel");
    // LUPNT_CHECK((P_post.diagonal().array() >= 0).all(), "P_post has negative diagonal",
    //             "AdaptiveModel");
    if ((Sigma_dx.diagonal().array() < 0).any()) {
      Logger::Warn("Sigma_dx has negative diagonal elements, skipping update", "AdaptiveModel");
      Reset();
      return;
    }
    if ((P_bar.diagonal().array() < 0).any()) {
      Logger::Warn("P_bar has negative diagonal elements, skipping update", "AdaptiveModel");
      Reset();
      return;
    }
    if ((P_post.diagonal().array() < 0).any()) {
      Logger::Warn("P_post has negative diagonal elements, skipping update", "AdaptiveModel");
      Reset();
      return;
    }

    if (algorithm_ == ProcessNoiseAlgorithm::SNC) return;

    step_++;
    if (step_ < N_wait_) return;
    if (step_ == N_wait_) Logger::Info("Starting adaptive filtering", "AdaptiveModel");

    dx_list_.push_back(dx);
    P_post_list_.push_back(P_post);
    P_bar_list_.push_back(P_bar);
    Sigma_dx_list_.push_back(Sigma_dx);

    if (dx_list_.size() > N_window_) {
      dx_list_.pop_front();
      P_post_list_.pop_front();
      P_bar_list_.pop_front();
      Sigma_dx_list_.pop_front();
    }
  }

  std::tuple<MatXd, MatXd> ProcessNoise::ComputeSlidingWindowStats() {
    int N_steps = dx_list_.size();
    int N_state = dx_list_[0].size();

    MatXd Q_hat = MatXd::Zero(N_state, N_state);
    MatXd Sigma_bar_sum = MatXd::Zero(N_state, N_state);

    LUPNT_CHECK(N_steps == N_window_, "Window size must be equal to the number of steps",
                "AdaptiveModel");

    for (int i = 0; i < N_steps; ++i) {
      const auto& P_post = P_post_list_[i];
      const auto& P_bar = P_bar_list_[i];
      const auto& dx = dx_list_[i];
      const auto& Sigma_dx = Sigma_dx_list_[i];

      Q_hat += P_post - P_bar + dx * dx.transpose();
      VecXd Sigma_dx_diag = Sigma_dx.diagonal();
      Sigma_bar_sum
          += (Sigma_dx.array().square()).matrix() + Sigma_dx_diag * Sigma_dx_diag.transpose();
    }
    Q_hat /= N_steps;

    LUPNT_CHECK((Sigma_bar_sum.diagonal().array() >= 0).all(),
                "Sigma_bar_sum has negative diagonal", "AdaptiveModel");
    return std::make_tuple(Q_hat, Sigma_bar_sum);
  }

  ArrXd ProcessNoise::SolveForQtilde(const MatXd& Q_hat, const MatXd& Sigma_bar_sum,
                                     const MatXd& C_coeffs) {
    // Check dimensions
    int n = 3;
    LUPNT_CHECK(Q_hat.rows() == 2 * n && Q_hat.cols() == 2 * n, "Q_hat has incorrect dimensions",
                "AdaptiveModel");
    LUPNT_CHECK(Sigma_bar_sum.rows() == 2 * n && Sigma_bar_sum.cols() == 2 * n,
                "Sigma_bar_sum has incorrect dimensions", "AdaptiveModel");

    MatXd Q_a = MatXd::Zero(n, n);

    for (int i = 0; i < n; ++i) {
      Arr3d Xbar_i;  // [C_11, C_21, C_22]
      if (algorithm_ == ProcessNoiseAlgorithm::ASNC) {
        Xbar_i = Arr3d{C_coeffs(0), C_coeffs(1), C_coeffs(2)};
      } else {
        Xbar_i = Arr3d{C_coeffs(0, i), C_coeffs(1, i), C_coeffs(2, i)};
      }
      Arr3d b_i{Q_hat(i, i), Q_hat(i + n, i), Q_hat(i + n, i + n)};
      Arr3d W_i_diag{Sigma_bar_sum(i, i), Sigma_bar_sum(i + n, i), Sigma_bar_sum(i + n, i + n)};

      Arr3d W_i_inv_diag = 1.0 / W_i_diag;

      double denom = (Xbar_i * W_i_inv_diag * Xbar_i).sum();
      double num = (Xbar_i * W_i_inv_diag * b_i).sum();
      if (denom < EPS) denom = EPS;

      // Q_a(i, i) = std::clamp(num / denom, Q_diag_min_, Q_diag_max_);
      Q_a(i, i) = num / denom;
      if (Q_a(i, i) < Q_diag_min_) Q_a(i, i) = Q_diag_min_;
      if (Q_a(i, i) > Q_diag_max_) Q_a(i, i) = Q_diag_max_;
    }
    LUPNT_CHECK((Q_a.array() >= 0).all(), "Q_a has negative diagonal", "AdaptiveModel");

    // std::cout << std::scientific << std::setprecision(3) << "Q_a: \n" << Q_a << std::endl;
    // std::cout << "Q_hat: \n" << Q_hat << std::endl;
    // std::cout << "Sigma_bar_sum: \n" << Sigma_bar_sum << std::endl;
    // std::cout << "C_coeffs: \n" << C_coeffs << std::endl;
    // std::cout << "Q_min: " << Q_diag_min_ << ", Q_max: " << Q_diag_max_ << std::endl;
    return Q_a;
  }

  MatXd ProcessNoise::ComputeAccNoise() {
    LUPNT_CHECK(dt_ > 0.0, "dt must be set", "Asnc");
    LUPNT_CHECK(Q_a_.size() > 0, "Qa_initial must be set", "Asnc");
    if (algorithm_ == ProcessNoiseAlgorithm::SNC) {
      return Q_a_;
    }

    if (dx_list_.size() < N_window_ || step_ < N_wait_) {
      return Q_a_;
    }

    auto [Q_hat, Sigma_bar_sum] = ComputeSlidingWindowStats();

    MatXd C_coeffs;
    if (algorithm_ == ProcessNoiseAlgorithm::ASNC) {
      C_coeffs = ProcessNoisePosVelCoeffs(dt_);
    } else {
      C_coeffs = ProcessNoisePosVelAccCoeffs(dt_, beta_);
    }

    MatXd Q_a = SolveForQtilde(Q_hat, Sigma_bar_sum, C_coeffs);

    if (alpha_ < 1.0) {
      Q_a_ = (1.0 - alpha_) * Q_a_ + alpha_ * Q_a;
    } else {
      Q_a_ = Q_a;
    }
    return Q_a_;
  }

  MatXd ProcessNoise::ComputeProcessNoise() {
    MatXd Q_a = ComputeAccNoise();

    if (algorithm_ == ProcessNoiseAlgorithm::ASNC || algorithm_ == ProcessNoiseAlgorithm::SNC) {
      return ProcessNoisePosVel(Q_a_, dt_);
    } else {
      return ProcessNoisePosVelAcc(Q_a_, dt_, beta_);
    }
  }

  ProcessNoiseFunction ProcessNoise::GetProcessNoiseFunction() {
    return [this](const State& x, Real t0, Real tf) {
      (void)x;
      (void)t0;
      (void)tf;
      return ComputeProcessNoise();
    };
  }

}  // namespace lupnt

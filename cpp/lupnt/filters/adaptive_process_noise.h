#pragma once
#include <deque>

#include "lupnt/core/definitions.h"
#include "lupnt/filters/filter.h"

namespace lupnt {

  enum class ProcessNoiseAlgorithm { SNC, ASNC, ADMC };

  class ProcessNoise {
  public:
    ProcessNoise() = default;

    void Update(const VecXd& dx, const MatXd& Sigma_dx, const MatXd& P_bar, const MatXd& P_post);

    MatXd ComputeAccNoise();
    MatXd ComputeProcessNoise();

    void SetAlgorithm(ProcessNoiseAlgorithm algorithm) { algorithm_ = algorithm; }
    void SetWaitSteps(int N_wait) { N_wait_ = N_wait; }
    void SetWindowSize(int N_window) { N_window_ = N_window; }
    void SetNoiseLimits(double sigma_min, double sigma_max) {
      Q_diag_min_ = sigma_min * sigma_min;
      Q_diag_max_ = sigma_max * sigma_max;
    }
    void SetTimeStep(double dt) { dt_ = dt; }

    std::tuple<MatXd, MatXd> ComputeSlidingWindowStats();
    ArrXd SolveForQtilde(const MatXd& Q_hat, const MatXd& Sigma_bar, const MatXd& C_coeffs);

    void SetProcessNoise(const MatXd& Q_a) { Q_a_ = Q_a; }
    void SetBeta(const VecXd& beta) { beta_ = beta; }
    void SetAlpha(double alpha) { alpha_ = alpha; }
    void Reset();
    ProcessNoiseFunction GetProcessNoiseFunction();

  protected:
    int step_ = 0;
    int N_wait_ = 0;
    int N_window_ = 30;
    double Q_diag_min_ = 0.0;
    double Q_diag_max_ = 1.0;
    double dt_ = -1.0;
    MatXd Q_a_;
    VecXd beta_;

    ProcessNoiseAlgorithm algorithm_ = ProcessNoiseAlgorithm::ASNC;

    double alpha_ = 1.0;
    MatXd Q_tilde_diag_;

    std::deque<VecXd> dx_list_;
    std::deque<MatXd> P_post_list_;
    std::deque<MatXd> P_bar_list_;
    std::deque<MatXd> Sigma_dx_list_;
  };

}  // namespace lupnt

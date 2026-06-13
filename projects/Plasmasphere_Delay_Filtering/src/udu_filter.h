#pragma once

#include <lupnt/lupnt.h>

namespace filtering_sim {
  using namespace lupnt;

  class UDUEKF : public EKF {
  protected:
    MatXd U_minus_;  // Time-updated U matrix (2N)
    MatXd U_plus_;   // Measurement-updated U matrix (2N)
    MatXd U_;        // Current U matrix. (2N after time update, N after measurement update)

    VecXd D_minus_diag_;  // Time-updated D diagonal (2N)
    VecXd D_plus_diag_;   // Measurement-updated D diagonal (2N)
    VecXd D_diag_;        // Current D diagonal (2N after time update, N after measurement update)

    // Smoother Logging
    std::vector<MatXd> U_sm_;
    std::vector<VecXd> D_sm_;
    std::vector<VecXd> x_prior_log_;  // Logged prior states
    std::vector<VecXd> x_pos_log_;    // Logged posterior states

  public:
    UDUEKF() : EKF() {}
    UDUEKF(bool is_augmented) : EKF() {}
    ~UDUEKF() override = default;

    void SetUd();

    void SetCovariance(const MatXd& P) override;

    void CarlsonUpdate();

    void Predict(Real t, const State* u = nullptr) override;
    void Update(const VecX& z_true) override;

    // Smoothers
    void InitializeLogger(int max_tidx) override;
    void InitializeSmootherState() override;
    void LogFilterEstimate(int tidx) override;
    void UpdateSmoother(int tidx) override;
  };

  class UDUDelayedEKF : public UDUEKF {
  protected:
    // Submatrices
    MatXd U11_;
    MatXd U12_;
    MatXd U22_;
    VecXd D11_diag_;
    VecXd D22_diag_;
    MatXd P11_;
    MatXd P12_;
    MatXd P22_;

    // Smoother logging
    std::vector<MatXd> U12_plus_log_;  // Logged prior covariances
    std::vector<VecXd> D22_plus_log_;  // Logged posterior states
    std::vector<MatXd> U22_plus_log_;  // Logged posterior covariances
    std::vector<MatXd> P22_log_;       // Logged prior covariance
    std::vector<MatXd> P12_log_;       // Logged prior covariance
    std::vector<MatXd> P11_log_;       // Logged prior covariance
    std::vector<MatXd> Ukk_log_;       // Logged posterior covariance
    std::vector<VecXd> Dkk_log_;       // Logged posterior covariance

  public:
    UDUDelayedEKF() : UDUEKF() {}
    ~UDUDelayedEKF() override = default;

    void SetUd();

    void SetCovariance(const MatXd& P) override;

    void StoreSubmatrices(int N_state);

    /**
     * Update the U
     */
    void UpdateCurrStepUD();

    void DelayedCarlsonUpdate();

    void Predict(Real t, const State* u = nullptr) override;
    void Update(const VecX& z_true) override;

    // Smoothers
    void InitializeLogger(int max_tidx) override;
    void InitializeSmootherState() override;
    void LogFilterEstimate(int tidx) override;
    void UpdateSmoother(int tidx) override;
  };

}  // namespace filtering_sim

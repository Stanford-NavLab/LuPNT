#pragma once

#include <lupnt/lupnt.h>

namespace filtering_sim {
  using namespace lupnt;

  class DelayedEKF : public EKF {
  private:
    VecXd x_aug_;  // Augmented state for delayed measurements
    MatXd P_aug_;  // Augmented state and covariance for delayed measurements

  public:
    DelayedEKF() : EKF() {}
    ~DelayedEKF() override = default;

    void Update(const VecX& z_true) override;

    void Predict(Real t, const State* u = nullptr) override;

    VecXd GetAugmentedState() { return x_aug_; }
    MatXd GetAugmentedCovariance() { return P_aug_; }
  };

}  // namespace filtering_sim

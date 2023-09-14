#pragma once

#include <lupnt/numerics/eigenmvn.h>

#include <Eigen/Dense>
#include <tuple>

#include "State.h"

namespace LPT {

enum class ClockModel {
  kMicrosemiCsac,
  kRafs,
  kUso,
  kMiniRafs,
};

static std::tuple<double, double, double> GetClockSigma(ClockModel clk_model) {
  double sigma1 = NAN;
  double sigma2 = 1e-25;
  double sigma3 = 1e-40;

  switch (clk_model) {
    case ClockModel::kMicrosemiCsac:
      // file:///Users/keidaiiiyama/Downloads/Microchip_CSAC_Space_Datasheet_900-00744-007C%20(1).pdf
      sigma1 = 3e-10;
      sigma2 = 3e-12;
      sigma3 = 0;
      break;
    case ClockModel::kRafs:
      // 'Lunar Relay Onboard Navigation Performance and Effects on Lander
      // 'Decent to Surface
      sigma1 = 1.923538e-12;
      sigma2 = 4.324350e-17;
      sigma3 = 8.694826e-30;
      break;
    case ClockModel::kUso:
      // 'Lunar Relay Onboard Navigation Performance and Effects on Lander
      // 'Decent to Surface'
      sigma1 = sqrt(2.53e-23);
      sigma2 = sqrt(4.22e-24);
      sigma3 = sqrt(1.00e-38);
      break;
    case ClockModel::kMiniRafs:
      sigma1 = 1.0e-11;
      sigma2 = 1.1e-15;
      break;
    default:
      throw std::runtime_error("Unknown clock model");
      break;
  }
  return std::make_tuple(sigma1, sigma2, sigma3);
}

static Eigen::Matrix2d GetClockProcessNoise(ClockModel clk_model, double dt) {
  double sig1, sig2, sig3;
  std::tie(sig1, sig2, sig3) = GetClockSigma(clk_model);
  Eigen::Matrix2d Q;
  Q(0, 0) = std::pow(sig1, 2) * dt + std::pow(sig2, 2) * std::pow(dt, 3) / 3;
  Q(0, 1) = std::pow(sig2, 2) * std::pow(dt, 2) / 2;
  Q(1, 0) = std::pow(sig2, 2) * std::pow(dt, 2) / 2;
  Q(1, 1) = std::pow(sig2, 2) * dt;
  return Q;
};

/**
 * @brief Clock State Class
 *
 */
class ClockState : public IState {
 private:
  ad::VectorXreal x_;
  int state_size_;  // it can be 2 or 3

 public:
  ClockState(int state_size) {
    state_size_ = state_size;
    x_ = ad::VectorXreal::Zero(state_size_);
  }
  ClockState(ad::VectorXreal clock_vec) {
    x_ = clock_vec;
    state_size_ = clock_vec.size();
  }
  ad::VectorXreal GetVector() const { return x_; }
  ad::real GetValue(int i) const { return x_(i); }
  int GetStateSize() { return state_size_; };
  inline void SetValue(const ad::real val, const int idx) { x_(idx) = val; }
  void SetVector(const ad::VectorXreal& x) { x_ = x; }
};

class ClockDynamics {
 private:
  ClockModel clk_model_;
  std::unique_ptr<Eigen::EigenMultivariateNormal<double>> mvrnd_clk_;
  double dt_;

 public:
  ClockDynamics(ClockModel clk_model) {
    clk_model_ = clk_model;
    mvrnd_clk_ = std::make_unique<Eigen::EigenMultivariateNormal<double>>(
        Eigen::Vector2d::Zero(), Eigen::Matrix2d::Identity());
    dt_ = 0;
  }

  ClockDynamics(const ClockDynamics& other) {
    clk_model_ = other.clk_model_;
    mvrnd_clk_ = std::make_unique<Eigen::EigenMultivariateNormal<double>>(
        Eigen::Vector2d::Zero(), Eigen::Matrix2d::Identity());
    dt_ = 0;
  }

  void Propagate(ad::Vector2real& clk, ad::real dt) {
    if (dt.val() != dt_) {
      dt_ = dt.val();
      auto Q_clk = GetClockProcessNoise(clk_model_, dt.val());
      mvrnd_clk_->setCovar(Q_clk);
    }
    Eigen::Matrix2d Phi_clk_{{1, dt.val()}, {0, 1}};
    clk = Phi_clk_ * clk + mvrnd_clk_->samples(1);
  }
};
}  // namespace LPT
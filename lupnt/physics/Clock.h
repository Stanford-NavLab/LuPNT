/**
 * @file clock.h
 * @author Stanford NAV LAB
 * @brief Clock class
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <tuple>

#include "state.h"

namespace lupnt {

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

static Matrix2d GetClockProcessNoise(ClockModel clk_model, double dt) {
  double sig1, sig2, sig3;
  std::tie(sig1, sig2, sig3) = GetClockSigma(clk_model);
  Matrix2d Q;
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
  VectorX x_;
  int state_size_;  // it can be 2 or 3

 public:
  ClockState(int state_size) {
    state_size_ = state_size;
    x_ = VectorX::Zero(state_size_);
  }
  ClockState(VectorX clock_vec) {
    x_ = clock_vec;
    state_size_ = clock_vec.size();
  }
  VectorX GetVector() const { return x_; }
  real GetValue(int i) const { return x_(i); }
  int GetSize() { return state_size_; };
  inline void SetValue(const real val, const int idx) { x_(idx) = val; }
  void SetVector(const VectorX& x) { x_ = x; }
};

class ClockDynamics {
 private:
  ClockModel clk_model_;
  double dt_;

 public:
  ClockDynamics(ClockModel clk_model) {
    clk_model_ = clk_model;
    dt_ = 0;
  }

  ClockDynamics(const ClockDynamics& other) {
    clk_model_ = other.clk_model_;
    dt_ = 0;
  }

  void Propagate(Vector2& clk, real dt) {
    if (dt.val() != dt_) {
      dt_ = dt.val();
      auto Q_clk = GetClockProcessNoise(clk_model_, dt.val());
      Matrix2d Phi_clk_{{1, dt.val()}, {0, 1}};
      clk = Phi_clk_ * clk + SampleMVN(Vector2d::Zero(), Q_clk, 1);
    }
  }
};
}  // namespace lupnt
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

#include "lupnt/dynamics/dynamics.h"
#include "lupnt/numerics/math_utils.h"
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

  static Mat2d GetClockProcessNoise(ClockModel clk_model, double dt) {
    double sig1, sig2, sig3;
    std::tie(sig1, sig2, sig3) = GetClockSigma(clk_model);
    Mat2d Q;
    Q(0, 0) = std::pow(sig1, 2) * dt + std::pow(sig2, 2) * std::pow(dt, 3) / 3;
    Q(0, 1) = std::pow(sig2, 2) * std::pow(dt, 2) / 2;
    Q(1, 0) = std::pow(sig2, 2) * std::pow(dt, 2) / 2;
    Q(1, 1) = std::pow(sig2, 2) * dt;
    return Q;
  };

  static Mat3d GetClockProcessNoise3(ClockModel clk_model, double dt) {
    double sig1, sig2, sig3, t2, t3, t4, t5;
    std::tie(sig1, sig2, sig3) = GetClockSigma(clk_model);

    t2 = std::pow(dt, 2);
    t3 = std::pow(dt, 3);
    t4 = std::pow(dt, 4);
    t5 = std::pow(dt, 5);

    sig1 = std::pow(sig1, 2);
    sig2 = std::pow(sig2, 2);
    sig3 = std::pow(sig3, 2);

    Mat3d Q;
    Q(0, 0) = sig1 * dt + std::pow(sig2, 2) * t3 / 3 + sig3 * t5 / 20;
    Q(0, 1) = sig2 * t2 / 2 + sig3 * t4 / 8;
    Q(0, 2) = sig3 * t3 / 6;
    Q(1, 0) = Q(0, 1);
    Q(1, 1) = sig2 * dt + sig3 * t3 / 3;
    Q(1, 2) = sig3 * t2 / 2;
    Q(2, 0) = Q(0, 2);
    Q(2, 1) = Q(1, 2);
    Q(2, 2) = sig3 * dt;
    return Q;
  };

  /**
   * @brief Clock State Class
   *
   */
  class ClockState : public IState {
  private:
    VecX x_;
    int state_size_;  // it can be 2 or 3

  public:
    ClockState(int state_size) {
      state_size_ = state_size;
      x_ = VecX::Zero(state_size_);
    }
    ClockState(VecX clock_vec) {
      x_ = clock_vec;
      state_size_ = clock_vec.size();
    }
    VecX GetVec() const { return x_; }
    Real GetValue(int i) const { return x_(i); }
    inline int GetSize() const { return state_size_; };
    inline void SetValue(const Real val, const int idx) { x_(idx) = val; }
    void SetVec(const VecX& x) { x_ = x; }
  };

  class ClockDynamics : public IDynamics {
  private:
    ClockModel clk_model_;
    double dt_;

  public:
    ClockDynamics(ClockModel clk_model) { clk_model_ = clk_model; }

    ClockDynamics(const ClockDynamics& other) { clk_model_ = other.clk_model_; }

    Mat2 TwoStatePhi(Real dt) {
      Mat2 Phi_clk;
      Phi_clk << 1, dt, 0, 1;
      return Phi_clk;
    }

    Mat2d TwoStatePhi(double dt) {
      Mat2d Phi_clk;
      Phi_clk << 1, dt, 0, 1;
      return Phi_clk;
    }

    Mat3 ThreeStatePhi(Real dt) {
      Mat3 Phi_clk{{1, dt, dt * dt / 2}, {0, 1, dt}, {0, 0, 1}};
      return Phi_clk;
    }

    Mat3d ThreeStatePhi(double dt) {
      Mat3d Phi_clk{{1, dt, dt * dt / 2}, {0, 1, dt}, {0, 0, 1}};
      return Phi_clk;
    }

    // two state clock
    void Propagate(Vec2& clk, Real t0, Real tf) {
      Real dt = tf - t0;
      Mat2 Phi_clk = TwoStatePhi(dt);
      clk = Phi_clk * clk;
    }

    void Propagate(Vec2& clk, Real dt) {
      Mat2 Phi_clk = TwoStatePhi(dt);
      clk = Phi_clk * clk;
    }

    void PropagateWithNoise(Vec2& clk, Real t0, Real tf) {
      Real dt = tf - t0;
      auto Q_clk = GetClockProcessNoise(clk_model_, dt.val());
      Mat2 Phi_clk_ = TwoStatePhi(dt);
      clk = SampleMVN(Phi_clk_ * clk, Q_clk, 1);
    }

    void PropagateWithStm(Vec2& clk, Real t0, Real tf, MatXd& stm) {
      Real dt = tf - t0;
      Mat2 Phi_clk = TwoStatePhi(dt);
      Mat2d Phi_clk_d = Phi_clk.cast<double>();
      clk = Phi_clk * clk;
      stm.block(0, 0, 2, 2) = Phi_clk_d;
    }

    // three state clock
    void Propagate(Vec3& clk, Real t0, Real tf) {
      Real dt = tf - t0;
      Mat3 Phi_clk = ThreeStatePhi(dt);
      clk = Phi_clk * clk;
    }

    void Propagate(Vec3& clk, Real dt) {
      Mat3 Phi_clk = ThreeStatePhi(dt);
      clk = Phi_clk * clk;
    }

    void PropagateWithNoise(Vec3& clk, Real t0, Real tf) {
      Real dt = tf - t0;
      auto Q_clk = GetClockProcessNoise3(clk_model_, dt.val());
      Mat3 Phi_clk_ = ThreeStatePhi(dt);
      clk = Phi_clk_ * clk + SampleMVN(Vec3d::Zero(), Q_clk, 1);
    }

    void PropagateWithStm(Vec3& clk, Real t0, Real tf, MatXd& stm) {
      Real dt = tf - t0;
      Mat3 Phi_clk = ThreeStatePhi(dt);
      Mat3d Phi_clk_d = Phi_clk.cast<double>();
      clk = Phi_clk * clk;
      stm.block(0, 0, 3, 3) = Phi_clk_d;
    }

    // arbitrary state clock
    void PropagateX(VecX& clk, Real t0, Real tf) {
      if (clk.size() == 2) {
        Vec2 clk2 = clk;
        Propagate(clk2, t0, tf);
        clk = clk2;
      } else if (clk.size() == 3) {
        Vec3 clk3 = clk;
        Propagate(clk3, t0, tf);
        clk = clk3;
      } else {
        assert(false && "Invalid clock state size");
      }
    }

    void PropagateWithNoiseX(VecX& clk, Real t0, Real tf) {
      if (clk.size() == 2) {
        Vec2 clk2 = clk;
        PropagateWithNoise(clk2, t0, tf);
        clk = clk2;
      } else if (clk.size() == 3) {
        Vec3 clk3 = clk;
        PropagateWithNoise(clk3, t0, tf);
        clk = clk3;
      } else {
        assert(false && "Invalid clock state size");
      }
    }

    void PropagateWithStmX(VecX& clk, Real t0, Real tf, MatXd& stm) {
      if (clk.size() == 2) {
        Vec2 clk2 = clk;
        stm.resize(2, 2);
        PropagateWithStm(clk2, t0, tf, stm);
        clk = clk2;
      } else if (clk.size() == 3) {
        Vec3 clk3 = clk;
        stm.resize(3, 3);
        PropagateWithStm(clk3, t0, tf, stm);
        clk = clk3;
      } else {
        assert(false && "Invalid clock state size");
      }
    }

    // Propagate using clockState
    void Propagate(ClockState& clk, Real t0, Real tf) {
      VecX clk_vec = clk.GetVec();
      PropagateX(clk_vec, t0, tf);
      clk.SetVec(clk_vec);
    }

    void PropagateWithNoise(ClockState& clk, Real t0, Real tf) {
      VecX clk_vec = clk.GetVec();
      PropagateWithNoiseX(clk_vec, t0, tf);
      clk.SetVec(clk_vec);
    }
  };
}  // namespace lupnt
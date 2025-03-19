#include "lupnt/physics/clock.h"

#include <tuple>

#include "lupnt/core/progress_bar.h"

namespace lupnt {

  std::tuple<Real, Real, Real> GetClockSigma(ClockModel clk_model) {
    Real sigma1 = NAN;
    Real sigma2 = 1e-25;
    Real sigma3 = 1e-40;

    switch (clk_model) {
      case ClockModel::kMicrosemiCsac:
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
      default: throw std::runtime_error("Unknown clock model"); break;
    }
    return std::make_tuple(sigma1, sigma2, sigma3);
  }

  Mat2 ClockDynamics::TwoStateNoise(ClockModel clk_model, Real dt) {
    Real sig1, sig2, sig3;
    std::tie(sig1, sig2, sig3) = GetClockSigma(clk_model);
    Mat2 Q;
    Q(0, 0) = pow(sig1, 2) * dt + pow(sig2, 2) * pow(dt, 3) / 3;
    Q(0, 1) = pow(sig2, 2) * pow(dt, 2) / 2;
    Q(1, 0) = pow(sig2, 2) * pow(dt, 2) / 2;
    Q(1, 1) = pow(sig2, 2) * dt;
    return Q;
  };

  Mat3 ClockDynamics::ThreeStateNoise(ClockModel clk_model, Real dt) {
    Real sig1, sig2, sig3;
    std::tie(sig1, sig2, sig3) = GetClockSigma(clk_model);

    Real t2 = pow(dt, 2);
    Real t3 = pow(dt, 3);
    Real t4 = pow(dt, 4);
    Real t5 = pow(dt, 5);

    sig1 = pow(sig1, 2);
    sig2 = pow(sig2, 2);
    sig3 = pow(sig3, 2);

    Mat3 Q;
    Q(0, 0) = sig1 * dt + pow(sig2, 2) * t3 / 3 + sig3 * t5 / 20;
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

  ClockState::ClockState() = default;
  ClockState::ClockState(int state_size) {
    state_size_ = state_size;
    x_ = VecX::Zero(state_size_);
  }
  ClockState::ClockState(VecX clock_vec) {
    x_.resize(clock_vec.size());
    x_ = clock_vec;
    state_size_ = clock_vec.size();
  }

  // Overrides
  int ClockState::GetSize() const { return state_size_; };
  VecX ClockState::GetVec() const { return x_; }
  void ClockState::SetVec(const VecX& x) { x_ = x; }
  Real ClockState::GetValue(int i) const { return x_(i); }
  void ClockState::SetValue(int idx, Real val) { x_(idx) = val; }
  StateType ClockState::GetStateType() const { return static_cast<StateType>(-1); }

  ClockDynamics::ClockDynamics() = default;
  ClockDynamics::ClockDynamics(ClockModel clk_model) { clk_model_ = clk_model; }

  Mat2 ClockDynamics::TwoStatePhi(Real dt) {
    Mat2 Phi_clk{{1, dt}, {0, 1}};
    return Phi_clk;
  }

  Mat3 ClockDynamics::ThreeStatePhi(Real dt) {
    Mat3 Phi_clk{{1, dt, dt * dt / 2}, {0, 1, dt}, {0, 0, 1}};
    return Phi_clk;
  }

  ClockState ClockDynamics::PropagateState(const ClockState state, Real t0, Real tf, MatXd* stm) {
    VecX x0 = state.GetVec();
    VecX xf = Propagate(x0, t0, tf, stm);
    return ClockState(xf);
  }

  Ptr<IState> ClockDynamics::PropagateState(const Ptr<IState>& state, Real t0, Real tf,
                                            MatXd* stm) {
    ClockState* clk_state = dynamic_cast<ClockState*>(state.get());
    if (clk_state == nullptr) {
      throw std::runtime_error("Invalid state type");
    }
    return std::make_shared<ClockState>(PropagateState(*clk_state, t0, tf, stm));
  }

  Vec2 ClockDynamics::Propagate(const Vec2& x0, Real t0, Real tf, Mat2* stm) {
    Mat2 Phi_clk = TwoStatePhi(tf - t0);
    Vec2 clkf = Phi_clk * x0;
    if (noise_) {
      Mat2 Q_clk = TwoStateNoise(clk_model_, tf - t0);
      Vec2 noise = SampleMVN(Vec2::Zero(), Q_clk, 1);
      clkf += noise;
    }
    if (stm != nullptr) *stm = Phi_clk;
    return clkf;
  }

  Vec3 ClockDynamics::Propagate(const Vec3& x0, Real t0, Real tf, Mat3* stm) {
    Mat3 Phi_clk = ThreeStatePhi(tf - t0);
    Vec3 clkf = Phi_clk * x0;
    if (noise_) {
      Mat3 Q_clk = ThreeStateNoise(clk_model_, tf - t0);
      Vec3 noise = SampleMVN(Vec3::Zero(), Q_clk, 1);
      clkf += noise;
    }
    if (stm != nullptr) *stm = Phi_clk;
    return clkf;
  }

  VecX ClockDynamics::Propagate(const VecX& x0, Real t0, Real tf, MatXd* stm) {
    int state_size = x0.size();
    if (state_size == 2) {
      Vec2 clk0 = ClockState(x0).GetVec();
      Mat2 stm_dum;
      Vec2 clkf = Propagate(clk0, t0, tf, &stm_dum);
      if (stm != nullptr) {
        stm->resize(2, 2);
        *stm = stm_dum.cast<double>();
      }
      return clkf;
    } else if (state_size == 3) {
      Vec3 clk0 = ClockState(x0).GetVec();
      Mat3 stm_dum;
      Vec3 clkf = Propagate(clk0, t0, tf, &stm_dum);
      if (stm != nullptr) {
        stm->resize(3, 3);
        *stm = stm_dum.cast<double>();
      }
      return clkf;
    } else {
      throw std::runtime_error("Invalid clock state size");
    }
  }

}  // namespace lupnt

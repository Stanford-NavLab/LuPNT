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

  enum class ClockModel { kMicrosemiCsac, kRafs, kUso, kMiniRafs, kUnknown };

  std::tuple<Real, Real, Real> GetClockSigma(ClockModel clk_model);

  class ClockState : public IState {
  private:
    VecX x_;
    int state_size_;  // it can be 2 or 3

  public:
    ClockState();
    ClockState(int state_size);
    ClockState(VecX clock_vec);

    // Overrides
    int GetSize() const;
    VecX GetVec() const;
    void SetVec(const VecX& x);
    Real GetValue(int i) const;
    void SetValue(int idx, Real val);
    StateType GetStateType() const;
  };

  class ClockDynamics : public IDynamics {
  private:
    ClockModel clk_model_ = ClockModel::kUnknown;
    bool noise_ = false;

  public:
    ClockDynamics();
    ClockDynamics(ClockModel clk_model);

    inline void SetNoise(bool noise) { noise_ = noise; }

    static Mat2 TwoStatePhi(Real dt);
    static Mat3 ThreeStatePhi(Real dt);

    static Mat2 TwoStateNoise(ClockModel clk_model, Real dt);
    static Mat3 ThreeStateNoise(ClockModel clk_model, Real dt);

    using IDynamics::Propagate;
    ClockState PropagateState(const ClockState state, Real t0, Real tf, MatXd* stm = nullptr);
    Ptr<IState> PropagateState(const Ptr<IState>& state, Real t0, Real tf,
                               MatXd* stm = nullptr) override;
    VecX Propagate(const VecX& x0, Real t0, Real tf, MatXd* stm = nullptr) override;
    Vec2 Propagate(const Vec2& x0, Real t0, Real tf, Mat2* stm = nullptr);
    Vec3 Propagate(const Vec3& x0, Real t0, Real tf, Mat3* stm = nullptr);
  };
}  // namespace lupnt

/**
 * @file agent.cpp
 * @author Stanford NAV Lab
 * @brief Base class fo
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#include "lupnt/agents/agent.h"

#include "lupnt/physics/body.h"

namespace lupnt {

  int Agent::id_counter_ = 0;

  /********************** Agent  ********************************/

  void Agent::Propagate(const Real epoch) {
    if (epoch == epoch_) return;

    // Update orbit state
    VecX x0 = rv_->GetVec();
    VecX xf = dynamics_->Propagate(x0, epoch_, epoch);
    rv_->SetVec(xf);

    // Update clock state
    if (clock_dynamics_ != nullptr) {
      ClockState clkf = GetClockStateAtEpoch(epoch, true);
      VecX clk0_vec = clock_.GetVec();
      VecX clkf_vec = clkf.GetVec();
      VecX dclk = clkf_vec - clk0_vec;
      clock_ = clkf;
    }

    // TBD: Update attitude state
    epoch_ = epoch;
  }

  VecX Agent::GetRvStateAtEpoch(const Real epoch) {
    if (epoch == epoch_) return rv_->GetVec();

    // Get the state at epoch without changing the agent's epoch and state
    VecX x = rv_->GetVec();
    VecX xf = dynamics_->Propagate(x, epoch_, epoch);

    return xf;
  }

  ClockState Agent::GetClockStateAtEpoch(const Real epoch, bool with_noise) {
    if (epoch == epoch_) return clock_;

    // Get the state at epoch without changing the agent's epoch and clock state
    // Create a deep copy of the clock state
    ClockState clk = clock_;
    ClockState clk_new = ClockState(clk.GetVec());
    MatXd stm;
    if (clock_dynamics_ != nullptr && with_noise) {
      clock_dynamics_->SetNoise(true);
      clk_new = clock_dynamics_->PropagateState(clk, epoch_, epoch, &stm);
    } else {
      clock_dynamics_->SetNoise(false);
      clk_new = clock_dynamics_->PropagateState(clk, epoch_, epoch, &stm);
    }
    return clk_new;
  }

  VecX Agent::GetClockStateVecAtEpoch(const Real epoch, bool with_noise) {
    if (epoch == epoch_) return clock_.GetVec();

    // Get the state at epoch without changing the agent's epoch and clock state
    // Create a deep copy of the clock state
    ClockState clk = GetClockStateAtEpoch(epoch, with_noise);
    return clk.GetVec();
  }

  /********************** SpaceCraft  ********************************/

  CartesianOrbitState Spacecraft::GetCartesianGCRFStateAtEpoch(Real epoch) {
    Ptr<OrbitState> state = GetOrbitState();

    Real current_epoch = GetEpoch();
    Ptr<NumericalOrbitDynamics> dynamics
        = std::dynamic_pointer_cast<NumericalOrbitDynamics>(GetDynamics());

    if (epoch != current_epoch) {
      // set dt
      dynamics->PropagateState(*state, current_epoch, epoch);
    }
    // TODO
    Real GM = GetBodyData(GetBodyId()).GM;
    Ptr<CartesianOrbitState> cartOrbitState = std::static_pointer_cast<CartesianOrbitState>(
        ConvertOrbitStateRepresentation(state, OrbitStateRepres::CARTESIAN, GM));
    return ConvertOrbitStateFrame(*cartOrbitState, epoch, Frame::GCRF);
  }

};  // namespace lupnt

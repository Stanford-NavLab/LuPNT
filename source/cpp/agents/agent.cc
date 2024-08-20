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

namespace lupnt {

  int Agent::id_counter_ = 0;

  /********************** Agent  ********************************/

  void Agent::Propagate(const Real epoch) {
    if (epoch == epoch_) return;

    // Update orbit state
    VecX x = rv_->GetVecX();
    dynamics_->PropagateX(x, epoch_, epoch);
    rv_->SetVecX(x);

    // Update clock state
    if (clock_dynamics_ != nullptr) {
      clock_dynamics_->PropagateWithNoise(clock_, epoch_, epoch);
    }

    // TBD: Update attitude state

    epoch_ = epoch;
  }

  VecX Agent::GetRvStateAtEpoch(const Real epoch) {
    if (epoch == epoch_) return rv_->GetVecX();

    // Get the state at epoch without changing the agent's epoch and state
    VecX x = rv_->GetVecX();
    dynamics_->PropagateX(x, epoch_, epoch);

    return x;
  }

  ClockState Agent::GetClockStateAtEpoch(const Real epoch, bool with_noise) {
    if (epoch == epoch_) return clock_;

    // Get the state at epoch without changing the agent's epoch and clock state
    // Create a deep copy of the clock state
    ClockState clk = clock_;
    if (clock_dynamics_ != nullptr && with_noise) {
      clock_dynamics_->PropagateWithNoise(clk, epoch_, epoch);
    } else {
      clock_dynamics_->Propagate(clk, epoch_, epoch);
    }
    return clk;
  }

  VecX Agent::GetClockStateVecAtEpoch(const Real epoch, bool with_noise) {
    if (epoch == epoch_) return clock_.GetVecX();

    // Get the state at epoch without changing the agent's epoch and clock state
    // Create a deep copy of the clock state
    ClockState clk = clock_;
    if (clock_dynamics_ != nullptr && with_noise) {
      clock_dynamics_->PropagateWithNoise(clk, epoch_, epoch);
    } else {
      clock_dynamics_->Propagate(clk, epoch_, epoch);
    }
    return clk.GetVecX();
  }

  /********************** SpaceCraft  ********************************/

  std::shared_ptr<CartesianOrbitState> Spacecraft::GetCartesianGCRFStateAtEpoch(Real epoch) {
    std::shared_ptr<OrbitState> state = GetOrbitState();

    Real current_epoch = GetEpoch();
    std::shared_ptr<NumericalOrbitDynamics> dynamics
        = std::dynamic_pointer_cast<NumericalOrbitDynamics>(GetDynamics());

    if (epoch != current_epoch) {
      // set dt
      Real dt = (epoch - current_epoch) / 10;
      dynamics->Propagate(*state, current_epoch, epoch, dt);
    }
    // TODO
    Real GM = 0.0;  //  GetBodyData(bodyId_).GM;
    auto cartOrbitState = std::static_pointer_cast<CartesianOrbitState>(
        ConvertOrbitStateRepresentation(state, OrbitStateRepres::CARTESIAN, GM));
    return ConvertOrbitStateFrame(cartOrbitState, epoch, Frame::GCRF);
  }

};  // namespace lupnt

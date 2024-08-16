#pragma once

// lupnt includes
#include "comm_device.h"
#include "lupnt/core/constants.h"
#include "lupnt/dynamics/dynamics.h"
#include "lupnt/physics/clock.h"
#include "lupnt/physics/frame_converter.h"
#include "lupnt/physics/orbit_state.h"

namespace lupnt {

  class ICommDevice;

  class Agent {
  private:
    static int id_counter_;
    const int id_;
    std::string name_;

    NaifId bodyId_;
    Real epoch_;
    Ptr<IState> state_;
    Ptr<IDynamics> dynamics_;
    std::vector<Ptr<ICommDevice>> devices_;

    ClockState clock_;
    ClockDynamics clock_dynamics_;

  public:
    Agent() : id_(id_counter_++) {}
    int GetId() { return id_; }
    std::string GetName() { return name_; }
    void SetName(std::string name) { name_ = name; }
    // Getters
    Real GetEpoch() { return epoch_; }
    NaifId GetBodyId() { return bodyId_; }
    Ptr<IState> GetState() { return state_; }
    Ptr<IDynamics> GetDynamics() { return dynamics_; }
    ClockState GetClockState() { return clock_; }
    VecX GetStateVec() {
      Vec6 rv = GetState()->GetVec();
      Vec2 clk = GetClockState().GetVec();
      VecX state(8);
      state << rv, clk;
      return state;
    }

    // Setters
    void SetState(Ptr<IState> state) { state_ = static_cast<Ptr<IState>>(state); }
    void SetDynamics(Ptr<IDynamics> dyn) { dynamics_ = dyn; }
    void SetEpoch(Real epoch) { epoch_ = epoch; }
    void SetBodyId(NaifId bodyId) { bodyId_ = bodyId; }
    void SetClock(ClockState clk) { clock_ = clk; }
    void SetClockDynamics(ClockDynamics clock_dyn) { clock_dynamics_ = clock_dyn; }

    void AddDevice(Ptr<ICommDevice> device) { devices_.push_back(device); }

    // Cartesian State at epoch in GCRF frame
    CartesianOrbitState GetCartesianGCRFStateAtEpoch(Real epoch);

    void Propagate(const Real epoch) {
      if (epoch == epoch_) return;
      state_ = dynamics_->PropagateState(state_, epoch_, epoch);
      clock_ = clock_dynamics_.PropagateState(clock_, epoch_, epoch);
    }
  };

  class Spacecraft : public Agent {
  public:
    Spacecraft() : Agent() {}
  };

  class Rover : public Agent {
  public:
    Rover() : Agent() {}
  };

  class GroundStation : public Agent {
  public:
    GroundStation() : Agent() {}

    void SetPosition(Vec3d pos) { pos_ = pos; }
    Vec3d GetPosition() { return pos_; }

  private:
    Vec3d pos_;
  };

};  // namespace lupnt

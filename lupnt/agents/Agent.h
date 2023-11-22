/**
 * @file agent.h
 * @author Stanford NAV LAB
 * @brief List of agents
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

// C++ includes
#include <memory>





// lupnt includes
#include "../core/constants.h"
#include "../dynamics/dynamics.h"
#include "../physics/clock.h"
#include "../physics/coord_converter.h"
#include "../physics/orbit_state.h"
#include "../physics/orbit_state_utils.h"
#include "comm_device.h"

namespace lupnt {

class ICommDevice;

/**
 * @brief Agent base class
 *
 */
class Agent {
 private:
  static int id_counter_;
  const int id_;
  std::string name_;

  BodyId bodyId_;
  real epoch_;
  std::shared_ptr<OrbitState> state_;
  std::shared_ptr<IOrbitDynamics> dynamics_;
  std::vector<std::shared_ptr<ICommDevice>> devices_;

  Vector2real clock_;
  std::unique_ptr<ClockDynamics> clock_dynamics_;

 public:
  Agent() : id_(id_counter_++) {}

  // Getters
  real GetEpoch() { return epoch_; }
  BodyId GetBodyId() { return bodyId_; }
  std::shared_ptr<OrbitState> GetOrbitState() { return state_; }
  std::shared_ptr<IOrbitDynamics> GetDynamics() { return dynamics_; }
  Vector2real GetClock() { return clock_; }

  // Setters
  void SetOrbitState(std::shared_ptr<OrbitState> state) { state_ = state; }
  void SetDynamics(std::shared_ptr<IOrbitDynamics> dyn) { dynamics_ = dyn; }
  void SetEpoch(real epoch) { epoch_ = epoch; }
  void SetBodyId(BodyId bodyId) { bodyId_ = bodyId; }
  void SetClock(Vector2real clk) { clock_ = clk; }
  void SetClockDynamics(ClockDynamics& clock_dyn) {
    clock_dynamics_ = std::make_unique<ClockDynamics>(clock_dyn);
  }

  void AddDevice(std::shared_ptr<ICommDevice> device) {
    devices_.push_back(device);
  }

  // Cartesian OrbitState at epoch in GCRF frame
  std::shared_ptr<CartesianOrbitState> GetCartesianGCRFStateAtEpoch(
      real epoch, CoordSystem coord_sys = CoordSystem::GCRF);

  void Propagate(const real epoch) {
    if (epoch == epoch_) return;

    dynamics_->Propagate(*state_, epoch_, epoch, 1.0 * SECS_PER_MINUTE);
    if (clock_dynamics_ != nullptr) {
      clock_dynamics_->Propagate(clock_, epoch - epoch_);
    }

    epoch_ = epoch;
  }
};

/**
 * @brief Spacecraft Agent
 *
 */
class Spacecraft : public Agent {
 public:
  Spacecraft() : Agent() {}
};

/**
 * @brief Rover Agent
 *
 */
class Rover : public Agent {
 public:
  Rover() : Agent() {}
};

};  // namespace lupnt

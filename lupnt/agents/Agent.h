/**
 * @file Agent.h
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

// autodiff includes
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

// lupnt includes
#include "lupnt/agents/CommDevice.h"
#include "lupnt/core/Constants.h"
#include "lupnt/dynamics/Dynamics.h"
#include "lupnt/physics/Clock.h"
#include "lupnt/physics/CoordConverter.h"
#include "lupnt/physics/OrbitState.h"
#include "lupnt/physics/OrbitStateUtils.h"

namespace LPT {

class ICommDevice;

/**
 * @brief Agent base class
 *
 */
class Agent {
 private:
  static int IDCounter;
  const int ID;
  std::string name;

  BodyId bodyId_;
  ad::real epoch_;
  std::shared_ptr<OrbitState> state_;
  std::shared_ptr<IOrbitDynamics> dynamics_;
  std::vector<std::shared_ptr<ICommDevice>> devices_;

  ad::Vector2real clk_;
  std::unique_ptr<ClockDynamics> clk_dyn_;

 public:
  Agent() : ID(IDCounter++) {}

  // Getters
  ad::real GetEpoch() { return epoch_; }
  BodyId GetBodyId() { return bodyId_; }
  std::shared_ptr<OrbitState> GetOrbitState() { return state_; }
  std::shared_ptr<IOrbitDynamics> GetDynamics() { return dynamics_; }
  ad::Vector2real GetClock() { return clk_; }

  // Setters
  void SetOrbitState(std::shared_ptr<OrbitState> state) { state_ = state; }
  void SetDynamics(std::shared_ptr<IOrbitDynamics> dyn) { dynamics_ = dyn; }
  void SetEpoch(ad::real epoch) { epoch_ = epoch; }
  void SetBodyId(BodyId bodyId) { bodyId_ = bodyId; }
  void SetClock(ad::Vector2real clk) { clk_ = clk; }
  void SetClockDynamics(ClockDynamics& clk_dyn) {
    clk_dyn_ = std::make_unique<ClockDynamics>(clk_dyn);
  }

  void AddDevice(std::shared_ptr<ICommDevice> device) {
    devices_.push_back(device);
  }

  // Cartesian OrbitState at epoch in GCRF frame
  std::shared_ptr<CartesianOrbitState> GetCartesianGCRFStateAtEpoch(
      ad::real epoch, CoordSystem coord_sys = CoordSystem::GCRF);

  void Propagate(const ad::real epoch) {
    if (epoch == epoch_) return;

    dynamics_->Propagate(*state_, epoch_, epoch, 1.0 * SECS_PER_MINUTE);
    if (clk_dyn_ != nullptr) {
      clk_dyn_->Propagate(clk_, epoch - epoch_);
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

};  // namespace LPT

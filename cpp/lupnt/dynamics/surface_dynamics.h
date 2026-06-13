#pragma once

#include <yaml-cpp/yaml.h>

#include <random>

#include "lupnt/dynamics/dynamics.h"
#include "lupnt/states/state.h"

namespace lupnt {
  // ****************************************************************************
  // Dynamics for Surface Objects (ground stations, rovers, etc.)
  // ****************************************************************************
  class StaticDynamics : public Dynamics {
  public:
    StaticDynamics();
    StaticDynamics(Config& config);
    State Propagate(const State& x0, Real t0, Real tf, const State* u = nullptr) override;
    StateType GetStateType() const override { return State::TYPE; }
  };

  class SurfaceDynamics2D : public Dynamics {
  protected:
    int seed_ = 0;
    bool add_noise_ = false;
    Mat3d Q_ = Mat3d::Zero();  // [km^2/s^2, km^2/s^2, rad^2/s^2]
    Ptr<std::mt19937> rng_ = nullptr;

  public:
    SurfaceDynamics2D();
    SurfaceDynamics2D(Config& config);

    using Dynamics::Propagate;
    State Propagate(const State& x0, Real t0, Real tf, const State* u = nullptr) override;

    void SetSeed(int seed) {
      seed_ = seed;
      rng_ = MakePtr<std::mt19937>(seed_);
      add_noise_ = true;
    }
    int GetSeed() const { return seed_; }

    void SetQ(const Mat3d& Q) {
      Q_ = Q;
      add_noise_ = true;
    }
    Mat3d GetQ() const { return Q_; }

    StateType GetStateType() const override { return SurfaceState2D::TYPE; }
  };

}  // namespace lupnt

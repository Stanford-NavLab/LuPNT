#pragma once

#include "lupnt/core/definitions.h"
#include "lupnt/devices/imu.h"
#include "lupnt/dynamics/dynamics.h"

namespace lupnt {

  enum class ImuModel;

  /// @brief IMU bias random-walk propagation model.
  ///
  /// Propagates an `ImuState` (gyro bias `b_w` and accelerometer bias `b_a`, each
  /// 3-element) as a random walk driven by the bias-instability noise parameters of
  /// the configured `ImuModel` (see `devices/imu.h::Imu::GetParameters`). Used by
  /// `devices/imu.h::Imu` to evolve the "true" IMU bias state that corrupts
  /// simulated IMU measurements consumed by attitude/inertial-navigation filters.
  class ImuDynamics : public Dynamics {
  public:
    /// @brief Construct with an undefined IMU model.
    ImuDynamics();

    /// @brief Construct from a YAML configuration node.
    ///
    /// Reads the optional `model` (`ImuModel`) field.
    ImuDynamics(Config& config);

    /// @brief Set the IMU noise/bias model used to generate process noise.
    void SetModel(ImuModel imu_model);
    /// @brief Return the configured IMU noise/bias model.
    ImuModel GetModel() const;

    /// @brief Set the random-number-generator seed used for bias-noise sampling.
    ///
    /// Reseeds the internal Mersenne Twister RNG (`rng_`).
    void SetSeed(int seed) {
      seed_ = seed;
      rng_ = MakePtr<std::mt19937>(seed_);
    }

    using Dynamics::Propagate;

    /// @brief Propagate the IMU gyro/accelerometer bias state as a random walk over
    /// `[t0, tf]`.
    ///
    /// For each of the 3 gyro-bias and 3 accelerometer-bias components, draws a new
    /// value from a normal distribution centered on the previous bias with standard
    /// deviation `sigma_*_bias / sqrt(1/(tf-t0))` (i.e. scaled by `sqrt(tf-t0)`),
    /// using `imu_model_`'s noise parameters (`Imu::GetParameters`) and the internal
    /// RNG (`rng_`).
    ///
    /// @param x0  Initial IMU bias state (`ImuState`, 6-element: `b_w` then `b_a`)
    ///            at time `t0`.
    /// @param t0  Initial epoch [s].
    /// @param tf  Final epoch [s].
    /// @param u   Unused; must be nullptr (IMU bias dynamics has no control input).
    /// @return    Propagated IMU bias state at time `tf`.
    State Propagate(const State& x0, Real t0, Real tf, const State* u = nullptr) override;

    /// @brief Return `ImuState::TYPE`, the gyro/accelerometer bias state type
    /// propagated by this model.
    StateType GetStateType() const override { return ImuState::TYPE; }

  protected:
    ImuModel imu_model_;
    int seed_ = 0;
    Ptr<std::mt19937> rng_ = nullptr;
  };

}  // namespace lupnt

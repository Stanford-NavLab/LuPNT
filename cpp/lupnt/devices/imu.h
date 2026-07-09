#pragma once

#include "lupnt/devices/device.h"
#include "lupnt/dynamics/imu_dynamics.h"
#include "lupnt/measurements/measurement.h"
#include "lupnt/measurements/measurement_utils.h"
#include "lupnt/states/state.h"

namespace lupnt {

  class ImuDynamics;

  // ImuModel
  enum class ImuModel { LN200S, UNDEFINED };

  // ImuParameters
  struct ImuParameters {
    Real sigma_a;       // [m/s^2 1/sqrt(Hz))] Accelerometer white noise
    Real sigma_a_bias;  // [m/s^2 sqrt(Hz))] Accelerometer bias random walk
    Real sigma_w;       // [rad/s 1/sqrt(Hz))] Gyroscope white noise
    Real sigma_w_bias;  // [rad/s sqrt(Hz))] Gyroscope bias random walk
  };

  /// @brief Onboard inertial measurement unit (IMU) device: simulates gyro
  /// (angular rate) and accelerometer measurements with bias and noise
  /// driven by `ImuDynamics`.
  ///
  /// Attached to a satellite/rover `Agent` (e.g. as the `"imu"` device) to
  /// generate `ImuMeasurement`s that feed `ImuDynamics`-based attitude/state
  /// propagation and navigation filters (EKF/UKF) for inertial state
  /// estimation.
  class Imu : public Device {
  protected:
    ImuModel model_;
    Real time_ = 0.0;  // [s] Time
    ImuState state_;

    std::vector<ImuMeasurement> data_;
    Ptr<ImuDynamics> dynamics_;

    Vec3 position_ = Vec3::Zero();         // [m] Position of the IMU in the body frame
    Mat3 orientation_ = Mat3::Identity();  // Orientation of the IMU in the body frame

  public:
    Imu();
    Imu(Config& config);

    /// @brief Propagate the IMU's bias/state from the last step to time `t`
    /// and generate a new `ImuMeasurement` (gyro + accelerometer reading).
    ///
    /// Called periodically by the `Simulation` event scheduler (per
    /// `Device::Setup`) at the IMU's configured `frequency_`. Propagates
    /// `state_` via `dynamics_->Propagate`, looks up noise parameters for
    /// `model_` via `GetParameters`, samples gyro angular-rate `w` [rad/s]
    /// and accelerometer specific-force `a` [m/s^2] about the current bias
    /// with white-noise scaled by `1/sqrt(t - time_)`, updates `time_`, and
    /// logs the result. Overrides `Device::Step`.
    ///
    /// @param t  Simulation time of this step [s]; must be `>= ` the time of
    ///           the previous `Step`
    void Step(Real t) override;
    /// @brief IMU device setup; currently a no-op.
    void Setup() override;

    /// @brief Get the IMU measurement (gyro angular rate [rad/s] and
    /// accelerometer specific force [m/s^2]) at time `t`.
    ImuMeasurement GetMeasurement(Real t);

    /// @brief Look up the noise/bias-stability parameters (accelerometer and
    /// gyroscope white-noise and bias-random-walk standard deviations) for a
    /// given `ImuModel` (e.g. `LN200S`), used by `Step` to sample
    /// measurement noise.
    ///
    /// @param imu_model  IMU hardware model to look up
    /// @return           Noise parameters: `sigma_a` [m/s^2 /sqrt(Hz)],
    ///                   `sigma_a_bias` [m/s^2 sqrt(Hz)], `sigma_w` [rad/s
    ///                   /sqrt(Hz)], `sigma_w_bias` [rad/s sqrt(Hz)]
    static ImuParameters GetParameters(ImuModel imu_model);

    /// @brief Get all `ImuMeasurement`s accumulated since the last
    /// `EmptyData` call.
    std::vector<ImuMeasurement> GetData() const { return data_; }
    /// @brief Clear the buffered IMU measurement history.
    void EmptyData() { data_.clear(); }

    /// @brief Log the IMU's current state/measurement to the simulation's
    /// `DataLogger` output. Overrides `Device::Log`; currently a no-op.
    void Log(Real time) override;
  };

}  // namespace lupnt

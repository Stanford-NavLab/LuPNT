#include "lupnt/devices/imu.h"

#include <fmt/format.h>

#include "lupnt/core/asset_factory.h"
#include "lupnt/core/logger.h"
#include "lupnt/measurements/measurement.h"
#include "lupnt/numerics/math_utils.h"

namespace lupnt {

  // Imu
  Imu::Imu() : Device() {
    Logger::Debug(fmt::format("Creating {}", name_), "Imu");
    model_ = ImuModel::UNDEFINED;
    dynamics_ = MakePtr<ImuDynamics>();
  }

  Imu::Imu(Config& config) : Device(config) {
    Logger::Debug(fmt::format("Creating {}", name_), "Imu");
    model_ = ImuModel::UNDEFINED;
    dynamics_ = MakePtr<ImuDynamics>();
    if (config["model"]) {
      model_ = enum_cast<ImuModel>(config["model"].as<std::string>()).value();
      dynamics_->SetModel(model_);
    }
    if (config["frequency"]) {
      frequency_ = config["frequency"].as<double>();
      // sim_.Schedule(0.0, [this]() { this->Step(); }, Event::Priority::kSensor, frequency_);
    }
  }

  // Setup
  void Imu::Setup() { Logger::Debug(fmt::format("Setting up {}", name_), "Imu"); }

  // Step
  void Imu::Step(Real time) {
    Logger::Debug("Step", "Imu", time);
    LUPNT_CHECK(time >= time_, "Time must be greater than or equal to the current time", "Imu");
    Real sqrt_f = 1.0 / sqrt(time - time_);
    state_ = dynamics_->Propagate(state_, time_, time);
    ImuParameters params = GetParameters(model_);
    LUPNT_CHECK(false, "Not implemented", "Imu");
    ImuMeasurement measurement;
    for (int i = 0; i < 3; i++) {
      measurement.w(i) = SampleNormal(state_.b_w()(i), params.sigma_w * sqrt_f);
      measurement.a(i) = SampleNormal(state_.b_a()(i), params.sigma_a * sqrt_f);
    }
    time_ = time;
    Log(time);
  }

  // GetImuParameters
  ImuParameters Imu::GetParameters(ImuModel imu_model) {
    // https://stechschulte.net/2023/10/11/imu-specs.html
    constexpr double m_s2_per_ug = 9.80665e-6;  // [m/s^2 / µg]
    switch (imu_model) {
      case ImuModel::LN200S: {
        // https://www.northropgrumman.com/what-we-do/mission-solutions/assured-navigation/ln-200s-inertial-measurement-unit

        // Accelerometer
        double sigma_a = 35.0;                      // [ug / sqrt(Hz)]
        sigma_a = sigma_a * m_s2_per_ug;            // [m/s^2 1/sqrt(Hz)]
        double sigma_a_bias = 300.0;                // [ug sqrt(Hz)]
        sigma_a_bias = sigma_a_bias * m_s2_per_ug;  // [m/s^2 sqrt(Hz)]

        // Gyroscope
        double sigma_w = 0.07;            // [deg/hour 1/sqrt(Hz)]
        sigma_w *= RAD / SECS_HOUR;       // [rad/s 1/sqrt(Hz)]
        double sigma_w_bias = 1.0;        // [deg/hour 1/sqrt(Hz)]
        sigma_w_bias *= RAD / SECS_HOUR;  // [rad/s sqrt(Hz)]

        return ImuParameters{sigma_a, sigma_a_bias, sigma_w, sigma_w_bias};
      }
      default: LUPNT_CHECK(false, "Invalid IMU model", "ImuDynamics");
    }
  }

  // Log
  void Imu::Log(Real time) {}

  // Factory
  REGISTER_FACTORY_CLASS(Device, Imu);
}  // namespace lupnt

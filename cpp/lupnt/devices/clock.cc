#include "lupnt/devices/clock.h"

#include "lupnt/agents/agent.h"
#include "lupnt/conversions/frame_converter.h"
#include "lupnt/core/asset_factory.h"
#include "lupnt/core/data_logger.h"
#include "lupnt/core/definitions.h"
#include "lupnt/core/error.h"
#include "lupnt/core/logger.h"
#include "lupnt/environment/body.h"

namespace lupnt {

  // Clock
  Clock::Clock() : Device() {
    Logger::Debug(fmt::format("Creating {}", name_), "Clock");
    dynamics_ = MakePtr<ClockDynamics>();
    ApplyClockStateUnits();
  }

  // Constructor
  Clock::Clock(Config& config) : Device(config) {
    Logger::Debug(fmt::format("Creating {}", name_), "Clock");
    dynamics_ = MakePtr<ClockDynamics>();
    if (config["model"]) {
      model_ = enum_cast<ClockModel>(config["model"].as<std::string>()).value();
      dynamics_->SetModel(model_);
    }
    if (config["clock_bias_unit"]) {
      SetClockBiasUnit(
          enum_cast<ClockBiasUnit>(config["clock_bias_unit"].as<std::string>()).value());
    } else if (config["bias_unit"]) {
      SetClockBiasUnit(enum_cast<ClockBiasUnit>(config["bias_unit"].as<std::string>()).value());
    } else {
      ApplyClockStateUnits();
    }
    if (config["use_relativity"]) use_relativity_ = config["use_relativity"].as<bool>();
    if (config["relativity_center_body"]) {
      SetRelativityCenterBody(
          enum_cast<BodyId>(config["relativity_center_body"].as<std::string>()).value());
    }
  }

  void Clock::ApplyClockStateUnits() {
    state_.SetUnits(ClockDynamics::GetStateUnits(state_.size(), bias_unit_));
  }

  void Clock::SetModel(ClockModel model) {
    model_ = model;
    if (dynamics_) dynamics_->SetModel(model_);
  }

  void Clock::SetClockBiasUnit(ClockBiasUnit unit) {
    bias_unit_ = unit;
    if (dynamics_) dynamics_->SetClockBiasUnit(unit);
    ApplyClockStateUnits();
  }

  void Clock::SetDynamics(Ptr<ClockDynamics> dynamics) {
    dynamics_ = std::move(dynamics);
    if (dynamics_) {
      dynamics_->SetModel(model_);
      dynamics_->SetClockBiasUnit(bias_unit_);
    }
  }

  ClockRelativityContext Clock::BuildRelativityContext(Real t0, Real tf) const {
    LUPNT_CHECK(agent_ != nullptr, "Relativistic clock propagation requires an owning agent",
                "Clock");

    Cart6 rv0 = agent_->GetStateAt(t0);
    Cart6 rvf = agent_->GetStateAt(tf);
    BodyId center_body
        = relativity_center_body_set_ ? relativity_center_body_ : GetFrameCenter(rvf.GetFrame());
    BodyData center_data = GetBodyData(center_body, relativity_units_);
    Frame center_inertial_frame = center_data.inertial_frame;

    Cart6 rv0_centered = rv0;
    if (rv0.GetFrame() != center_inertial_frame) {
      rv0_centered = Cart6(ConvertFrame(t0, rv0, center_inertial_frame));
    }
    Cart6 rvf_centered = rvf;
    if (rvf.GetFrame() != center_inertial_frame) {
      rvf_centered = Cart6(ConvertFrame(tf, rvf, center_inertial_frame));
    }

    PhysicalConstants constants = GetPhysicalConstants(relativity_units_);
    ClockRelativityContext context;
    context.rv0_centered = rv0_centered;
    context.rvf_centered = rvf_centered;
    context.GM = center_data.GM;
    context.c = constants.C;
    return context;
  }

  // ReadTime
  Real Clock::Read(Real t) {
    Logger::Debug(fmt::format("Reading time {}", t), "Clock");
    LUPNT_CHECK(t >= time_, "Time must be greater than or equal to the current time", "Clock");
    if (use_relativity_) {
      ClockRelativityContext relativity = BuildRelativityContext(time_, t);
      state_ = dynamics_->PropagateWithRelativity(state_, time_, t, relativity);
    } else {
      state_ = dynamics_->Propagate(state_, time_, t);
    }
    time_ = t;
    return time_ + ClockDynamics::BiasUnitsToSeconds(state_.b(), bias_unit_);
  }

  // Step
  void Clock::Step(Real t) { Logger::Debug("Step", "Clock", t); }

  // Setup
  void Clock::Setup() { Logger::Debug(fmt::format("Setting up {}", name_), "Clock"); }

  // Log
  void Clock::Log(Real time) {
    DataLogger::Log(fmt::format("{}/time", name_), time);
    for (int i = 0; i < state_.size(); i++) {
      auto name = state_.GetNames()[i];
      auto unit = state_.GetUnits()[i];
      std::replace(unit.begin(), unit.end(), '/', '_');
      DataLogger::Log(fmt::format("{}/state/{}_{}", name_, name, unit), state_(i));
    }
  }

  // Factory
  REGISTER_FACTORY_CLASS(Device, Clock);
}  // namespace lupnt

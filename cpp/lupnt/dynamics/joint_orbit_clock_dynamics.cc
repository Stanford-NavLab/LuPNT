#include "lupnt/dynamics/joint_orbit_clock_dynamics.h"

#include <algorithm>
#include <limits>

#include "lupnt/core/asset_factory.h"
#include "lupnt/core/constants.h"
#include "lupnt/core/error.h"
#include "lupnt/data/kernels.h"
#include "lupnt/devices/clock.h"
#include "lupnt/environment/solar_system.h"

namespace lupnt {

  namespace {
    Ptr<NumericalDynamics> CreateNumericalDynamics(Config& config) {
      std::string dyn_class = config["class"].as<std::string>();
      Ptr<Dynamics> dynamics = AssetFactory<Dynamics, Config&>::Create(dyn_class, config);
      Ptr<NumericalDynamics> numerical = std::dynamic_pointer_cast<NumericalDynamics>(dynamics);
      LUPNT_CHECK(numerical != nullptr, "Orbit dynamics must derive from NumericalDynamics",
                  "JointOrbitClockDynamics");
      return numerical;
    }
  }  // namespace

  JointOrbitClockDynamics::JointOrbitClockDynamics() : NumericalDynamics() {
    SetODE([this](Real t, const VecX& x) { return ComputeRates(t, State(x)); });
    clock_dynamics_ = MakePtr<ClockDynamics>();
  }

  JointOrbitClockDynamics::JointOrbitClockDynamics(Config& config) : NumericalDynamics(config) {
    SetODE([this](Real t, const VecX& x) { return ComputeRates(t, State(x)); });
    clock_dynamics_ = MakePtr<ClockDynamics>();

    if (config["orbit_dynamics"]) {
      Config orbit_config(config["orbit_dynamics"]);
      SetOrbitDynamics(CreateNumericalDynamics(orbit_config));
    }
    if (config["clock_dynamics"]) {
      Config clock_config(config["clock_dynamics"]);
      SetClockDynamics(MakePtr<ClockDynamics>(clock_config));
    }
    if (config["clock_bias_unit"]) {
      clock_dynamics_->SetClockBiasUnit(
          enum_cast<ClockBiasUnit>(config["clock_bias_unit"].as<std::string>()).value());
    } else if (config["bias_unit"]) {
      clock_dynamics_->SetClockBiasUnit(
          enum_cast<ClockBiasUnit>(config["bias_unit"].as<std::string>()).value());
    }
    if (config["use_clock_relativity"]) {
      use_clock_relativity_ = config["use_clock_relativity"].as<bool>();
    } else if (config["use_relativity"]) {
      use_clock_relativity_ = config["use_relativity"].as<bool>();
    }
    if (config["relativity_center_body"]) {
      SetRelativityCenterBody(
          enum_cast<BodyId>(config["relativity_center_body"].as<std::string>()).value());
    }
    if (config["reference_rate_offset"]) {
      reference_rate_offset_ = config["reference_rate_offset"].as<double>();
    }
    if (config["add_clock_noise"]) add_clock_noise_ = config["add_clock_noise"].as<bool>();
    if (config["frame"]) frame_ = enum_cast<Frame>(config["frame"].as<std::string>()).value();
  }

  void JointOrbitClockDynamics::SetOrbitDynamics(Ptr<NumericalDynamics> orbit_dynamics) {
    LUPNT_CHECK(orbit_dynamics != nullptr, "Orbit dynamics cannot be null",
                "JointOrbitClockDynamics");
    orbit_dynamics_ = std::move(orbit_dynamics);
    SetParams(orbit_dynamics_->GetParams());

    if (auto nbody = std::dynamic_pointer_cast<NBodyDynamics>(orbit_dynamics_)) {
      frame_ = nbody->GetFrame();
      units_ = nbody->GetUnits();
    }
  }

  void JointOrbitClockDynamics::SetClockDynamics(Ptr<ClockDynamics> clock_dynamics) {
    LUPNT_CHECK(clock_dynamics != nullptr, "Clock dynamics cannot be null",
                "JointOrbitClockDynamics");
    clock_dynamics_ = std::move(clock_dynamics);
  }

  BodyId JointOrbitClockDynamics::SelectRelativityCenter(Real t, const State& orbit_state) const {
    if (relativity_center_body_set_) return relativity_center_body_;
    if (frame_ != Frame::UNDEFINED) {
      BodyId frame_center = GetFrameCenter(frame_);
      if (frame_center != BodyId::SSB && frame_center != BodyId::SUN) return frame_center;
    }

    auto nbody = std::dynamic_pointer_cast<NBodyDynamics>(orbit_dynamics_);
    if (!nbody) return BodyId::EARTH;

    BodyId closest = BodyId::EARTH;
    Real closest_distance = std::numeric_limits<double>::infinity();
    for (const Body& body : nbody->GetBodies()) {
      if (body.id == BodyId::SSB || body.id == BodyId::SUN) continue;
      Vec6 rv_center = GetCenterState(t, body.id);
      Real distance = (orbit_state.head(3) - rv_center.head(3)).norm();
      if (distance < closest_distance) {
        closest_distance = distance;
        closest = body.id;
      }
    }
    return closest;
  }

  Vec6 JointOrbitClockDynamics::GetCenterState(Real t, BodyId center_body) const {
    LUPNT_CHECK(frame_ != Frame::UNDEFINED, "Frame must be set for relativistic clock coupling",
                "JointOrbitClockDynamics");
    BodyId frame_center = GetFrameCenter(frame_);
    if (frame_center == center_body) return Vec6::Zero();
    Real t_tdb = t + GetLupntEpoch();
    return GetBodyPosVel(t_tdb, center_body, frame_, units_);
  }

  Real JointOrbitClockDynamics::RelativisticRate(Real t, const State& orbit_state) const {
    BodyId center_body = SelectRelativityCenter(t, orbit_state);
    BodyData center_data = GetBodyData(center_body, units_);
    PhysicalConstants constants = GetPhysicalConstants(units_);
    Vec6 rv_center = GetCenterState(t, center_body);
    return ClockDynamics::RelativisticRateCorrection(
        orbit_state.head(3) - rv_center.head(3), orbit_state.tail(3) - rv_center.tail(3),
        center_data.GM, constants.C, reference_rate_offset_);
  }

  State JointOrbitClockDynamics::BuildState(const State& x0, const VecX& values) const {
    JointOrbitClockState xf(x0);
    xf.head(values.size()) = values;
    std::vector<std::string> units = xf.GetUnits();
    std::vector<std::string> clock_units = ClockDynamics::GetStateUnits(
        values.size() - ORBIT_STATE_SIZE, clock_dynamics_->GetClockBiasUnit());
    std::copy(clock_units.begin(), clock_units.end(), units.begin() + ORBIT_STATE_SIZE);
    xf.SetUnits(units);
    if (frame_ != Frame::UNDEFINED) xf.SetFrame(frame_);
    return xf;
  }

  State JointOrbitClockDynamics::Propagate(const State& x0, Real t0, Real tf, const State* u) {
    LUPNT_CHECK(orbit_dynamics_ != nullptr, "Orbit dynamics must be set",
                "JointOrbitClockDynamics");
    LUPNT_CHECK(clock_dynamics_ != nullptr, "Clock dynamics must be set",
                "JointOrbitClockDynamics");
    LUPNT_CHECK(u == nullptr, "JointOrbitClockDynamics does not support control inputs",
                "JointOrbitClockDynamics");
    if (std::abs((tf - t0).val()) < EPS) return BuildState(x0, x0);

    State x_prop = NumericalDynamics::Propagate(x0, t0, tf, nullptr);
    VecX values = x_prop;
    int clock_state_size = x0.size() - ORBIT_STATE_SIZE;
    if (add_clock_noise_ && clock_dynamics_->GetModel() != ClockModel::UNDEFINED) {
      values.tail(clock_state_size) += clock_dynamics_->GetProcessNoise(tf - t0, clock_state_size);
    }
    return BuildState(x0, values);
  }

  State JointOrbitClockDynamics::PropagateWithParams(const State& x0, Real t0, Real tf,
                                                     const ParamState& params, const State* u) {
    LUPNT_CHECK(orbit_dynamics_ != nullptr, "Orbit dynamics must be set",
                "JointOrbitClockDynamics");
    orbit_dynamics_->SetParams(params);
    SetParams(params);
    return Propagate(x0, t0, tf, u);
  }

  VecX JointOrbitClockDynamics::ComputeRates(Real t, const State& x) const {
    LUPNT_CHECK(orbit_dynamics_ != nullptr, "Orbit dynamics must be set",
                "JointOrbitClockDynamics");
    LUPNT_CHECK(clock_dynamics_ != nullptr, "Clock dynamics must be set",
                "JointOrbitClockDynamics");
    LUPNT_CHECK(x.size() == 8 || x.size() == 9, "State size must be 8 or 9",
                "JointOrbitClockDynamics");

    int clock_state_size = x.size() - ORBIT_STATE_SIZE;
    State orbit_state = x.head(ORBIT_STATE_SIZE);
    orbit_state.SetFrame(frame_ == Frame::UNDEFINED ? x.GetFrame() : frame_);

    VecX orbit_rates = orbit_dynamics_->ComputeRates(t, orbit_state);
    Real rel_rate = 0.0;
    if (use_clock_relativity_) {
      rel_rate = RelativisticRate(t, orbit_state);
    }
    Real rel_bias_rate
        = ClockDynamics::SecondsToBiasUnits(rel_rate, clock_dynamics_->GetClockBiasUnit());

    VecX clock_rates(clock_state_size);
    clock_rates(0) = x(ORBIT_STATE_SIZE + 1) + rel_bias_rate;
    if (clock_state_size == 2) {
      clock_rates(1) = 0.0;
    } else {
      clock_rates(1) = x(ORBIT_STATE_SIZE + 2);
      clock_rates(2) = 0.0;
    }

    VecX rates(x.size());
    rates.head(ORBIT_STATE_SIZE) = orbit_rates;
    rates.tail(clock_state_size) = clock_rates;
    return rates;
  }

  REGISTER_FACTORY_CLASS(Dynamics, JointOrbitClockDynamics)

}  // namespace lupnt

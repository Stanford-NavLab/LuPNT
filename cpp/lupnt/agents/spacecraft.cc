#include "lupnt/agents/spacecraft.h"

#include "lupnt/applications/application.h"
#include "lupnt/conversions/state_converter.h"
#include "lupnt/core/asset_factory.h"
#include "lupnt/core/data_logger.h"
#include "lupnt/devices/clock.h"  // ClockModel, ClockBiasUnit
#include "lupnt/dynamics/clock_dynamics.h"
#include "lupnt/dynamics/numerical_orbit_dynamics.h"
#include "lupnt/environment/body.h"

namespace lupnt {

  Spacecraft::Spacecraft(Config& config) : Agent(config) {
    Logger::Debug(fmt::format("Creating Spacecraft {}", name_), "Spacecraft");
    ConfigureTruth(config);
  }

  void Spacecraft::ConfigureTruth(Config& config) {
    // --- Orbit truth: from the agent's `dynamics:` block (the Simulation injects
    // the shared world force model here when the agent declares none). ---
    LUPNT_CHECK(config["dynamics"],
                "Spacecraft needs a `dynamics:` block (or a `world.force_model` to inherit)",
                "Spacecraft");
    Config dyn(config["dynamics"]);
    std::string dyn_class = dyn["class"].as<std::string>("NBodyDynamics");
    if (!dyn["name"]) dyn["name"] = name_ + "/orbit_truth";
    auto base = AssetFactory<Dynamics, Config&>::Create(dyn_class, dyn);
    auto orbit = std::dynamic_pointer_cast<NBodyDynamics>(base);
    LUPNT_CHECK(orbit, "Spacecraft `dynamics:` must be an NBodyDynamics-style force model",
                "Spacecraft");
    orbit->SetAutodiff(false);  // truth needs no state-transition matrix
    Frame frame = orbit->GetFrame();
    double step = dyn["dt"].as<double>(60.0);

    // --- Clock truth: from the optional `clock:` block. ---
    auto clock = MakePtr<ClockDynamics>();
    Config clk = config["clock"] ? Config(config["clock"]) : Config(YAML::Node());
    ClockModel model = ClockModel::OCXO;
    if (clk["model"]) model = enum_cast<ClockModel>(clk["model"].as<std::string>()).value();
    int clock_seed = clk["seed"].as<int>(config["seed"].as<int>(0));
    bool add_noise = clk["add_noise"].as<bool>(true);
    clock->SetModel(model);
    clock->SetClockBiasUnit(ClockBiasUnit::SECONDS);
    clock->SetAddNoise(add_noise);
    clock->SetSeed(clock_seed);

    truth_dyn_ = MakePtr<JointOrbitClockDynamics>();
    truth_dyn_->SetOrbitDynamics(orbit);
    truth_dyn_->SetClockDynamics(clock);
    truth_dyn_->SetUseClockRelativity(true);
    truth_dyn_->SetRelativityCenterBody(BodyId::MOON);
    truth_dyn_->SetAddClockNoise(add_noise);
    truth_dyn_->SetFrame(frame);
    truth_dyn_->SetTimeStep(step);
    truth_dyn_->SetIntegrator(IntegratorType::RKF45);
    truth_dyn_->SetIntegratorParams(IntegratorParams(20, 1.0e-12, 1.0e-12));

    // --- Initial truth 8-state [r, v, clock_bias, clock_drift] (Frame::MOON_CI). ---
    Config init(config["initial_state"]);
    Cart6 orbit0;
    if (init["r0_m"]) {  // explicit Cartesian in MOON_CI
      Vec3 r0(init["r0_m"][0].as<Real>(), init["r0_m"][1].as<Real>(), init["r0_m"][2].as<Real>());
      Vec3 v0(init["v0_mps"][0].as<Real>(), init["v0_mps"][1].as<Real>(),
              init["v0_mps"][2].as<Real>());
      orbit0 = Cart6(r0, v0, Frame::MOON_CI);
    } else if (init["class"] && init["class"].as<std::string>() == "ClassicalOE") {
      Real a = init["a"].as<Real>(), e = init["e"].as<Real>(), i = RAD * init["i"].as<Real>();
      Real Om = RAD * init["Omega"].as<Real>(), om = RAD * init["omega"].as<Real>();
      Real M = RAD * init["M"].as<Real>();
      Frame f = enum_cast<Frame>(init["frame"].as<std::string>()).value();
      State coe = ClassicalOE(Vec6(a, e, i, Om, om, M), f);
      State rv
          = ConvertFrame(GetLupntEpoch(), ConvertState(coe, orbit->GetStateType()), Frame::MOON_CI);
      orbit0 = Cart6(Vec6(rv.head(6)), Frame::MOON_CI);
    } else {
      LUPNT_CHECK(false, "Spacecraft `initial_state:` needs `r0_m`/`v0_mps` or a ClassicalOE block",
                  "Spacecraft");
    }
    ClockState2 clock0;
    clock0.b() = init["clock_bias_s"].as<Real>(0.0);
    clock0.d() = init["clock_drift_sps"].as<Real>(0.0);
    truth_ = JointOrbitClockState(orbit0, clock0);
    truth_time_ = 0.0;
  }

  void Spacecraft::Step(Real t) {
    if (truth_dyn_ && abs(t - truth_time_) > EPS) {
      truth_ = truth_dyn_->Propagate(truth_, truth_time_, t, nullptr);
      truth_time_ = t;
    }
    DataLogger::Log(fmt::format("{}/truth", name_), truth_);
    if (application_) application_->Log(t);
  }

  State Spacecraft::GetTruthStateAt(Real t) const {
    if (abs(t - truth_time_) < EPS) return truth_;
    LUPNT_CHECK(truth_dyn_, "Spacecraft truth dynamics not set", "Spacecraft");
    return truth_dyn_->Propagate(truth_, truth_time_, t, nullptr);
  }

  Cart6 Spacecraft::GetStateAt(Real t) const {
    State x = GetTruthStateAt(t);
    return Cart6(Vec6(x.head(6)), Frame::MOON_CI);
  }

  REGISTER_FACTORY_CLASS(Agent, Spacecraft)

}  // namespace lupnt

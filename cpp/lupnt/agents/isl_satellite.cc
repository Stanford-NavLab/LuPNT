#include "lupnt/agents/isl_satellite.h"

#include "lupnt/applications/application.h"
#include "lupnt/core/asset_factory.h"
#include "lupnt/core/data_logger.h"
#include "lupnt/devices/clock.h"  // ClockModel, ClockBiasUnit (clock_dynamics.h only forward-declares)
#include "lupnt/dynamics/clock_dynamics.h"
#include "lupnt/dynamics/numerical_orbit_dynamics.h"
#include "lupnt/environment/body.h"

namespace lupnt {

  IslSatellite::IslSatellite(Config& config) : Agent(config) {
    // --- Truth force model (orbit) + clock, forming the joint truth dynamics. ---
    Config truth = config["truth"];
    int deg = truth["moon_gravity_degree"].as<int>(16);
    int ord = truth["moon_gravity_order"].as<int>(16);
    bool earth = truth["include_earth"].as<bool>(true);
    bool sun = truth["include_sun"].as<bool>(true);
    bool use_rel = truth["use_relativity"].as<bool>(true);
    double step = truth["integration_step_s"].as<double>(60.0);
    int clock_seed = truth["clock_seed"].as<int>(0);

    auto orbit = MakePtr<NBodyDynamics>();
    orbit->SetIntegrator(IntegratorType::RKF45);
    orbit->SetIntegratorParams(IntegratorParams(20, 1.0e-12, 1.0e-12));
    orbit->AddBody(Body::Moon(deg, ord));
    if (earth) orbit->AddBody(Body::Earth());
    if (sun) orbit->AddBody(Body::Sun());
    orbit->SetFrame(Frame::MOON_CI);
    orbit->SetTimeStep(step);
    orbit->SetAutodiff(false);  // truth needs no STM
    orbit->SetUseRelativity(use_rel);

    auto clock = MakePtr<ClockDynamics>();
    clock->SetModel(ClockModel::OCXO);
    clock->SetClockBiasUnit(ClockBiasUnit::SECONDS);
    clock->SetAddNoise(true);
    clock->SetSeed(clock_seed);

    truth_dyn_ = MakePtr<JointOrbitClockDynamics>();
    truth_dyn_->SetOrbitDynamics(orbit);
    truth_dyn_->SetClockDynamics(clock);
    truth_dyn_->SetUseClockRelativity(true);
    truth_dyn_->SetRelativityCenterBody(BodyId::MOON);
    truth_dyn_->SetAddClockNoise(true);
    truth_dyn_->SetFrame(Frame::MOON_CI);
    truth_dyn_->SetTimeStep(step);
    truth_dyn_->SetIntegrator(IntegratorType::RKF45);
    truth_dyn_->SetIntegratorParams(IntegratorParams(20, 1.0e-12, 1.0e-12));

    // --- Initial truth 8-state [r, v, clock_bias, clock_drift] (Frame::MOON_CI). ---
    Config init = config["initial_state"];
    Vec3 r0(init["r0_m"][0].as<Real>(), init["r0_m"][1].as<Real>(), init["r0_m"][2].as<Real>());
    Vec3 v0(init["v0_mps"][0].as<Real>(), init["v0_mps"][1].as<Real>(),
            init["v0_mps"][2].as<Real>());
    Cart6 orbit0(r0, v0, Frame::MOON_CI);
    ClockState2 clock0;
    clock0.b() = init["clock_bias_s"].as<Real>(0.0);
    clock0.d() = init["clock_drift_sps"].as<Real>(0.0);
    truth_ = JointOrbitClockState(orbit0, clock0);
    truth_time_ = 0.0;
  }

  void IslSatellite::Step(Real t) {
    if (truth_dyn_ && abs(t - truth_time_) > EPS) {
      truth_ = truth_dyn_->Propagate(truth_, truth_time_, t, nullptr);
      truth_time_ = t;
    }
    DataLogger::Log(fmt::format("{}/truth", name_), truth_);
    if (application_) application_->Log(t);
  }

  State IslSatellite::GetTruthStateAt(Real t) const {
    if (abs(t - truth_time_) < EPS) return truth_;
    LUPNT_CHECK(truth_dyn_, "IslSatellite truth dynamics not set", "IslSatellite");
    return truth_dyn_->Propagate(truth_, truth_time_, t, nullptr);
  }

  Cart6 IslSatellite::GetStateAt(Real t) const {
    State x = GetTruthStateAt(t);
    return Cart6(Vec6(x.head(6)), Frame::MOON_CI);
  }

  REGISTER_FACTORY_CLASS(Agent, IslSatellite)

}  // namespace lupnt

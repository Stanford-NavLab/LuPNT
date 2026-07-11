#include "lupnt/applications/lunar_gnss_odts/lunar_gnss_odts_app.h"

#include <filesystem>

#include "lupnt/agents/agent.h"
#include "lupnt/agents/spacecraft.h"
#include "lupnt/core/asset_factory.h"
#include "lupnt/lupnt.h"
#include "lupnt/simulations/simulation.h"
#include "lupnt/simulations/world.h"

namespace lupnt {

  LunarGnssODTSConfig ConfigToLunarGnssODTSConfig(Config& config) {
    // The `application:` block carries the same nested sections as the scenario file; relative
    // data/output paths are resolved against the current working directory (run scripts inject
    // absolute paths).
    return ParseLunarGnssODTSConfig(config, std::filesystem::current_path());
  }

  LunarGnssOdtsApp::LunarGnssOdtsApp(Config& config) : Application(config) {
    cfg_ = ConfigToLunarGnssODTSConfig(config);
  }

  LunarGnssOdtsApp::LunarGnssOdtsApp(LunarGnssODTSConfig config) : cfg_(std::move(config)) {
    ResolveLunarGnssODTSConfigForRun(cfg_);
  }

  void LunarGnssOdtsApp::Setup() {
    LUPNT_CHECK(agent_, "Agent not set", "LunarGnssOdtsApp");
    // Plasma is a shared truth-environment property: prefer the `world: plasma:` block over any
    // legacy `application: plasma:` block (parsed into cfg_ at construction). The World holds the
    // raw node; parse it into the config so the truth-generation + filter paths see it.
    if (World* world = agent_->GetWorld(); world && world->HasPlasma()) {
      ParsePlasmaDelayConfig(world->GetPlasma(), cfg_.plasma);
    }
    // Defer the heavy engine to the first Step (a single event); the base Run() drains it.
    Simulation* sim = agent_->GetSimulation();
    sim->Schedule(
        0.0, [this](Real t) { Step(t); }, Event::SINGLE_EVENT, Event::Priority::APPLICATION);
  }

  void LunarGnssOdtsApp::Precompute() {
    // The link precompute time-keeps in ABSOLUTE TDB (like RunAll), and the orbit propagation
    // adds `GetLupntEpoch()` to its time argument. When this app is hosted in a `Simulation`
    // (which sets the global epoch from the scenario `epoch:`), a cold recompute would then
    // double-count the epoch and request ephemerides far outside the loaded kernels. Reset the
    // epoch to 0 for the precompute (matching RunAll and the standalone struct path), then
    // restore it so the rest of the run is unaffected.
    Real saved_epoch = GetLupntEpoch();
    SetLupntEpoch(0.0);
    PrecomputeLunarGnssODTSLinks(cfg_);
    SetLupntEpoch(saved_epoch);
  }

  void LunarGnssOdtsApp::RunAll() {
    // Agent-driven truth: when hosted on a physical `Spacecraft`, sample its self-propagated
    // orbit+clock truth on the engine's receiver epoch grid (at the real epoch, *before* the
    // reset below) and hand it to the engine, which then reads the agent's truth instead of
    // building its own. Falls back to the internal build on a non-Spacecraft host (legacy).
    const std::vector<State>* truth_ptr = nullptr;
    std::vector<State> truth;
    if (auto* sc = dynamic_cast<Spacecraft*>(agent_)) {
      VecXd elapsed = LunarGnssODTSReceiverElapsedTimes(cfg_);
      auto* dyn = sc->GetTruthDynamics();
      LUPNT_CHECK(dyn, "Spacecraft has no truth dynamics", "LunarGnssOdtsApp");
      State x = sc->GetTruthState();
      truth.reserve(elapsed.size());
      for (int i = 0; i < elapsed.size(); ++i) {
        if (i > 0) x = dyn->Propagate(JointOrbitClockState(x), elapsed(i - 1), elapsed(i), nullptr);
        truth.push_back(x);
      }
      truth_ptr = &truth;
    }

    // The engine time-keeps in ABSOLUTE TDB seconds past J2000, passing those directly to the
    // orbit dynamics; orbit propagation adds `GetLupntEpoch()` to its time argument, so we
    // reset the global epoch to 0 here (the base `Simulation` sets it from the scenario
    // `epoch:`) to avoid double-counting. The sampled truth is already at the correct absolute
    // epochs, so this reset only affects the filter's estimate propagation.
    SetLupntEpoch(0.0);
    summaries_ = RunLunarGnssODTSMonteCarlo(cfg_, truth_ptr);
  }

  void LunarGnssOdtsApp::Step(Real /*t*/) {
    if (ran_) return;
    ran_ = true;
    RunAll();
  }

  REGISTER_FACTORY_CLASS(Application, LunarGnssOdtsApp)

}  // namespace lupnt

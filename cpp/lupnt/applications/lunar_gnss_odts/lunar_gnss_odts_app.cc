#include "lupnt/applications/lunar_gnss_odts/lunar_gnss_odts_app.h"

#include <filesystem>

#include "lupnt/agents/agent.h"
#include "lupnt/core/asset_factory.h"
#include "lupnt/lupnt.h"
#include "lupnt/simulations/simulation.h"

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
    // Defer the heavy engine to the first Step (a single event); the base Run() drains it.
    Simulation* sim = agent_->GetSimulation();
    sim->Schedule(
        0.0, [this](Real t) { Step(t); }, Event::SINGLE_EVENT, Event::Priority::APPLICATION);
  }

  void LunarGnssOdtsApp::Precompute() { PrecomputeLunarGnssODTSLinks(cfg_); }

  void LunarGnssOdtsApp::RunAll() {
    // The engine time-keeps in ABSOLUTE TDB seconds past J2000, passing those directly to the
    // orbit dynamics; orbit propagation adds `GetLupntEpoch()` to its time argument, so we
    // reset the global epoch to 0 here (the base `Simulation` sets it from the scenario
    // `epoch:`) to avoid double-counting and reproduce the former monolith bit-for-bit.
    SetLupntEpoch(0.0);
    summaries_ = RunLunarGnssODTSMonteCarlo(cfg_);
  }

  void LunarGnssOdtsApp::Step(Real /*t*/) {
    if (ran_) return;
    ran_ = true;
    RunAll();
  }

  REGISTER_FACTORY_CLASS(Application, LunarGnssOdtsApp)

}  // namespace lupnt

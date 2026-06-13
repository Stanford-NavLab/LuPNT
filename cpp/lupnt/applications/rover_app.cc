#include "lupnt/applications/rover_app.h"

#include "lupnt/agents/agent.h"
#include "lupnt/core/asset_factory.h"
#include "lupnt/core/logger.h"
#include "lupnt/devices/clock.h"
#include "lupnt/filters/filter.h"
#include "lupnt/states/state.h"

namespace lupnt {

  // Constructor
  RoverApp::RoverApp(Config& config) : Application(config) {
    Logger::Debug(fmt::format("Creating {}", name_), "RoverApp");

    // Filter dynamics
    LUPNT_CHECK(config["dynamics"], "No dynamics found", "RoverApp");
    Config dyn_config(config["dynamics"]);
    if (!dyn_config["name"]) dyn_config["name"] = "dynamics";
    dyn_config["name"] = name_ + "/" + dyn_config["name"].as<std::string>();
    std::string dynamics_class = dyn_config["class"].as<std::string>();
    dynamics_ = AssetFactory<Dynamics, Config&>::Create(dynamics_class, dyn_config);

    // Filter
    LUPNT_CHECK(config["filter"], "No filter found", "RoverApp");
    Config filter_config(config["filter"]);
    if (!filter_config["name"]) filter_config["name"] = "filter";
    filter_config["name"] = name_ + "/" + filter_config["name"].as<std::string>();
    std::string filter_class = filter_config["class"].as<std::string>();
    filter_ = AssetFactory<Filter, Config&>::Create(filter_class, filter_config);
  }

  // Setup
  void RoverApp::Setup() {
    Logger::Debug(fmt::format("Setting up {}", name_), "RoverApp");

    // Checks
    LUPNT_CHECK(filter_, "Filter not set", "RoverApp");
    LUPNT_CHECK(config_["initial_state"], "No initial state found", "RoverApp");
    LUPNT_CHECK(agent_, "Agent not set", "RoverApp");

    // Initial state
    Real sigma_r = config_["initial_state"]["sigma_r"].as<Real>();
    Real sigma_theta = RAD * config_["initial_state"]["sigma_theta"].as<Real>();
    MatXd P0 = Mat3d::Zero();
    P0(0, 0) = sigma_r * sigma_r;
    P0(1, 1) = sigma_r * sigma_r;
    P0(2, 2) = sigma_theta * sigma_theta;

    // Filter
    Real t = agent_->GetTime();
    Cart6 rv = agent_->GetState();

    // Dynamics function
    FilterDynamicsFunction f_dyn
        = [this](const State& x, Real t0, Real tf, const State* u, MatXd* F) {
            return dynamics_->Propagate(x, t0, tf, u, F);
          };

    // Process noise function
    MatXd Q = P0;
    ProcessNoiseFunction f_proc = [Q](const State& x, Real t0, Real tf) {
      (void)x;
      (void)t0;
      (void)tf;
      return Q;
    };

    filter_->SetTime(t);
    filter_->SetState(rv);
    filter_->SetCovariance(P0);
    filter_->SetDynamicsFunction(f_dyn);
    filter_->SetProcessNoiseFunction(f_proc);

    // Log after setup
    Log(0.0);
  }

  // Step
  void RoverApp::Step(Real t) {
    Logger::Debug("Step", "RoverApp", t);
    LUPNT_CHECK(agent_, "Agent not set", "RoverApp");
    LUPNT_CHECK(filter_, "Filter not set", "RoverApp");

    // auto clk = dynamic_cast<Clock*>(agent_->GetDevice("clock"))->Read(t);

    // Control
    // auto u = agent_->GetControl();
    State u = Vec2::Zero();
    if (u.size() == 0) {
      u = State::Zero(2);
      u.SetNames({"v", "w"});
      u.SetUnits({"m/s", "rad/s"});
    }
    // filter_->Predict(t, &u);

    // Control
    Real R = config_["trajectory"]["R"].as<Real>();
    Real v = config_["trajectory"]["v"].as<Real>();
    Real w = v / R;
    State u_new(2);
    u_new << v, w;
    u_new.SetNames({"v", "w"});
    u_new.SetUnits({"m/s", "rad/s"});
    // agent_->SetControl(u_new);
    Logger::Debug(fmt::format("Control: {}", u_new), "RoverApp");

    Log(t);
  }

  // Log
  void RoverApp::Log(Real t) {
    Logger::Debug(fmt::format("Logging {}", name_), "RoverApp", t);
    filter_->Log(t);

    // Errors
  }

  REGISTER_FACTORY_CLASS(Application, RoverApp)

}  // namespace lupnt

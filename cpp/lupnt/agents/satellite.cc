#include "lupnt/agents/satellite.h"

#include "lupnt/conversions/state_converter.h"
#include "lupnt/core/constants.h"
#include "lupnt/core/definitions.h"
#include "lupnt/core/logger.h"
#include "lupnt/core/simulation.h"
#include "lupnt/states/state.h"

namespace lupnt {

  // Satellite
  Satellite::Satellite(Config& agent_config) : AgentWithDynamics(agent_config) {
    Logger::Debug(fmt::format("Creating Satellite {}", name_), "Satellite");

    // Initial state
    if (agent_config["initial_state"]) {
      auto cfg = agent_config["initial_state"];
      if (cfg["class"].as<std::string>() == "ClassicalOE") {
        Real a = cfg["a"].as<Real>();
        Real e = cfg["e"].as<Real>();
        Real i = RAD * cfg["i"].as<Real>();
        Real Omega = RAD * cfg["Omega"].as<Real>();
        Real omega = RAD * cfg["omega"].as<Real>();
        Real M = RAD * cfg["M"].as<Real>();

        Frame frame = enum_cast<Frame>(cfg["frame"].as<std::string>()).value();
        State state = ClassicalOE(Vec6(a, e, i, Omega, omega, M), frame);
        Logger::Debug(fmt::format("Initial state:\n{}", state), "Satellite");

        state_ = ConvertFrame(GetLupntEpoch(), ConvertState(state, dynamics_->GetStateType()),
                              Frame::MOON_CI);
        Logger::Debug(fmt::format("Setting state:\n{}", state_), "Satellite");
      } else {
        LUPNT_CHECK(false, "Unsupported initial state", "Satellite");
      }
    }
  }

  void Satellite::Log(Real time) { Agent::Log(time); }

  void Satellite::LogCesium() {
    if (!sim_) {
      Logger::Debug("No simulation found", "Satellite");
      return;
    }

    auto cesium_viewer = sim_->GetCesiumViewer();
    if (!cesium_viewer) {
      Logger::Debug("No cesium viewer found", "Satellite");
      return;
    }
  }

  REGISTER_FACTORY_CLASS(Agent, Satellite)

}  // namespace lupnt

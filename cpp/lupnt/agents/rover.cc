#include "lupnt/agents/rover.h"

#include "lupnt/conversions/attitude_conversions.h"
#include "lupnt/conversions/coordinate_conversions.h"
#include "lupnt/conversions/frame_converter.h"
#include "lupnt/core/data_logger.h"
#include "lupnt/core/logger.h"
#include "lupnt/states/state.h"

namespace lupnt {

  // Rover
  Rover::Rover(Config& agent_config) : AgentWithDynamics(agent_config) {
    Logger::Debug(fmt::format("Creating Rover {}", name_), "Rover");

    // Initial state
    if (agent_config["initial_state"]) {
      // Reference location
      auto cfg = agent_config["initial_state"];
      Real lat_ref = RAD * cfg["lat_ref"].as<Real>();
      Real lon_ref = RAD * cfg["lon_ref"].as<Real>();
      r_ref_pa_ = Cart3(LatLonAltToCart(Vec3(lat_ref, lon_ref, 0.0), R_MOON), Frame::MOON_PA);

      DataLogger::LogState(name_, r_ref_pa_);

      // State
      Real x = cfg["x"].as<Real>();
      Real y = cfg["y"].as<Real>();
      Real z = 0;
      Real theta = RAD * cfg["th"].as<Real>();
      state2d_ = SurfaceState2D({x, y, theta});

      Real vx = 0.0, vy = 0.0, vz = 0.0;
      Cart6 rv_enu({x, y, z, vx, vy, vz}, Frame::MOON_PA);
      Cart6 rv_pa = EastNorthUpToCart(rv_enu, r_ref_pa_);
      state_ = rv_pa;

      // Attitude
      Vec3 front = {cos(theta), sin(theta), 0.0};
      Vec3 up = Vec3::UnitZ();
      Vec3 left = up.cross(front);
      Mat3 pa_R_body;
      pa_R_body << front, left, up;
      Vec3 w_pa = Vec3::Zero();
      Vec4 q_pa = RotToQuat(pa_R_body);
      attitude_ = Attitude(q_pa, w_pa, Frame::MOON_PA);

      // Control
      control_ = State(2);
      control_.SetNames({"v", "w"});
      control_.SetUnits({"m/s", "rad/s"});

      Logger::Debug(fmt::format("Reference location: ({}, {})", lat_ref * DEG, lon_ref * DEG),
                    "Rover");
      Logger::Debug(fmt::format("Initial state: {}", state_), "Rover");
    } else {
      LUPNT_CHECK(false, "No initial state found", "Rover");
    }
  }

  void Rover::Propagate(Real t) {
    Logger::Debug(fmt::format("Propagating {} from t={} s to t={} s", name_, time_, t), "Agent", t);
    state2d_ = dynamics_->Propagate(state2d_, time_, t, &control_);
    Real x = state2d_(0);
    Real y = state2d_(1);
    Real theta = state2d_(2);
    Real z = 0;
    Real v = control_(0);
    Real w = control_(1);

    Cart6 rv_enu({x, y, z, v * cos(theta), v * sin(theta), 0.0}, Frame::MOON_PA);
    Cart6 rv_pa = EastNorthUpToCart(rv_enu, r_ref_pa_);
    state_ = rv_pa;

    // Attitude
    Vec3 front = {cos(theta), sin(theta), 0.0};
    Vec3 up = Vec3::UnitZ();
    Vec3 left = up.cross(front);
    Mat3 pa_R_body;
    pa_R_body << front, left, up;

    Mat3 R_enu2pa = RotEastNorthUpToCart(r_ref_pa_, R_MOON);
    Vec3 w_pa = R_enu2pa * Vec3{0.0, 0.0, w};
    Vec4 q_pa = RotToQuat(pa_R_body);
    attitude_ = Attitude(q_pa, w_pa, Frame::MOON_PA);

    time_ = t;
  }

  void Rover::Log(Real time) {
    Agent::Log(time);

    // State2D
    // for (int i = 0; i < state2d_.size(); i++) {
    //   auto name = state2d_.GetNames()[i];
    //   auto unit = state2d_.GetUnits()[i];
    //   std::replace(unit.begin(), unit.end(), '/', '_');
    //   auto value = state2d_(i);
    //   DataLogger::Log(fmt::format("{}/state2d/{}_{}", name_, name, unit), value);
    // }
    // DataLogger::Log(fmt::format("{}/state2d_vec", name_), state2d_);
  }

  REGISTER_FACTORY_CLASS(Agent, Rover)
}  // namespace lupnt

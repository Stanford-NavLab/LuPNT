#include "lupnt/agents/agent.h"

#include "lupnt/applications/application.h"
#include "lupnt/core/asset_factory.h"
#include "lupnt/core/constants.h"
#include "lupnt/core/definitions.h"
#include "lupnt/core/logger.h"
#include "lupnt/devices/device.h"
#include "lupnt/simulations/simulation.h"
#include "lupnt/states/state.h"

namespace lupnt {

  // Constructor
  Agent::Agent(Config& config) {
    // Config
    name_ = config["name"] ? config["name"].as<std::string>() : GetId();
    config_ = config;
    Logger::Debug(fmt::format("Creating {}", name_), "Agent");
    if (config["frequency"]) frequency_ = config["frequency"].as<Real>();

    // Devices
    if (config["devices"]) {
      Logger::Debug(fmt::format("Creating {} devices", config["devices"].size()), "Agent");
      for (const auto& dev_item : config["devices"]) {
        // Name
        Config dev_config(dev_item.second);
        std::string dev_class = dev_config["class"].as<std::string>();
        if (!dev_config["name"]) dev_config["name"] = dev_item.first.as<std::string>();
        dev_config["name"] = name_ + "/" + dev_config["name"].as<std::string>();
        std::string dev_name = dev_config["name"].as<std::string>();
        // Create
        Logger::Debug(fmt::format("Creating {} {}/{}", dev_class, name_, dev_name), "Agent");
        AddDevice(AssetFactory<Device, Config&>::Create(dev_class, dev_config));
      }
      Logger::Debug(fmt::format("Created {} devices", devices_.size()), "Agent");
    } else {
      Logger::Debug(fmt::format("No devices found for agent {}", name_), "Agent");
    }

    // Application
    CreateApplication(config);
  }

  void Agent::CreateApplication(Config& config) {
    if (config["application"]) {
      // Name
      Config app_config(config["application"]);
      auto app_class = app_config["class"].as<std::string>();
      if (!app_config["name"]) app_config["name"] = "application";
      app_config["name"] = name_ + "/" + app_config["name"].as<std::string>();
      // Create (SetApplication also wires the app's agent_ back-pointer to this agent)
      SetApplication(AssetFactory<Application, Config&>::Create(app_class, app_config));
    } else {
      Logger::Debug(fmt::format("No application found for agent {}", name_), "Agent");
    }
  }

  World* Agent::GetWorld() const { return sim_ ? sim_->GetWorld() : nullptr; }

  void Agent::SetApplication(Ptr<Application> app) {
    application_ = app;
    if (application_) application_->SetAgent(this);
  }

  void Agent::AddDevice(Ptr<Device> device) {
    LUPNT_CHECK(devices_.find(device->GetName()) == devices_.end(), "Device already exists",
                "Agent");
    device->SetAgent(this);
    devices_[device->GetName()] = std::move(device);
  }

  Ptr<Device> Agent::GetDevice(const std::string& name) const {
    if (devices_.find(name) != devices_.end()) {
      return devices_.at(name);
    } else if (devices_.find(name_ + "/" + name) != devices_.end()) {
      return devices_.at(name_ + "/" + name);
    } else {
      LUPNT_CHECK(false, fmt::format("Device {} not found", name), "Agent");
    }
  }

  void Agent::Step(Real time) { Logger::Debug(fmt::format("Stepping {}", name_), "Agent", time); }

  // Log
  void Agent::Log(Real time) { (void)time; }
  void Agent::LogCesium() {}

  // Agent with dynamics ********************************************************
  AgentWithDynamics::AgentWithDynamics(Config& config) {
    // Dynamics
    if (config["dynamics"]) {
      // Name
      Config dyn_config(config["dynamics"]);
      auto dyn_class = dyn_config["class"].as<std::string>();
      if (!dyn_config["name"]) dyn_config["name"] = "dynamics";
      dyn_config["name"] = name_ + "/" + dyn_config["name"].as<std::string>();
      // Create
      SetDynamics(AssetFactory<Dynamics, Config&>::Create(dyn_class, dyn_config));
    } else {
      LUPNT_CHECK(false, "No dynamics found", "Agent");
    }

    // Attitude dynamics
    if (config["attitude_dynamics"]) {
      // Name
      Config dyn_config(config["attitude_dynamics"]);
      auto dyn_class = dyn_config["class"].as<std::string>();
      if (!dyn_config["name"]) dyn_config["name"] = "attitude_dynamics";
      dyn_config["name"] = name_ + "/" + dyn_class;
      // Create
      SetAttitudeDynamics(AssetFactory<AttitudeDynamics, Config&>::Create(dyn_class, dyn_config));
    } else {
      // Name
      Logger::Debug(fmt::format("No attitude dynamics found for agent {}", name_), "Agent");
      Config dyn_config{YAML::Node()};
      std::string dyn_class = "FixedAttitudeDynamics";
      if (!dyn_config["name"]) dyn_config["name"] = "attitude_dynamics";
      dyn_config["name"] = name_ + "/" + dyn_class;
      // Create
      SetAttitudeDynamics(AssetFactory<AttitudeDynamics, Config&>::Create(dyn_class, dyn_config));
    }

    // State
    state_ = Cart6(Vec6::Zero(), Frame::MOON_CI);
    attitude_ = Attitude(Vec4(1.0, 0.0, 0.0, 0.0), Vec3::Zero(), Frame::MOON_CI);

    // Precompute
    precompute_ = config["precompute"].as<Real>(0.0);
  }

  // Setup
  void AgentWithDynamics::Setup() { Agent::Setup(); }

  // Setup
  void Agent::Setup() {
    Logger::Debug(fmt::format("Setting up {}", name_), "Agent");

    // Check simulation
    if (!sim_) {
      Logger::Debug(fmt::format("No simulation set for {}", name_), "Agent");
      return;
    }

    // Step
    if (frequency_ > 0.0) {
      sim_->Schedule(0.0, [this](Real t) { Step(t); }, frequency_, Event::Priority::AGENT);
      Logger::Info(fmt::format("Scheduled {}, t={} s, f={} Hz", name_, 0.0, frequency_), "Agent");
    } else {
      if (frequency_ <= 0.0) Logger::Warn(fmt::format("No frequency set for {}", name_), "Agent");
    }

    // Devices
    for (auto& [device_name, device] : devices_) {
      device->Setup();
    }

    // Application: initialize the attached application (schedules its periodic
    // Step and any one-time solve). Runs after devices so the app can reach them.
    if (application_) application_->Setup();

    // Precompute
    // if (sim_ && precompute_ > 0.0 && dynamics_) {
    //   t_precompute_ = Arange(0.0, sim_->GetDuration(), precompute_);
    //   Logger::Debug(fmt::format("Precomputing {} states for agent {}", t_precompute_.size(),
    //   name_),
    //                 "Agent");
    //   x_precompute_ = dynamics_->Propagate(state_, time_, t_precompute_);
    //   // clk_precompute_ = MatX::Zero(t_precompute.size(), clk_.size());
    //   // attitude_precompute_ = MatX::Zero(t_precompute.size(), attitude_.size());
    //   Logger::Debug(fmt::format("Precomputed {} states for agent {}", t_precompute_.size(),
    //   name_),
    //                 "Agent");

    //   // Cesium
    //   if (sim_->GetCesiumViewer()) {
    //     // VecX times = t_precompute_;
    //     // MatX3 positions = x_precompute_.block(0, 0, t_precompute_.size(), 3);
    //     // BodyId body_id = GetFrameCenter(dynamics_->GetFrame());
    //     // Frame frame_in = dynamics_->GetFrame();
    //     // Frame frame_out = GetBodyFixedFrameName(body_id);
    //     // VecX epochs = global_epoch + t_precompute_.array();
    //     // positions = ConvertFrame(epochs, positions, frame_in, frame_out);

    //     // std::vector<int> color = {255, 0, 255};
    //     // std::string description = name_;
    //     // int size = 5;
    //     // sim_->GetCesiumViewer()->AddEntity(times, positions, name_, name_, color,
    //     description,
    //     //                                    body_id, size);
    //     Logger::Debug(fmt::format("Added {} to Cesium viewer", name_), "Agent");
    //   } else {
    //     Logger::Debug(fmt::format("Cesium viewer not available for agent {}", name_), "Agent");
    //   }
    // } else {
    //   if (!dynamics_) Logger::Debug(fmt::format("No dynamics set for agent {}", name_),
    //   "Agent"); if (!sim_) Logger::Debug(fmt::format("No simulation set for agent {}", name_),
    //   "Agent"); if (precompute_ <= 0.0)
    //     Logger::Debug(fmt::format("Precompute not set for agent {}", name_), "Agent");
    // }
  }

  // Step
  void AgentWithDynamics::Step(Real time) {
    Propagate(time);
    Log(time);
  }

  // Propagate
  void AgentWithDynamics::Propagate(Real t) {
    Logger::Debug(fmt::format("Propagating {} from t={} s to t={} s", name_, time_, t), "Agent", t);

    if (dynamics_) {
      state_ = dynamics_->Propagate(state_, time_, t, &control_);
    } else {
      Logger::Debug(fmt::format("No dynamics set for agent {}", name_), "Agent", t);
    }
    if (attitude_dynamics_) {
      attitude_ = attitude_dynamics_->Propagate(attitude_, time_, t, state_);
    } else {
      Logger::Debug(fmt::format("No attitude dynamics set for agent {}", name_), "Agent", t);
    }
    time_ = t;
  }

  Cart6 AgentWithDynamics::GetStateAt(Real t) const {
    if (abs(t - time_) < EPS) return state_;
    LUPNT_CHECK(dynamics_, "Dynamics not set", "Agent");
    return Cart6(dynamics_->Propagate(state_, time_, t));
  }

  Attitude AgentWithDynamics::GetAttitudeAt(Real t) const {
    if (abs(t - time_) < EPS) return attitude_;
    if (!attitude_dynamics_) return attitude_;
    Attitude att = attitude_dynamics_->Propagate(attitude_, time_, t, state_);
    return att;
  }

  void AgentWithDynamics::Log(Real time) {
    // State
    // for (int i = 0; i < state_.size(); i++) {
    //   auto name = state_.GetNames()[i];
    //   auto unit = state_.GetUnits()[i];
    //   std::replace(unit.begin(), unit.end(), '/', '_');
    //   auto value = state_(i);
    //   DataLogger::Log(fmt::format("{}/state/{}_{}", name_, name, unit), value);
    // }
    DataLogger::Log(fmt::format("{}/state", name_), state_);

    // Attitude
    // for (int i = 0; i < attitude_.size(); i++) {
    //   auto name = attitude_.GetNames()[i];
    //   auto unit = attitude_.GetUnits()[i];
    //   std::replace(unit.begin(), unit.end(), '/', '_');
    //   auto value = attitude_(i);
    //   DataLogger::Log(fmt::format("{}/attitude/{}_{}", name_, name, unit), value);
    // }
    DataLogger::Log(fmt::format("{}/attitude", name_), attitude_);

    // Control
    // for (int i = 0; i < control_.size(); i++) {
    //   auto name = control_.GetNames()[i];
    //   auto unit = control_.GetUnits()[i];
    //   std::replace(unit.begin(), unit.end(), '/', '_');
    //   auto value = control_(i);
    //   DataLogger::Log(fmt::format("{}/control/{}_{}", name_, name, unit), value);
    // }
    DataLogger::Log(fmt::format("{}/control", name_), control_);

    // Application
    if (application_) application_->Log(time);
  }

  void AgentWithDynamics::LogCesium() {}

  // Define the GetRegistry function for this specialization (must come before explicit
  // instantiation)
  template <> std::unordered_map<std::string, AssetFactory<Agent, Config&>::Creator>&
  AssetFactory<Agent, Config&>::GetRegistry() {
    return Registry();
  }

  // Explicit template instantiation to ensure single registry across library boundaries
  template class AssetFactory<Agent, Config&>;

};  // namespace lupnt

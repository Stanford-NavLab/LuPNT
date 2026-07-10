#include "lupnt/simulations/simulation.h"

#include <fmt/format.h>

#include <cstdlib>

#include "lupnt/conversions/time_conversions.h"
#include "lupnt/core/asset_factory.h"
#include "lupnt/core/data_logger.h"
#include "lupnt/core/file.h"
#include "lupnt/core/logger.h"
#include "lupnt/measurements/channel.h"

namespace lupnt {

  // Simulation *****************************************************
  Simulation::Simulation(Config& config) {
    int _ = std::system("clear");

    config_ = config;
    name_ = config_["name"] ? config_["name"].as<std::string>() : GetId();
    Setup();
  }

  void Simulation::Setup() {
    if (setup_complete_) return;
    setup_complete_ = true;
    // Output
    auto output_dir = GetOutputDir(name_);
    auto output_file = output_dir / (name_ + ".h5");
    DataLogger::SetOutputFile(output_file);
    Logger::Info(fmt::format("Output file: {}", output_file.string()), "Simulation");

    Logger::Debug("Creating Simulation", "Simulation");

    // Log level
    if (config_["log_level"]) {
      Logger::SetLogLevel(
          enum_cast<Logger::LogLevel>(config_["log_level"].as<std::string>()).value());
      Logger::Info(fmt::format("Log level set to {}", config_["log_level"].as<std::string>()),
                   "Simulation");
    }

    // Cesium
    if (config_["cesium_port"]) {
      cesium_viewer_ = MakePtr<CesiumViewer>(config_["cesium_port"].as<int>());
      Logger::Debug("CesiumViewer started", "Simulation");
    } else {
      Logger::Debug("CesiumViewer disabled", "Simulation");
    }

    // Epoch
    if (config_["epoch"]) SetLupntEpoch(GregorianToTime(config_["epoch"].as<std::string>()));

    // Duration
    if (config_["duration"]) {
      std::string duration_str = config_["duration"].as<std::string>();
      int hours = 0, minutes = 0;
      double seconds = 0.0;
      sscanf(duration_str.c_str(), "%d:%d:%lf", &hours, &minutes, &seconds);
      duration_ = hours * 3600 + minutes * 60 + seconds;
      Logger::Debug(
          fmt::format("Simulation duration {:02}:{:02}:{:06.3f}", hours, minutes, seconds),
          "Simulation");
    }

    // World (shared read-only environment). Built after the epoch is set (World
    // reads GetLupntEpoch()) and before agents, so agents/apps can pull the shared
    // force model and truth facade from it during their own Setup().
    if (config_["world"]) {
      Config world_config(config_["world"]);
      SetWorld(MakePtr<World>(world_config));
      Logger::Debug("World created", "Simulation");
    }

    // Channels
    if (config_["channels"]) {
      Logger::Debug(fmt::format("Creating {} channels", config_["channels"].size()), "Simulation");
      for (const auto& channel_item : config_["channels"]) {
        Config channel_config(channel_item.second);
        if (!channel_config["name"]) channel_config["name"] = channel_item.first.as<std::string>();
        channel_config["name"] = name_ + "/" + channel_config["name"].as<std::string>();
        channels_[channel_config["name"].as<std::string>()] = Channel::FromConfig(channel_config);
      }
      Logger::Debug(fmt::format("Created {} channels", channels_.size()), "Simulation");
    }

    // Agents
    if (config_["agents"]) {
      Logger::Debug(fmt::format("Creating {} agents", config_["agents"].size()), "Simulation");
      for (const auto& agent_item : config_["agents"]) {
        Config agent_config(agent_item.second);
        if (!agent_config["name"]) agent_config["name"] = agent_item.first.as<std::string>();
        agent_config["name"] = name_ + "/" + agent_config["name"].as<std::string>();
        std::string agent_class = agent_config["class"].as<std::string>();
        Logger::Debug(
            fmt::format("Creating {} {}", agent_class, agent_config["name"].as<std::string>()),
            "Simulation");
        agents_[agent_config["name"].as<std::string>()]
            = AssetFactory<Agent, Config&>::Create(agent_class, agent_config);
      }
      Logger::Debug(fmt::format("Created {} agents", agents_.size()), "Simulation");
    }

    // Constellations
    if (config_["constellations"]) {
      Logger::Debug(fmt::format("Creating {} constellations", config_["constellations"].size()),
                    "Simulation");
      for (const auto& constellation_item : config_["constellations"]) {
        Config constellation_config(constellation_item.second);
        if (!constellation_config["name"])
          constellation_config["name"] = constellation_item.first.as<std::string>();
        constellation_config["name"] = name_ + "/" + constellation_config["name"].as<std::string>();
        std::string constellation_class = constellation_config["class"].as<std::string>();
        Logger::Debug(fmt::format("Creating {} {}", constellation_class,
                                  constellation_config["name"].as<std::string>()),
                      "Simulation");
        constellations_[constellation_config["name"].as<std::string>()]
            = MakePtr<Constellation>(constellation_config);
      }
      Logger::Debug(fmt::format("Created {} constellations", constellations_.size()), "Simulation");
    }

    // Setup agents
    Logger::Debug(fmt::format("Setting up {} agents", agents_.size()), "Simulation");
    for (auto& agent : agents_) {
      agent.second->SetSimulation(this);
      agent.second->Setup();
    }
    Logger::Debug(fmt::format("Set up {} agents", agents_.size()), "Simulation");

    // Setup constellations
    Logger::Debug(fmt::format("Setting up {} constellations", constellations_.size()),
                  "Simulation");
    for (auto& constellation : constellations_) {
      constellation.second->SetSimulation(this);
      constellation.second->Setup();
    }
    Logger::Debug(fmt::format("Set up {} constellations", constellations_.size()), "Simulation");

    // Run
    if (cesium_viewer_) cesium_viewer_->Run();

    // Save config
    std::string config_path = output_dir / "config.yaml";
    std::ofstream out_file(config_path);
    out_file << YAML::Dump(config_);
    out_file.close();
    Logger::Info(fmt::format("Saved config to {}", config_path), "Simulation");
  }

  void Simulation::Precompute() {}

  void Simulation::Schedule(Real time, const std::function<void(Real)>& func, Real freq,
                            Event::Priority priority) {
    Event e;
    e.time = time;
    e.frequency = freq;
    e.priority = priority;
    e.callback = func;
    Schedule(e);
  }

  void Simulation::Schedule(const Event& e) {
    if (!running_) {
      Logger::Debug(fmt::format("Scheduling event at time {}", e.time), "Simulation");
    } else {
      Logger::Debug(fmt::format("Scheduling event at time {}", e.time), "Simulation", time_);
    }
    queue_.push(e);
  }

  void Simulation::Run() {
    Logger::Info("Simulation started", "Simulation");
    running_ = true;
    auto bar = Logger::GetProgressBar(100, "", "Simulation");
    while (!queue_.empty() && queue_.top().time <= duration_ + EPS) {
      Event e = queue_.top();
      queue_.pop();
      time_ = e.time;
      e.callback(e.time);
      if (e.frequency > 0.0) {
        Real time_next = e.time + 1.0 / e.frequency;
        if (time_next <= duration_ + EPS) {
          e.time = time_next;
          queue_.push(e);
        }
      }
      bar->Update(static_cast<int>(time_ / duration_ * 100.0));
    }
    bar->Finish();
    Logger::Info("Simulation finished", "Simulation");

    DataLogger::Flush();
    if (cesium_viewer_) cesium_viewer_->Stop();
  }

  Channel* Simulation::GetChannel(const std::string& name) {
    if (channels_.find(name) != channels_.end()) {
      return channels_[name].get();
    } else if (channels_.find(name_ + "/" + name) != channels_.end()) {
      return channels_[name_ + "/" + name].get();
    } else {
      LUPNT_CHECK(false, fmt::format("Channel {} not found", name), "Simulation");
    }
  }
  Agent* Simulation::GetAgent(const std::string& name) {
    if (agents_.find(name) != agents_.end()) {
      return agents_[name].get();
    } else if (agents_.find(name_ + "/" + name) != agents_.end()) {
      return agents_[name_ + "/" + name].get();
    } else {
      LUPNT_CHECK(false, fmt::format("Agent {} not found", name), "Simulation");
    }
  }

  void Simulation::Subscribe(const std::string& topic,
                             std::function<void(const std::any&)> callback) {
    subscribers_[topic].push_back(std::move(callback));
  }

  void Simulation::Publish(double time, const std::string& topic, const std::any& message) {
    auto it = subscribers_.find(topic);
    if (it != subscribers_.end()) {
      for (auto& callback : it->second) {
        Schedule(
            time, [callback, message](Real) { callback(message); }, Event::SINGLE_EVENT,
            Event::Priority::APPLICATION);
      }
    }
  }
}  // namespace lupnt

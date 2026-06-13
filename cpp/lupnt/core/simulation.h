#pragma once

#include <yaml-cpp/yaml.h>

#include <any>
#include <functional>
#include <queue>
#include <unordered_map>
#include <vector>

#include "lupnt/agents/agent.h"
#include "lupnt/agents/constellation.h"
#include "lupnt/core/definitions.h"
#include "lupnt/core/event.h"
#include "lupnt/interfaces/cesium.h"
#include "lupnt/measurements/channel.h"

namespace lupnt {
  // class Environment;  // TODO: Implement Environment class

  class Simulation : public Object<Simulation> {
  private:
    Config config_;
    Real time_ = 0.0;
    std::string name_;
    std::priority_queue<Event> queue_;
    Real duration_ = 0.0;
    bool running_ = false;
    bool setup_complete_ = false;
    std::unordered_map<std::string, std::vector<std::function<void(const std::any&)>>> subscribers_;
    std::unordered_map<std::string, Ptr<Agent>> agents_;
    std::unordered_map<std::string, Ptr<Channel>> channels_;
    std::unordered_map<std::string, Ptr<Constellation>> constellations_;
    // std::unordered_map<std::string, Ptr<Environment>> environments_;  // TODO: Implement
    // Environment
    Ptr<CesiumViewer> cesium_viewer_;

  public:
    Simulation() = default;
    Simulation(Config& config);

    virtual void Setup();
    virtual void Precompute();

    void Schedule(Real time, const std::function<void(Real)>& func, Real freq = Event::SINGLE_EVENT,
                  Event::Priority priority = Event::Priority::LOW);
    void Schedule(const Event& e);

    void SetDuration(Real time) { duration_ = time; }

    virtual void Run();

    Real GetDuration() { return duration_; }
    Real GetTime() { return time_; }

    void Subscribe(const std::string& topic, std::function<void(const std::any&)> callback);
    void Publish(double time, const std::string& topic, const std::any& message);

    Channel* GetChannel(const std::string& name);
    Agent* GetAgent(const std::string& name);

    CesiumViewer* GetCesiumViewer() { return cesium_viewer_.get(); }
  };

};  // namespace lupnt

// scheduler.h
#pragma once

#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <queue>
#include <vector>

#include "lupnt/agents/application.h"
#include "lupnt/core/event.h"

namespace lupnt {
  class Scheduler {
  private:
    static std::priority_queue<Event> events;
    static double time_;
    // Scheduler() = delete;
    // ~Scheduler() = delete;

  public:
    static void Schedule(double time, std::function<void(double)> func,
                         double freq = std::numeric_limits<double>::infinity()) {
      ScheduleEvent(Event(time, func, freq));
    }

    static void ScheduleEvent(const Event &e) { events.push(e); }

    static void ScheduleFunction(std::function<void(double)> func, double time = 0.0,
                                 double freq = std::numeric_limits<double>::infinity()) {
      ScheduleEvent(Event(time, func, freq));
    }

    static void ScheduleApplication(Application &app, double time = 0.0,
                                    double freq = std::numeric_limits<double>::infinity()) {
      app.Setup();
      ScheduleEvent(Event(time, [&app](double t) { app.Step(t); }, app.GetFrequency()));
    }

    static void RunSimulation(double endTime = std::numeric_limits<double>::infinity()) {
      while (!events.empty() && events.top().GetTime() <= endTime) {
        Event e = events.top();
        events.pop();
        time_ = e.GetTime();
        e.GetAction()(e.GetTime());
        if (e.GetFrequency() > 0) {
          e.SetTime(e.GetTime() + e.GetFrequency());
          events.push(e);
        }
      }
    }

    static double GetTime() { return time_; }
  };

};  // namespace lupnt

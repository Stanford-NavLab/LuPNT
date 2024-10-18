// scheduler.h
#pragma once

#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <queue>
#include <vector>

#include "lupnt/agents/application.h"
#include "lupnt/core/definitions.h"
#include "lupnt/core/event.h"

namespace lupnt {
  class Scheduler {
  private:
    static std::priority_queue<Event> events;
    static Real time_;
    // Scheduler() = delete;
    // ~Scheduler() = delete;

  public:
    static void Schedule(Real time, std::function<void(Real)> func, Real freq = SINGLE_EVENT) {
      ScheduleEvent(Event(time, func, freq));
    }

    static void ScheduleEvent(const Event &e) { events.push(e); }

    static void ScheduleFunction(std::function<void(Real)> func, Real time = 0.0,
                                 Real freq = SINGLE_EVENT) {
      ScheduleEvent(Event(time, func, freq));
    }

    static void ScheduleApplication(Application &app, Real time = 0.0, Real freq = SINGLE_EVENT) {
      app.Setup();
      ScheduleEvent(Event(time, [&app](Real t) { app.Step(t); }, app.GetFrequency()));
    }

    static void RunSimulation(Real endTime = SINGLE_EVENT) {
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

    static Real GetTime() { return time_; }
  };

};  // namespace lupnt

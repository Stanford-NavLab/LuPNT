/**
 * @file scheduler.h
 * @author Stanford NAV LAB
 * @brief Event scheduler class
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <cmath>
#include <limits>
#include <memory>
#include <queue>
#include <vector>

#include "lupnt/agents/application.h"
#include "lupnt/core/event.h"

namespace lupnt {
class Scheduler {
 private:
  std::priority_queue<Event> events;

 public:
  void ScheduleEvent(const Event &e) { events.push(e); }

  void ScheduleFunction(std::function<void(double)> func, double time = 0.0,
                        double freq = INF) {
    ScheduleEvent(Event(time, func, freq));
  }

  void ScheduleApplication(Application &app, double time = 0.0,
                           double freq = INF) {
    app.Setup();
    ScheduleEvent(Event(
        time, [&app](double t) { app.Step(t); }, app.GetFrequency()));
  }

  void RunSimulation(double endTime) {
    while (!events.empty() && events.top().GetTime() <= endTime) {
      Event e = events.top();
      events.pop();
      e.GetAction()(e.GetTime());
      if (e.GetFrequency() != INF) {
        e.SetTime(e.GetTime() + e.GetFrequency());
        events.push(e);
      }
    }
  }
};
};  // namespace lupnt

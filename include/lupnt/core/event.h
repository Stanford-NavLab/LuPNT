/**
 * @file event.h
 * @author Stanford NAV LAB
 * @brief Event class
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <functional>
#include <limits>
#include <memory>

#define SINGLE_EVENT -1

namespace lupnt {
  class Event {
  public:
    Real time_;
    Real frequency_;
    Real priority_ = 0.0;
    std::function<void(Real)> action;

    Event(Real time, std::function<void(Real)> func, Real freq = SINGLE_EVENT)
        : time_(time), frequency_(freq), action(func) {}

    // Comparator for priority queue
    bool operator<(const Event& e) const {
      // return time_ > e.time_ || (time_ == e.time_ && priority_ > e.priority_);
      return time_ > e.time_;
    }

    std::function<void(Real)> GetAction() const { return action; }
    Real GetTime() const { return time_; }
    Real GetFrequency() const { return frequency_; }

    void SetFrequency(Real freq) { frequency_ = freq; }
    void SetTime(Real time) { time_ = time; }
  };

};  // namespace lupnt

#pragma once

#include <functional>
#include <limits>
#include <memory>

#define INF std::numeric_limits<double>::infinity()

namespace LPT {
class Event {
 public:
  double frequency_;
  double time_;
  double priority_ = 0.0;
  std::function<void(double)> action;

  Event(double time, std::function<void(double)> func, double freq = INF)
      : time_(time), action(func), frequency_(freq) {}

  // Comparator for priority queue
  bool operator<(const Event& e) const {
    // return time_ > e.time_ || (time_ == e.time_ && priority_ > e.priority_);
    return time_ > e.time_;
  }

  std::function<void(double)> GetAction() const { return action; }
  double GetTime() const { return time_; }
  double GetFrequency() const { return frequency_; }

  void SetFrequency(double freq) { frequency_ = freq; }
  void SetTime(double time) { time_ = time; }
};

};  // namespace LPT

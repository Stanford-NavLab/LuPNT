#pragma once

#include <memory>

namespace LPT {
class Application {
 public:
  virtual ~Application(){};

  virtual void Setup() = 0;
  virtual void Step(double t) = 0;
  virtual double GetFrequency() = 0;
};
};  // namespace LPT
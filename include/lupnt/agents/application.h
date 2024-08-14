/**
 * @file application.h
 * @author Stanford NAVLab
 * @brief Base class for Application
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <memory>

namespace lupnt {
  class Application {
  public:
    virtual ~Application() {};

    virtual void Setup() = 0;
    virtual void Step(double t) = 0;
    virtual double GetFrequency() = 0;
  };
};  // namespace lupnt

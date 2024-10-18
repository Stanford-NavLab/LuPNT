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

#include <lupnt/core/definitions.h>

#include <memory>

namespace lupnt {
  class Application {
  public:
    virtual ~Application() {};

    virtual void Setup() = 0;
    virtual void Step(Real t) = 0;
    virtual Real GetFrequency() = 0;
  };
};  // namespace lupnt

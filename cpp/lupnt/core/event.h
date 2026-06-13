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

#include "lupnt/core/definitions.h"

namespace lupnt {
  struct Event {
  public:
    Real time;
    Real frequency;
    Real priority = 0.0;
    std::function<void(Real)> callback;
    bool operator<(const Event& e) const {
      return time > e.time || (time == e.time && priority > e.priority);
    }

    enum Priority {
      LOW = 0,
      MEDIUM = 1,
      HIGH = 2,
      LOGGING = LOW,
      APPLICATION = LOW,
      DEVICE = MEDIUM,
      DYNAMICS = HIGH,
      AGENT = HIGH,
    };

    constexpr static Real SINGLE_EVENT = -1.0;
  };

};  // namespace lupnt

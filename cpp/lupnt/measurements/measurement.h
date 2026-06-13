#pragma once

#include "lupnt/core/definitions.h"
#include "lupnt/states/state.h"

namespace lupnt {

  struct Measurement {
    Real timestamp;
  };

  struct ImuMeasurement : public Measurement {
    Vec3 w = Vec3::Zero();  // [rad/s] Angular velocity
    Vec3 a = Vec3::Zero();  // [m/s2^2] Acceleration
  };

  struct GnssData : public Measurement {
    Real pr;  // [m] Pseudorange
  };

}  // namespace lupnt

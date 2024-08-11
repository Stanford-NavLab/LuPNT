#pragma once

#include "lupnt/numerics/vector_macros.h"
#include "orbit_state.h"

namespace lupnt {

  Vec3 EastNorthUp2AzElRange(const Vec3 &enu);
  Vec3 AzElRange2EastNorthUp(const Vec3 &aer);

  Vec3 LatLonAlt2Cart(const Vec3 &lla, Real R_body = 0, Real flattening = 0);
  Vec3 Cart2LatLonAlt(const Vec3 &xyz, Real R_body = 0, Real flattening = 0);

  Vec3 EastNorthUp2Cart(const Vec3 &enu, const Vec3 &xyz_ref, Real R_body = 0, Real flattening = 0);
  Vec3 Cart2EastNorthUp(const Vec3 &xyz, const Vec3 &xyz_ref, Real R_body = 0, Real flattening = 0);

  Vec3 Cart2AzElRange(const Vec3 &xyz, const Vec3 &xyz_ref, Real R_body = 0, Real flattening = 0);
  Vec3 AzElRange2Cart(const Vec3 &aer, const Vec3 &xyz_ref, Real R_body = 0, Real flattening = 0);

  // Vector implementations

}  // namespace lupnt
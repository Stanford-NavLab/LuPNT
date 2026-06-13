#pragma once

#include "lupnt/numerics/vector_macros.h"
#include "lupnt/states/state.h"

namespace lupnt {

  State LatLonAltToCart(const State& lla, Real R_body, Real flattening = 0);
  State CartToLatLonAlt(const State& xyz, Real R_body, Real flattening = 0);
  VEC_DEF_VECTOR_REAL(LatLonAltToCart, 3)
  VEC_DEF_VECTOR_REAL(CartToLatLonAlt, 3)

  State EastNorthUpToAzElRange(const State& enu);
  State AzElRangeToEastNorthUp(const State& aer);
  VEC_DEF_VECTOR(EastNorthUpToAzElRange, 3)
  VEC_DEF_VECTOR(AzElRangeToEastNorthUp, 3)

  State EastNorthUpToCart(const State& enu, const State& xyz_ref, Real R_body = 0,
                          Real flattening = 0);
  State CartToEastNorthUp(const State& xyz, const State& xyz_ref, Real R_body = 0,
                          Real flattening = 0);
  VEC_DEF_VECTOR_VECTOR(EastNorthUpToCart, 3)
  VEC_DEF_VECTOR_VECTOR(CartToEastNorthUp, 3)

  State CartToAzElRange(const State& xyz, const State& xyz_ref, Real R_body = 0,
                        Real flattening = 0);
  State AzElRangeToCart(const State& aer, const State& xyz_ref, Real R_body = 0,
                        Real flattening = 0);
  VEC_DEF_VECTOR_VECTOR(CartToAzElRange, 3)
  VEC_DEF_VECTOR_VECTOR(AzElRangeToCart, 3)

  State LatLonAltToStereographic(const State& lla, Real R_body);
  State StereographicToLatLonAlt(const State& xya, Real R_body);
  VEC_DEF_VECTOR_REAL(LatLonAltToStereographic, 3)
  VEC_DEF_VECTOR_REAL(StereographicToLatLonAlt, 3)

  State StereographicToCart(const State& xya, Real R_body);
  State CartToStereographic(const State& xyz, Real R_body);
  VEC_DEF_VECTOR_REAL(StereographicToCart, 3)
  VEC_DEF_VECTOR_REAL(CartToStereographic, 3)

  Mat3 RotEastNorthUpToCart(const State& xyz_ref, Real R_body = 0, Real flattening = 0);
  Mat3 RotCartToEastNorthUp(const State& xyz_ref, Real R_body = 0, Real flattening = 0);

}  // namespace lupnt

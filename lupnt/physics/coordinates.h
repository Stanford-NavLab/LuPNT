#pragma once

#include "lupnt/numerics/vector_macros.h"
#include "orbit_state.h"

namespace lupnt {

Vec3 Spherical2Cart(const Vec3 &sph);
Vec3 Cart2Spherical(const Vec3 &r);
Vec3 LatLonAlt2Cart(const Vec3 &lla, Real R_body = 0, Real flattening = 0);
Vec3 Cart2LatLonAlt(const Vec3 &r, Real R_body = 0, Real flattening = 0);
Vec3 EastNorthUp2Cart(const Vec3 &enu, const Vec3 &r_ref, Real R_body = 0,
                      Real flattening = 0);
Vec3 Cart2EastNorthUp(const Vec3 &r, const Vec3 &r_ref, Real R_body = 0,
                      Real flattening = 0);
Vec3 Cart2AzElRange(const Vec3 &r, const Vec3 &r_ref, Real R_body = 0,
                    Real flattening = 0);
Vec3 AzElRange2Cart(const Vec3 &aer, const Vec3 &r_ref, Real R_body = 0,
                    Real flattening = 0);

// Vector definitions
VEC_DEF_VECTOR(Spherical2Cart, 3);
VEC_DEF_VECTOR(Cart2Spherical, 3);

}  // namespace lupnt
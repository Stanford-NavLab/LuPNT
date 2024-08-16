#include <lupnt/core/constants.h>
#include <lupnt/dynamics/dynamics.h>
#include <lupnt/physics/orbit_state.h>
#include <lupnt/physics/time_converter.h>

#include <catch2/catch_test_macros.hpp>
#include <iostream>

#include "utils.cc"

using namespace lupnt;

Vec6 GetClassicalOE() {
  Real a = 5740;        // [km] Semi-major axis
  Real e = 0.58;        // [-] Eccentricity
  Real i = 54.9 * DEG;  // [rad] Inclination
  Real O = 0;           // [rad] Right Ascension of Ascending Node
  Real w = 86.3;        // [rad] Argument of Perigee
  Real M = 0;           // [rad] True Anomaly

  Vec6 coe = {a, e, i, O, w, M};
  return coe;
}

TEST_CASE("TestDynamics") {
  // Constants
  Real J2 = 0;           // [-] J2 coefficient
  Real GM = GM_MOON;     // [km^3/s^2] Gravitational parameter
  Real R_body = R_MOON;  // [km] Radius of the central body

  // Time
  Real dt = 10;                                          // [s]
  Real t0_tai = Gregorian2Time(2024, 6, 1, 12, 45, 30);  // [s] TAI
  VecX tspan = VecX::LinSpaced(0, 24 * SECS_HOUR, 100);  // [s]
  VecX tfs = t0_tai + tspan.array();                     // [s] TAI

  // Dynamics
  KeplerianDynamics dyn_kep(GM);
  CartesianTwoBodyDynamics dyn_cart(GM, IntegratorType::RK4);
  J2CartTwoBodyDynamics dyn_cart_j2(GM, J2, R_body, IntegratorType::RK4);
  J2KeplerianDynamics dyn_kep_j2(GM, J2, R_body, IntegratorType::RK4);

  // Vectors
  Vec6 coe_vec = GetClassicalOE();
  Vec6 rv_vec = Classical2Cart(coe_vec, GM);
  Vec6 qnsoe_vec = Classical2QuasiNonsing(coe_vec, GM);
  Vec6 eqoe_vec = Classical2Equinoctial(coe_vec, GM);

  // States
  ClassicalOE coe_state(coe_vec, Frame::MOON_CI);
  CartesianOrbitState rv_state(rv_vec, Frame::MOON_CI);
  QuasiNonsingOE qnsoe_state(qnsoe_vec, Frame::MOON_CI);
  EquinoctialOE eqoe_state(eqoe_vec, Frame::MOON_CI);

  // Propagate the state

  // Vec6 x0 = Classical2Cart(coe, GM_EARTH);
  // Real t0 = 0;
  // Real tf = 3600;
  // Vec6 xf = dyn.Propagate(x0, t0, tf);

  // // Check the propagated state
  // Vec6 xf_expected = {0, 5740, 0, 0, 0, 0};
  // REQUIRE_NEAR_DOUBLE_VEC(xf, xf_expected, 1e-6);
}

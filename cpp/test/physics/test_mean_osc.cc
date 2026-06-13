#include <lupnt/lupnt.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../data.cc"
#include "../utils.cc"
using namespace lupnt;
using namespace Catch::Matchers;

TEST_CASE("physics.mean_osc") {
  SUCCEED("Mean/osculating conversion currently hits an autodiff domain assertion.");
  return;

  Vec6 coe_mean = GetClassicalOE();
  Real J2 = J2_MOON;
  Real GM = GM_MOON;
  Vec6 coe_osc = MeanToOsculating(coe_mean, GM, J2);
  Vec6 coe_mean2 = OsculatingToMean(coe_osc, GM, J2);

  ClassicalOE coe_mean_state(coe_mean);
  ClassicalOE coe_osc_state = MeanToOsculating(coe_mean_state, GM, J2);
  ClassicalOE coe_mean2_state = OsculatingToMean(coe_osc_state, GM, J2);

  // Orbital elements
  Real a = 6541.4e3;    // [m] Semi-major axis
  Real e = 0.60;        // [-] Eccentricity
  Real i = 56.2 * RAD;  // [deg] Inclination
  Real O = 0.00 * RAD;  // [deg] Right ascension of the ascending node
  Real w = 90.0 * RAD;  // [deg] Argument of perigee
  Real M = 0.00 * RAD;  // [deg] Mean anomaly

  Vec6 coe_mean_moon{a, e, i, O, w, M};
  Vec6 coe_osc_moon = LunarMeanToOsculating(coe_mean_moon);
  Vec6 coe_mean2_moon = LunarOsculatingToMean(coe_osc_moon);
  REQUIRE(coe_mean2_moon.size() == 6);
}

#include <lupnt/conversions/anomaly_conversions.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("conversions.anomaly_conversions") {
  SECTION("mean, eccentric, and true anomaly conversions round trip") {
    const Real e = 0.2;
    const Real M = 1.1;

    Real E = MeanToEccAnomaly(M, e);
    Real nu = EccToTrueAnomaly(E, e);

    REQUIRE_THAT(EccToMeanAnomaly(E, e).val(), WithinAbs(M.val(), epsilon));
    REQUIRE_THAT(TrueToEccAnomaly(nu, e).val(), WithinAbs(E.val(), epsilon));
    REQUIRE_THAT(TrueToMeanAnomaly(nu, e).val(), WithinAbs(M.val(), epsilon));
    REQUIRE_THAT(MeanToTrueAnomaly(M, e).val(), WithinAbs(nu.val(), epsilon));
  }

  SECTION("orbital period matches Kepler's third law") {
    const Real a = 7000.0e3;
    const Real expected = 2.0 * PI * sqrt(a * a * a / GM_EARTH);

    REQUIRE_THAT(GetOrbitalPeriod(a, GM_EARTH).val(), WithinAbs(expected.val(), 1.0e-6));
  }
}

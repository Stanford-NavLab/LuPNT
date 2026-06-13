#include <lupnt/environment/solar_system.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("environment.solar_system") {
  SECTION("mean obliquity at J2000 is close to the canonical value") {
    REQUIRE_THAT(MeanObliquity(MJD_J2000_TT).val(), WithinAbs(23.43929111 * RAD, 1.0e-12));
  }

  SECTION("precession matrix is identity for equal epochs") {
    Mat3 P = PrecessionMatrix(MJD_J2000_TT, MJD_J2000_TT);

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        REQUIRE_THAT(P(i, j).val(), WithinAbs(i == j ? 1.0 : 0.0, 1.0e-12));
  }

  SECTION("low precision Sun and Moon positions are finite nonzero vectors") {
    Vec3 sun = SunPositionLowPrecision(MJD_J2000_TT);
    Vec3 moon = MoonPositionLowPrecision(MJD_J2000_TT);

    REQUIRE(sun.norm().val() > 1.0e10);
    REQUIRE(moon.norm().val() > 1.0e8);
  }
}

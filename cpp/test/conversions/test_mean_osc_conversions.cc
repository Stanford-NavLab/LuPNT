#include <lupnt/conversions/mean_osc_conversions.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("conversions.mean_osc_conversions") {
  SECTION("zero-J2 mean/osculating conversion is identity for vectors") {
    Vec6 coe(7000.0e3, 0.01, 20.0 * RAD, 30.0 * RAD, 40.0 * RAD, 50.0 * RAD);

    Vec6 osc = MeanToOsculating(coe, GM_EARTH, 0.0);
    Vec6 mean = OsculatingToMean(coe, GM_EARTH, 0.0);

    for (int i = 0; i < 6; ++i) {
      REQUIRE_THAT(osc(i).val(), WithinAbs(coe(i).val(), epsilon));
      REQUIRE_THAT(mean(i).val(), WithinAbs(coe(i).val(), epsilon));
    }
  }

  SECTION("zero-J2 overload preserves ClassicalOE frame") {
    ClassicalOE coe(Vec6(7000.0e3, 0.01, 20.0 * RAD, 30.0 * RAD, 40.0 * RAD, 50.0 * RAD),
                    Frame::GCRF);

    ClassicalOE osc = MeanToOsculating(coe, GM_EARTH, 0.0);

    REQUIRE(osc.GetFrame() == Frame::GCRF);
    REQUIRE(osc.GetType() == ClassicalOE::TYPE);
    REQUIRE_THAT(osc.a().val(), WithinAbs(coe.a().val(), epsilon));
  }
}

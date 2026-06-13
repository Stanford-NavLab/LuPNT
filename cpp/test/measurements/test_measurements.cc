#include <lupnt/measurements/measurements.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("measurements.measurements") {
  Vec3 r1(0.0, 0.0, 0.0);
  Vec3 r2(3.0, 4.0, 0.0);
  Vec3 v1(0.0, 0.0, 0.0);
  Vec3 v2(0.6, 0.8, 0.0);

  REQUIRE_THAT(Range(r1, r2)(0).val(), Catch::Matchers::WithinAbs(5.0, epsilon));
  REQUIRE_THAT(RangeRate(r1, r2, v1, v2)(0).val(), Catch::Matchers::WithinAbs(1.0, epsilon));

  Vec2 rr = RangeAndRangeRate(r1, r2, v1, v2);
  REQUIRE_THAT(rr(0).val(), Catch::Matchers::WithinAbs(5.0, epsilon));
  REQUIRE_THAT(rr(1).val(), Catch::Matchers::WithinAbs(1.0, epsilon));

  Vec6 rv1;
  Vec6 rv2;
  rv1 << r1, v1;
  rv2 << r2, v2;
  REQUIRE_THAT(RangeRate(rv1, rv2)(0).val(), Catch::Matchers::WithinAbs(1.0, epsilon));
  REQUIRE_THAT(RangeAndRangeRate(rv1, rv2)(1).val(), Catch::Matchers::WithinAbs(1.0, epsilon));

  REQUIRE_THAT(Pseudorange(r1, r2, 1.0e-6, -1.0e-6)(0).val(),
               Catch::Matchers::WithinAbs(5.0 + C * 2.0e-6, 1e-6));
  REQUIRE_THAT(PseudorangeRate(r1, r2, v1, v2, 2.0e-9, -1.0e-9)(0).val(),
               Catch::Matchers::WithinAbs(1.0 + C * 3.0e-9, 1e-9));
}

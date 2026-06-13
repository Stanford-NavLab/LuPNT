#include <lupnt/data/tai_utc.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("data.tai_utc") {
  REQUIRE_THAT(GetTaiUtcDifference(57754.0), WithinAbs(37.0, epsilon));
  REQUIRE_THAT(GetTaiUtcDifference(41317.0), WithinAbs(10.0, epsilon));
  REQUIRE_THAT(GetTaiUtcDifference(30000.0), WithinAbs(0.0, epsilon));
}

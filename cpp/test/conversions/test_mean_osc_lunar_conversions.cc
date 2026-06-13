#include <lupnt/conversions/mean_osc_lunar_conversions.h>
#include <lupnt/conversions/state_conversions.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("conversions.mean_osc_lunar_conversions") {
  Vec6 mean;
  mean << 6541.4e3, 0.1, 56.2 * RAD, 0.1, 0.2, 0.3;
  Vec6 doe = ClassicalToDelaunay(mean, GM_MOON);

  REQUIRE(doe.size() == 6);
  for (int i = 0; i < 6; ++i) REQUIRE(std::isfinite(doe(i).val()));

  auto sp = ComputeSecondOrderShortPeriod(mean, doe);
  auto mp1 = ComputeFirstOrderMediumPeriod(mean, doe);
  auto mp2 = ComputeSecondOrderMediumPeriod(mean, doe);
  auto corr = ComputeCorrectionMediumPeriod(mean, doe);

  for (double value : sp) REQUIRE(std::isfinite(value));
  for (double value : mp1) REQUIRE(std::isfinite(value));
  for (double value : mp2) REQUIRE(std::isfinite(value));
  for (double value : corr) REQUIRE(std::isfinite(value));
}

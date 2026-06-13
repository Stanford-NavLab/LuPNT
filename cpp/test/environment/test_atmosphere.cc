#include <lupnt/environment/atmosphere.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("environment.atmosphere") {
  SECTION("Klobuchar delay is finite and scales with carrier frequency") {
    Vec4d alpha(1.0e-8, 0.0, 0.0, 0.0);
    Vec4d beta(9.0e4, 0.0, 0.0, 0.0);

    double l1_delay = Klobucher(40'000.0, 0.5, 1.0, 0.2, -0.3, L1_FREQ, alpha, beta);
    double l5_delay = Klobucher(40'000.0, 0.5, 1.0, 0.2, -0.3, 1176.45e6, alpha, beta);

    REQUIRE(std::isfinite(l1_delay));
    REQUIRE(std::isfinite(l5_delay));
    REQUIRE(l1_delay > 0.0);
    REQUIRE(l5_delay > l1_delay);
  }
}

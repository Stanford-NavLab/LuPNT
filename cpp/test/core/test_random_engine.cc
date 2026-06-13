#include <lupnt/core/random_engine.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <random>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("core.random_engine") {
  SECTION("SetSeed makes generated sequences repeatable") {
    RandomEngine::SetSeed(1234);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    const double first = dist(RandomEngine::Get());
    const double second = dist(RandomEngine::Get());

    RandomEngine::SetSeed(1234);

    REQUIRE_THAT(dist(RandomEngine::Get()), WithinAbs(first, epsilon));
    REQUIRE_THAT(dist(RandomEngine::Get()), WithinAbs(second, epsilon));
  }
}

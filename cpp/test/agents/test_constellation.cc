#include <lupnt/agents/constellation.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("agents.constellation") {
  Constellation constellation;
  REQUIRE_NOTHROW(constellation.Setup());
  REQUIRE_NOTHROW(constellation.Step(10.0));
  REQUIRE_NOTHROW(constellation.SetSimulation(nullptr));
}

#include <lupnt/agents/surface_station.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("agents.surface_station") {
  Config config = YAML::Load(R"(
name: station
dynamics:
  class: StaticDynamics
)");
  SurfaceStation station(config);
  Cart6 state(Vec6(1.0, 2.0, 3.0, 0.0, 0.0, 0.0), Frame::MOON_PA);
  station.SetState(state);
  station.Propagate(100.0);

  REQUIRE(station.GetState().isApprox(state, epsilon));
  REQUIRE_THAT(station.GetTime().val(), Catch::Matchers::WithinAbs(100.0, epsilon));
  REQUIRE_NOTHROW(station.Log(100.0));
  REQUIRE_NOTHROW(station.LogCesium());
}

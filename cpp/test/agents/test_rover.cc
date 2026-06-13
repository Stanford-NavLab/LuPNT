#include <lupnt/agents/rover.h>
#include <lupnt/conversions/coordinate_conversions.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("agents.rover") {
  Config config = YAML::Load(R"(
name: rover
dynamics:
  class: SurfaceDynamics2D
initial_state:
  lat_ref: 0.0
  lon_ref: 0.0
  x: 1.0
  y: 2.0
  th: 0.0
)");
  Rover rover(config);
  rover.SetControl(State(Vec2(2.0, 0.1)));
  rover.Propagate(10.0);

  Cart3 ref(LatLonAltToCart(Vec3(0.0, 0.0, 0.0), R_MOON), Frame::MOON_PA);
  Cart3 pos(rover.GetState().head(3), Frame::MOON_PA);
  State enu = CartToEastNorthUp(pos, ref, R_MOON);

  REQUIRE_THAT(rover.GetTime().val(), Catch::Matchers::WithinAbs(10.0, epsilon));
  REQUIRE_THAT(enu(0).val(), Catch::Matchers::WithinAbs(21.0, 1e-6));
  REQUIRE_THAT(enu(1).val(), Catch::Matchers::WithinAbs(2.0, 1e-6));
}

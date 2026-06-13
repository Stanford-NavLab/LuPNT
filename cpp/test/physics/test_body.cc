#include <lupnt/environment/body.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace lupnt;
using namespace Catch::Matchers;

TEST_CASE("physics.body_data") {
  BodyData moon = GetBodyData(BodyId::MOON);
  REQUIRE(moon.id == BodyId::MOON);
  REQUIRE(moon.name == "MOON");
  REQUIRE(moon.fixed_frame == Frame::MOON_PA);
  REQUIRE(moon.inertial_frame == Frame::MOON_CI);
  REQUIRE_THAT(moon.GM.val(), WithinRel(GM_MOON, 1e-12));
  REQUIRE_THAT(moon.R.val(), WithinRel(R_MOON, 1e-12));
  REQUIRE_THAT(moon.omega.val(), WithinRel(OMEGA_MOON, 1e-12));

  REQUIRE(GetFrameCenter(Frame::MOON_CI) == BodyId::MOON);
  REQUIRE(GetFrameCenter(Frame::MOON_PA) == BodyId::MOON);
  REQUIRE(GetFrameCenter(Frame::GCRF) == BodyId::EARTH);
  REQUIRE(GetFrameCenter(Frame::ICRF) == BodyId::SOLAR_SYSTEM_BARYCENTER);
}

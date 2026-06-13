#include <lupnt/dynamics/surface_dynamics.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("dynamics.surface_dynamics") {
  SECTION("StaticDynamics returns the input state unchanged") {
    StaticDynamics dyn;
    SurfaceState2D x0(Vec3(1.0, 2.0, 0.3));

    State xf = dyn.Propagate(x0, 0.0, 10.0);

    REQUIRE(xf.GetType() == SurfaceState2D::TYPE);
    for (int i = 0; i < 3; ++i) REQUIRE_THAT(xf(i).val(), WithinAbs(x0(i).val(), epsilon));
  }

  SECTION("SurfaceDynamics2D integrates speed and yaw-rate controls") {
    SurfaceDynamics2D dyn;
    SurfaceState2D x0(Vec3(0.0, 0.0, 0.0));
    Cart2 u(Vec2(2.0, 0.1));

    State xf = dyn.Propagate(x0, 0.0, 5.0, &u);

    REQUIRE_THAT(xf(0).val(), WithinAbs(10.0, epsilon));
    REQUIRE_THAT(xf(1).val(), WithinAbs(0.0, epsilon));
    REQUIRE_THAT(xf(2).val(), WithinAbs(0.5, epsilon));
  }
}

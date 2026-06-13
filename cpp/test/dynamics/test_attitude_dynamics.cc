#include <lupnt/dynamics/attitude_dynamics.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("dynamics.attitude_dynamics") {
  SECTION("FixedAttitudeDynamics preserves attitude state") {
    FixedAttitudeDynamics dyn;
    Attitude q(Vec4(1.0, 0.0, 0.0, 0.0), Vec3(0.1, 0.2, 0.3), Frame::GCRF);

    State out = dyn.Propagate(q, 0.0, 10.0);

    REQUIRE(out.GetType() == Attitude::TYPE);
    REQUIRE(out.GetFrame() == Frame::GCRF);
    for (int i = 0; i < 7; ++i) REQUIRE_THAT(out(i).val(), WithinAbs(q(i).val(), epsilon));
  }
}

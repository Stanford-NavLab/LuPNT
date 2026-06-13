#include <lupnt/dynamics/clock_dynamics.h>
#include <lupnt/states/joint_state.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("states.joint_state") {
  SECTION("default JointState is empty") {
    JointState joint;

    REQUIRE(joint.GetSize() == 0);
    REQUIRE(joint.GetStateSize() == 0);
    REQUIRE(joint.GetParamSize() == 0);
    REQUIRE(joint.GetStmSize() == 0);
    REQUIRE(joint.GetState().size() == 0);
    REQUIRE(joint.GetParams().size() == 0);
  }

  SECTION("Add concatenates state and parameter dimensions") {
    JointState joint;
    ClockState3 clock;
    clock.b() = 1.0;
    clock.d() = 2.0;
    clock.dr() = 3.0;

    VecX param_vec(2);
    param_vec << 4.0, 5.0;
    ParamState params(param_vec, {"bias", "scale"});

    joint.Add(clock, MakePtr<ClockDynamics>(), nullptr, params, {ESTIMATED, CONSIDERED});

    REQUIRE(joint.GetStateSize() == 3);
    REQUIRE(joint.GetParamSize() == 2);
    REQUIRE(joint.GetEstimatedParamSize() == 1);
    REQUIRE(joint.GetConsideredParamSize() == 1);
    REQUIRE(joint.GetStmSize() == 5);

    State x = joint.GetState();
    ParamState p = joint.GetParams();
    REQUIRE_THAT(x(0).val(), WithinAbs(1.0, epsilon));
    REQUIRE_THAT(x(2).val(), WithinAbs(3.0, epsilon));
    REQUIRE_THAT(p(0).val(), WithinAbs(4.0, epsilon));
    REQUIRE_THAT(p(1).val(), WithinAbs(5.0, epsilon));
  }
}

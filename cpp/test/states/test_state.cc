#include <lupnt/states/state.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("states.state") {
  SECTION("generic State stores metadata and vector values") {
    State x(Vec3(1.0, 2.0, 3.0));
    x.SetName("position");
    x.SetType("CustomState");
    x.SetFrame(Frame::GCRF);
    x.SetNames({"x", "y", "z"});
    x.SetUnits({"m", "m", "m"});

    REQUIRE(x.size() == 3);
    REQUIRE(x.GetName() == "position");
    REQUIRE(x.GetType() == "CustomState");
    REQUIRE(x.GetFrame() == Frame::GCRF);
    REQUIRE(x.GetNames() == std::vector<std::string>{"x", "y", "z"});
    REQUIRE(x.GetUnits() == std::vector<std::string>{"m", "m", "m"});
    REQUIRE_THAT(x(2).val(), WithinAbs(3.0, epsilon));
  }

  SECTION("Cart6 exposes position, velocity, names, units, and frame") {
    Cart6 rv(Vec3(1.0, 2.0, 3.0), Vec3(4.0, 5.0, 6.0), Frame::GCRF);

    REQUIRE(rv.GetType() == Cart6::TYPE);
    REQUIRE(rv.GetFrame() == Frame::GCRF);
    REQUIRE(rv.GetNames() == std::vector<std::string>{"r_x", "r_y", "r_z", "v_x", "v_y", "v_z"});
    REQUIRE(rv.GetUnits() == std::vector<std::string>{"m", "m", "m", "m/s", "m/s", "m/s"});
    REQUIRE_THAT(rv.r()(1).val(), WithinAbs(2.0, epsilon));
    REQUIRE_THAT(rv.v()(2).val(), WithinAbs(6.0, epsilon));
  }

  SECTION("specialized states expose typed accessors") {
    ClassicalOE coe(Vec6(7000.0, 0.1, 0.2, 0.3, 0.4, 0.5), Frame::MOON_CI);
    ClockState3 clock;
    ImuState imu(Vec6(1.0, 2.0, 3.0, 4.0, 5.0, 6.0));

    clock.b() = 1.0e-6;
    clock.d() = 2.0e-9;
    clock.dr() = 3.0e-12;

    REQUIRE(coe.GetType() == ClassicalOE::TYPE);
    REQUIRE_THAT(coe.a().val(), WithinAbs(7000.0, epsilon));
    REQUIRE_THAT(clock.b().val(), WithinAbs(1.0e-6, 1.0e-12));
    REQUIRE_THAT(imu.b_w()(0).val(), WithinAbs(1.0, epsilon));
    REQUIRE_THAT(imu.b_a()(2).val(), WithinAbs(6.0, epsilon));
  }

  SECTION("JointOrbitClockState combines Cartesian orbit and selected clock units") {
    Cart6 orbit(Vec3(1.0, 2.0, 3.0), Vec3(4.0, 5.0, 6.0), Frame::GCRF);
    ClockState3 clock;
    clock.b() = 7.0;
    clock.d() = 8.0;
    clock.dr() = 9.0;
    clock.SetUnits({"m", "m/s", "m/s^2"});

    JointOrbitClockState joint(orbit, clock);
    State split_clock = joint.GetClockState();

    REQUIRE(joint.GetType() == JointOrbitClockState::TYPE);
    REQUIRE(joint.GetFrame() == Frame::GCRF);
    REQUIRE(joint.size() == 9);
    REQUIRE(joint.GetClockStateSize() == 3);
    REQUIRE(joint.GetUnits()
            == std::vector<std::string>{"m", "m", "m", "m/s", "m/s", "m/s", "m", "m/s", "m/s^2"});
    REQUIRE_THAT(joint.r()(2).val(), WithinAbs(3.0, epsilon));
    REQUIRE_THAT(joint.v()(1).val(), WithinAbs(5.0, epsilon));
    REQUIRE_THAT(split_clock(0).val(), WithinAbs(7.0, epsilon));
    REQUIRE(split_clock.GetUnits() == std::vector<std::string>{"m", "m/s", "m/s^2"});
  }
}

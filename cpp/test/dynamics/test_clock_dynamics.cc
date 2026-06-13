#include <lupnt/core/constants.h>
#include <lupnt/devices/clock.h>
#include <lupnt/dynamics/clock_dynamics.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("dynamics.clock_dynamics") {
  SECTION("transition matrices match constant-drift kinematics") {
    Mat2 phi2 = ClockDynamics::TwoStatePhi(10.0);
    Mat3 phi3 = ClockDynamics::ThreeStatePhi(10.0);

    REQUIRE_THAT(phi2(0, 1).val(), WithinAbs(10.0, epsilon));
    REQUIRE_THAT(phi3(0, 1).val(), WithinAbs(10.0, epsilon));
    REQUIRE_THAT(phi3(0, 2).val(), WithinAbs(50.0, epsilon));
    REQUIRE_THAT(phi3(1, 2).val(), WithinAbs(10.0, epsilon));
  }

  SECTION("noise-free propagation applies the clock STM") {
    ClockDynamics dyn;
    dyn.SetAddNoise(false);
    ClockState3 x0;
    x0.b() = 1.0;
    x0.d() = 0.1;
    x0.dr() = 0.01;

    State xf = dyn.Propagate(x0, 0.0, 10.0);

    REQUIRE_THAT(xf(0).val(), WithinAbs(2.5, epsilon));
    REQUIRE_THAT(xf(1).val(), WithinAbs(0.2, epsilon));
    REQUIRE_THAT(xf(2).val(), WithinAbs(0.01, epsilon));
  }

  SECTION("clock process noise can be expressed in range-equivalent units") {
    Mat3 q_seconds = ClockDynamics::ThreeStateNoise(ClockModel::USO, 10.0);
    Mat3 q_meters = ClockDynamics::ThreeStateNoise(ClockModel::USO, 10.0, ClockBiasUnit::METERS);
    Mat3 q_kilometers
        = ClockDynamics::ThreeStateNoise(ClockModel::USO, 10.0, ClockBiasUnit::KILOMETERS);

    REQUIRE_THAT(q_meters(0, 0).val(), WithinRel((q_seconds(0, 0) * C * C).val(), 1.0e-12));
    REQUIRE_THAT(q_kilometers(0, 0).val(),
                 WithinRel((q_seconds(0, 0) * C * KM_M * C * KM_M).val(), 1.0e-12));
  }

  SECTION("instance propagation preserves selected clock state units") {
    ClockDynamics dyn;
    dyn.SetAddNoise(false);
    dyn.SetClockBiasUnit(ClockBiasUnit::KILOMETERS);
    ClockState3 x0;

    State xf = dyn.Propagate(x0, 0.0, 10.0);

    REQUIRE(dyn.GetClockBiasUnit() == ClockBiasUnit::KILOMETERS);
    REQUIRE(xf.GetUnits() == std::vector<std::string>{"km", "km/s", "km/s^2"});
  }

  SECTION("relativistic correction adds deterministic clock bias") {
    ClockDynamics dyn;
    dyn.SetAddNoise(false);
    ClockState3 x0;

    Vec6 rv_vec;
    rv_vec << R_EARTH, 0.0, 0.0, 0.0, 0.0, 0.0;
    ClockRelativityContext relativity;
    relativity.rv0_centered = Cart6(rv_vec, Frame::GCRF);
    relativity.rvf_centered = Cart6(rv_vec, Frame::GCRF);
    relativity.GM = GM_EARTH;
    relativity.c = C;

    State xf = dyn.PropagateWithRelativity(x0, 0.0, 10.0, relativity);
    Real expected_rate = ClockDynamics::RelativisticRateCorrection(
        relativity.rv0_centered.head<3>(), relativity.rv0_centered.tail<3>(), GM_EARTH, C);

    REQUIRE_THAT(xf(0).val(), WithinAbs((expected_rate * 10.0).val(), 1.0e-15));
    REQUIRE_THAT(xf(1).val(), WithinAbs(0.0, epsilon));
    REQUIRE_THAT(xf(2).val(), WithinAbs(0.0, epsilon));
  }

  SECTION("relativistic correction respects range-equivalent clock units") {
    ClockDynamics dyn;
    dyn.SetAddNoise(false);
    dyn.SetClockBiasUnit(ClockBiasUnit::METERS);
    ClockState3 x0;

    Vec6 rv_vec;
    rv_vec << R_EARTH, 0.0, 0.0, 0.0, 0.0, 0.0;
    ClockRelativityContext relativity;
    relativity.rv0_centered = Cart6(rv_vec, Frame::GCRF);
    relativity.rvf_centered = Cart6(rv_vec, Frame::GCRF);
    relativity.GM = GM_EARTH;
    relativity.c = C;

    State xf = dyn.PropagateWithRelativity(x0, 0.0, 10.0, relativity);
    Real expected_rate = ClockDynamics::RelativisticRateCorrection(
        relativity.rv0_centered.head<3>(), relativity.rv0_centered.tail<3>(), GM_EARTH, C);

    REQUIRE_THAT(xf(0).val(), WithinAbs((C * expected_rate * 10.0).val(), 1.0e-6));
    REQUIRE(xf.GetUnits() == std::vector<std::string>{"m", "m/s", "m/s^2"});
  }
}

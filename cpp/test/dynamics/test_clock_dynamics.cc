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

  SECTION("GetClockValues returns the tabulated oscillator PSDs") {
    auto [q1, q2, q3] = ClockDynamics::GetClockValues(ClockModel::USO);
    REQUIRE_THAT(q1, WithinRel(1.2e-22, 1e-12));
    REQUIRE_THAT(q2, WithinRel(1.58e-26, 1e-12));
    REQUIRE_THAT(q3, WithinRel(5.9e-75, 1e-12));
  }

  SECTION("bias-unit scaling and round-trip conversions") {
    REQUIRE_THAT(ClockDynamics::SecondsToBiasUnitScale(ClockBiasUnit::SECONDS),
                 WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(ClockDynamics::SecondsToBiasUnitScale(ClockBiasUnit::METERS),
                 WithinRel(C, 1e-12));

    // SecondsToBiasUnits and BiasUnitsToSeconds are mutual inverses.
    const Real value_s = 3.5e-9;  // 3.5 ns
    for (auto unit : {ClockBiasUnit::SECONDS, ClockBiasUnit::METERS, ClockBiasUnit::KILOMETERS}) {
      Real in_unit = ClockDynamics::SecondsToBiasUnits(value_s, unit);
      REQUIRE_THAT(ClockDynamics::BiasUnitsToSeconds(in_unit, unit).val(),
                   WithinRel(value_s.val(), 1e-12));
    }
    // Meters bias equals seconds bias times the speed of light.
    REQUIRE_THAT(ClockDynamics::SecondsToBiasUnits(value_s, ClockBiasUnit::METERS).val(),
                 WithinRel((value_s * C).val(), 1e-12));
  }

  SECTION("TwoStateNoise matches the CWNA closed form and scales with units") {
    const Real dt = 10.0;
    auto [q1, q2, q3] = ClockDynamics::GetClockValues(ClockModel::USO);
    Mat2 Q = ClockDynamics::TwoStateNoise(ClockModel::USO, dt);

    REQUIRE_THAT(Q(0, 0).val(), WithinRel((q1 * dt + q2 * dt * dt * dt / 3.0).val(), 1e-10));
    REQUIRE_THAT(Q(0, 1).val(), WithinRel((q2 * dt * dt / 2.0).val(), 1e-10));
    REQUIRE_THAT(Q(1, 0).val(), WithinRel(Q(0, 1).val(), 1e-12));
    REQUIRE_THAT(Q(1, 1).val(), WithinRel((q2 * dt).val(), 1e-10));

    // Range-equivalent (meters) noise is the seconds noise scaled by c^2.
    Mat2 Q_m = ClockDynamics::TwoStateNoise(ClockModel::USO, dt, ClockBiasUnit::METERS);
    REQUIRE_THAT(Q_m(0, 0).val(), WithinRel((Q(0, 0) * C * C).val(), 1e-10));
  }
}

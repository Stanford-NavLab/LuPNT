#include <lupnt/core/constants.h>
#include <lupnt/dynamics/joint_orbit_clock_dynamics.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  class StaticOrbitDynamics : public NumericalDynamics {
  public:
    StaticOrbitDynamics() : NumericalDynamics() {
      SetODE([this](Real t, const VecX& x) { return ComputeRates(t, State(x)); });
    }

    VecX ComputeRates(Real t, const State& x) const override {
      (void)t;
      return VecX::Zero(x.size());
    }
  };
}  // namespace

TEST_CASE("dynamics.joint_orbit_clock_dynamics") {
  constexpr double eps = 1.0e-12;

  SECTION("propagates relativistic clock bias with a static centered orbit") {
    Cart6 orbit(Vec3(R_EARTH, 0.0, 0.0), Vec3::Zero(), Frame::GCRF);
    ClockState3 clock;
    JointOrbitClockState x0(orbit, clock);

    auto orbit_dynamics = MakePtr<StaticOrbitDynamics>();
    auto clock_dynamics = MakePtr<ClockDynamics>();
    clock_dynamics->SetAddNoise(false);

    JointOrbitClockDynamics dynamics;
    dynamics.SetOrbitDynamics(orbit_dynamics);
    dynamics.SetClockDynamics(clock_dynamics);
    dynamics.SetFrame(Frame::GCRF);
    dynamics.SetRelativityCenterBody(BodyId::EARTH);

    JointOrbitClockState xf(dynamics.Propagate(x0, 0.0, 10.0));
    Real expected_rate
        = ClockDynamics::RelativisticRateCorrection(orbit.r(), orbit.v(), GM_EARTH, C);

    REQUIRE_THAT(xf.r()(0).val(), WithinAbs(R_EARTH, eps));
    REQUIRE_THAT(xf.b().val(), WithinAbs((expected_rate * 10.0).val(), 1.0e-15));
    REQUIRE_THAT(xf.d().val(), WithinAbs(0.0, eps));
    REQUIRE_THAT(xf.dr().val(), WithinAbs(0.0, eps));
  }

  SECTION("converts relativistic clock bias into range-equivalent meters") {
    Cart6 orbit(Vec3(R_EARTH, 0.0, 0.0), Vec3::Zero(), Frame::GCRF);
    ClockState2 clock;
    clock.SetUnits({"m", "m/s"});
    JointOrbitClockState x0(orbit, clock);

    auto orbit_dynamics = MakePtr<StaticOrbitDynamics>();
    auto clock_dynamics = MakePtr<ClockDynamics>();
    clock_dynamics->SetAddNoise(false);
    clock_dynamics->SetClockBiasUnit(ClockBiasUnit::METERS);

    JointOrbitClockDynamics dynamics;
    dynamics.SetOrbitDynamics(orbit_dynamics);
    dynamics.SetClockDynamics(clock_dynamics);
    dynamics.SetFrame(Frame::GCRF);
    dynamics.SetRelativityCenterBody(BodyId::EARTH);

    JointOrbitClockState xf(dynamics.Propagate(x0, 0.0, 10.0));
    Real expected_rate
        = ClockDynamics::RelativisticRateCorrection(orbit.r(), orbit.v(), GM_EARTH, C);

    REQUIRE(xf.GetClockStateSize() == 2);
    REQUIRE_THAT(xf.b().val(), WithinAbs((C * expected_rate * 10.0).val(), 1.0e-6));
    REQUIRE_THAT(xf.d().val(), WithinAbs(0.0, eps));
    REQUIRE(xf.GetClockState().GetUnits() == std::vector<std::string>{"m", "m/s"});
  }
}

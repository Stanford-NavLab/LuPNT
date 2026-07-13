#include <lupnt/core/constants.h>
#include <lupnt/dynamics/joint_orbit_clock_dynamics.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  // Straight-line kinematic orbit: rdot = v, vdot = 0.
  class KinematicOrbitDynamics : public NumericalDynamics {
  public:
    KinematicOrbitDynamics() : NumericalDynamics() {
      SetODE([this](Real t, const VecX& x) { return ComputeRates(t, State(x)); });
    }
    VecX ComputeRates(Real t, const State& x) const override {
      (void)t;
      VecX rates = VecX::Zero(x.size());
      rates.head(3) = x.segment(3, 3);
      return rates;
    }
  };
}  // namespace

TEST_CASE("dynamics.joint_orbit_clock_dynamics_more.units_accessor") {
  JointOrbitClockDynamics dyn;
  // Default is SI; setting a different unit system is echoed back.
  dyn.SetUnits(SI_UNITS);
  REQUIRE(dyn.GetUnits().length == SI_UNITS.length);
  REQUIRE(dyn.GetUnits().time == SI_UNITS.time);
}

TEST_CASE("dynamics.joint_orbit_clock_dynamics_more.zero_interval_is_identity") {
  // |tf - t0| < EPS skips integration; the state must come back unchanged
  // (aside from clock noise, which is disabled here).
  Cart6 orbit_state(Vec3(7000.0e3, 0.0, 0.0), Vec3(0.0, 7500.0, 0.0), Frame::GCRF);
  ClockState3 clock0;
  clock0.b() = 1.0e-6;
  clock0.d() = 2.0e-9;
  clock0.dr() = 3.0e-13;
  JointOrbitClockState x0(orbit_state, clock0);

  auto orbit = MakePtr<KinematicOrbitDynamics>();
  auto clock = MakePtr<ClockDynamics>();
  clock->SetAddNoise(false);
  JointOrbitClockDynamics dyn;
  dyn.SetOrbitDynamics(orbit);
  dyn.SetClockDynamics(clock);
  dyn.SetFrame(Frame::GCRF);
  dyn.SetUseClockRelativity(false);

  JointOrbitClockState xf(dyn.Propagate(x0, 100.0, 100.0));  // tf == t0
  for (int i = 0; i < x0.size(); ++i) REQUIRE_THAT(xf(i).val(), WithinAbs(x0(i).val(), 1e-12));
}

TEST_CASE("dynamics.joint_orbit_clock_dynamics_more.relativity_rate_term") {
  // With relativity enabled and a fixed Earth center in the GCRF frame (whose
  // natural center is Earth), the relativistic clock-rate correction is
  // evaluated locally (no ephemeris lookup) and adds a small term onto the
  // drift in the clock-bias derivative (rates(6)).
  Cart6 orbit_state(Vec3(7000.0e3, 0.0, 0.0), Vec3(0.0, 7546.0, 0.0), Frame::GCRF);
  ClockState3 clock0;
  clock0.b() = 0.0;
  clock0.d() = 5.0e-9;
  clock0.dr() = 0.0;
  JointOrbitClockState x0(orbit_state, clock0);

  auto orbit = MakePtr<KinematicOrbitDynamics>();
  auto clock = MakePtr<ClockDynamics>();

  JointOrbitClockDynamics dyn;
  dyn.SetOrbitDynamics(orbit);
  dyn.SetClockDynamics(clock);
  dyn.SetFrame(Frame::GCRF);
  dyn.SetRelativityCenterBody(BodyId::EARTH);

  dyn.SetUseClockRelativity(false);
  VecX rates_off = dyn.ComputeRates(0.0, x0);

  dyn.SetUseClockRelativity(true);
  VecX rates_on = dyn.ComputeRates(0.0, x0);

  // Relativity off: rates(6) == drift exactly.
  REQUIRE_THAT(rates_off(6).val(), WithinAbs(5.0e-9, 1e-18));

  // Relativity on: a small, finite correction is added.
  REQUIRE(std::isfinite(rates_on(6).val()));
  double correction = rates_on(6).val() - rates_off(6).val();
  REQUIRE(correction != 0.0);
  REQUIRE(std::abs(correction) < 1e-8);  // relativistic scale ~ (v^2/2 + GM/r)/c^2

  // The orbit and drift-rate blocks are unaffected by the relativity toggle.
  for (int i = 0; i < 6; ++i) REQUIRE_THAT(rates_on(i).val(), WithinAbs(rates_off(i).val(), 1e-15));
  REQUIRE_THAT(rates_on(7).val(), WithinAbs(rates_off(7).val(), 1e-24));
}

TEST_CASE("dynamics.joint_orbit_clock_dynamics_more.two_state_relativity_shape") {
  // A 2-state clock [b, d] with relativity on: the derivative is
  // [drift + relativistic_rate, 0] (no drift-rate slot).
  Cart6 orbit_state(Vec3(6800.0e3, 100.0e3, 0.0), Vec3(0.0, 7600.0, 10.0), Frame::GCRF);
  ClockState2 clock0;
  clock0.b() = 1.0e-6;
  clock0.d() = 4.0e-9;
  JointOrbitClockState x0(orbit_state, clock0);

  auto orbit = MakePtr<KinematicOrbitDynamics>();
  auto clock = MakePtr<ClockDynamics>();
  JointOrbitClockDynamics dyn;
  dyn.SetOrbitDynamics(orbit);
  dyn.SetClockDynamics(clock);
  dyn.SetFrame(Frame::GCRF);
  dyn.SetRelativityCenterBody(BodyId::EARTH);
  dyn.SetUseClockRelativity(true);

  VecX rates = dyn.ComputeRates(0.0, x0);
  REQUIRE(rates.size() == x0.size());  // 8 elements for a 2-state clock
  REQUIRE(x0.GetClockStateSize() == 2);
  // Clock-bias rate is drift plus a small relativistic correction.
  REQUIRE(std::isfinite(rates(6).val()));
  REQUIRE(std::abs(rates(6).val() - 4.0e-9) < 1e-8);
  // Last clock slot (drift term of the 2-state model) is zero.
  REQUIRE_THAT(rates(7).val(), WithinAbs(0.0, 1e-18));
}

#include <lupnt/core/constants.h>
#include <lupnt/dynamics/joint_orbit_clock_dynamics.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  // Kinematic straight-line orbit dynamics: rdot = v, vdot = 0. This gives a
  // closed-form propagation (r = r0 + v0*dt, v = v0) that both the standalone
  // model and the joint model must reproduce for the orbit sub-block.
  class KinematicOrbitDynamics : public NumericalDynamics {
  public:
    KinematicOrbitDynamics() : NumericalDynamics() {
      SetODE([this](Real t, const VecX& x) { return ComputeRates(t, State(x)); });
    }
    VecX ComputeRates(Real t, const State& x) const override {
      (void)t;
      VecX rates = VecX::Zero(x.size());
      rates.head(3) = x.segment(3, 3);  // rdot = v
      return rates;                     // vdot = 0
    }
  };

  JointOrbitClockDynamics MakeDecoupledJoint(const Ptr<NumericalDynamics>& orbit,
                                             const Ptr<ClockDynamics>& clock) {
    JointOrbitClockDynamics dyn;
    dyn.SetOrbitDynamics(orbit);
    dyn.SetClockDynamics(clock);
    dyn.SetFrame(Frame::GCRF);
    dyn.SetUseClockRelativity(false);  // pure bias/drift integrator
    return dyn;
  }
}  // namespace

TEST_CASE("dynamics.joint_orbit_clock_dynamics_extra.accessors") {
  JointOrbitClockDynamics dyn;
  auto orbit = MakePtr<KinematicOrbitDynamics>();
  auto clock = MakePtr<ClockDynamics>();

  dyn.SetOrbitDynamics(orbit);
  dyn.SetClockDynamics(clock);
  REQUIRE(dyn.GetOrbitDynamics() == orbit);
  REQUIRE(dyn.GetClockDynamics() == clock);

  dyn.SetUseClockRelativity(false);
  REQUIRE_FALSE(dyn.GetUseClockRelativity());
  dyn.SetUseClockRelativity(true);
  REQUIRE(dyn.GetUseClockRelativity());

  REQUIRE_FALSE(dyn.HasRelativityCenterBody());
  dyn.SetRelativityCenterBody(BodyId::MOON);
  REQUIRE(dyn.HasRelativityCenterBody());
  REQUIRE(dyn.GetRelativityCenterBody() == BodyId::MOON);
  dyn.ClearRelativityCenterBody();
  REQUIRE_FALSE(dyn.HasRelativityCenterBody());

  dyn.SetReferenceRateOffset(1.5e-10);
  REQUIRE_THAT(dyn.GetReferenceRateOffset().val(), WithinAbs(1.5e-10, 1e-20));

  dyn.SetAddClockNoise(true);
  REQUIRE(dyn.GetAddClockNoise());
  dyn.SetAddClockNoise(false);
  REQUIRE_FALSE(dyn.GetAddClockNoise());

  dyn.SetFrame(Frame::MOON_CI);
  REQUIRE(dyn.GetFrame() == Frame::MOON_CI);

  REQUIRE(dyn.GetStateType() == JointOrbitClockState::TYPE);
}

TEST_CASE("dynamics.joint_orbit_clock_dynamics_extra.decoupled_clock_evolution") {
  // With relativity disabled the clock block is an independent integrator, so
  // the propagated clock must follow the exact polynomial law while the orbit
  // block follows the straight-line kinematic law.
  Cart6 orbit_state(Vec3(7000.0e3, 100.0e3, -50.0e3), Vec3(10.0, -3.0, 2.0), Frame::GCRF);

  auto orbit = MakePtr<KinematicOrbitDynamics>();
  auto clock = MakePtr<ClockDynamics>();
  clock->SetAddNoise(false);

  const Real dt = 120.0;

  SECTION("3-state clock [b, d, dr] follows the double-integrator solution") {
    ClockState3 clock0;
    clock0.b() = 5.0e-6;
    clock0.d() = 2.0e-9;
    clock0.dr() = 1.0e-13;
    JointOrbitClockState x0(orbit_state, clock0);

    JointOrbitClockDynamics dyn = MakeDecoupledJoint(orbit, clock);
    JointOrbitClockState xf(dyn.Propagate(x0, 0.0, dt));

    // Clock: b(t)=b0+d0*t+0.5*dr0*t^2, d(t)=d0+dr0*t, dr(t)=dr0.
    double b0 = 5.0e-6, d0 = 2.0e-9, dr0 = 1.0e-13, t = dt.val();
    REQUIRE_THAT(xf.b().val(), WithinAbs(b0 + d0 * t + 0.5 * dr0 * t * t, 1e-15));
    REQUIRE_THAT(xf.d().val(), WithinAbs(d0 + dr0 * t, 1e-18));
    REQUIRE_THAT(xf.dr().val(), WithinAbs(dr0, 1e-24));

    // Orbit: straight-line kinematics.
    for (int i = 0; i < 3; ++i)
      REQUIRE_THAT(xf.r()(i).val(),
                   WithinAbs((orbit_state.r()(i) + orbit_state.v()(i) * dt).val(), 1e-3));
    for (int i = 0; i < 3; ++i)
      REQUIRE_THAT(xf.v()(i).val(), WithinAbs(orbit_state.v()(i).val(), 1e-9));
  }

  SECTION("2-state clock [b, d] evolves as [b0 + d0*dt, d0]") {
    ClockState2 clock0;
    clock0.b() = -3.0e-6;
    clock0.d() = 4.0e-9;
    JointOrbitClockState x0(orbit_state, clock0);

    JointOrbitClockDynamics dyn = MakeDecoupledJoint(orbit, clock);
    JointOrbitClockState xf(dyn.Propagate(x0, 0.0, dt));

    REQUIRE(xf.GetClockStateSize() == 2);
    REQUIRE_THAT(xf.b().val(), WithinAbs((-3.0e-6 + 4.0e-9 * dt).val(), 1e-15));
    REQUIRE_THAT(xf.d().val(), WithinAbs(4.0e-9, 1e-18));
  }
}

TEST_CASE("dynamics.joint_orbit_clock_dynamics_extra.orbit_block_matches_standalone") {
  // The joint model's orbit sub-block must match propagating the underlying
  // orbit dynamics on its own over the same interval.
  Cart6 orbit_state(Vec3(6800.0e3, -200.0e3, 90.0e3), Vec3(-4.0, 7.5, 1.2), Frame::GCRF);
  ClockState3 clock0;
  clock0.b() = 1.0e-6;
  clock0.d() = 3.0e-9;
  JointOrbitClockState x0(orbit_state, clock0);

  auto orbit = MakePtr<KinematicOrbitDynamics>();
  auto clock = MakePtr<ClockDynamics>();
  clock->SetAddNoise(false);
  JointOrbitClockDynamics dyn = MakeDecoupledJoint(orbit, clock);

  Real dt = 300.0;
  JointOrbitClockState xf(dyn.Propagate(x0, 0.0, dt));

  Cart6 orbit_ref(orbit->Propagate(orbit_state, 0.0, dt));
  for (int i = 0; i < 6; ++i)
    REQUIRE_THAT(xf.head(6)(i).val(), WithinAbs(orbit_ref(i).val(), 1e-6));
}

TEST_CASE("dynamics.joint_orbit_clock_dynamics_extra.compute_rates_structure") {
  // ComputeRates must place the orbit rates in the head(6) block and the clock
  // derivative [drift(+rel), drift_rate, 0] in the tail, with relativity off.
  Cart6 orbit_state(Vec3(8000.0e3, 0.0, 0.0), Vec3(0.0, 7000.0, 0.0), Frame::GCRF);
  ClockState3 clock0;
  clock0.b() = 2.0e-6;
  clock0.d() = 5.0e-9;
  clock0.dr() = 3.0e-13;
  JointOrbitClockState x0(orbit_state, clock0);

  auto orbit = MakePtr<KinematicOrbitDynamics>();
  auto clock = MakePtr<ClockDynamics>();
  JointOrbitClockDynamics dyn = MakeDecoupledJoint(orbit, clock);

  VecX rates = dyn.ComputeRates(0.0, x0);
  REQUIRE(rates.size() == x0.size());

  // Orbit block: rdot = v, vdot = 0.
  for (int i = 0; i < 3; ++i)
    REQUIRE_THAT(rates(i).val(), WithinAbs(orbit_state.v()(i).val(), 1e-12));
  for (int i = 3; i < 6; ++i) REQUIRE_THAT(rates(i).val(), WithinAbs(0.0, 1e-12));

  // Clock block: [drift, drift_rate, 0] (relativity off).
  REQUIRE_THAT(rates(6).val(), WithinAbs(5.0e-9, 1e-18));
  REQUIRE_THAT(rates(7).val(), WithinAbs(3.0e-13, 1e-22));
  REQUIRE_THAT(rates(8).val(), WithinAbs(0.0, 1e-24));
}

TEST_CASE("dynamics.joint_orbit_clock_dynamics_extra.stm_shape") {
  // The inherited Dynamics STM entry point must return a square matrix matching
  // the joint state size, with finite, non-degenerate entries.
  Cart6 orbit_state(Vec3(7000.0e3, 0.0, 0.0), Vec3(0.0, 7500.0, 0.0), Frame::GCRF);
  ClockState3 clock0;
  clock0.d() = 1.0e-9;
  JointOrbitClockState x0(orbit_state, clock0);

  auto orbit = MakePtr<KinematicOrbitDynamics>();
  auto clock = MakePtr<ClockDynamics>();
  clock->SetAddNoise(false);
  JointOrbitClockDynamics dyn = MakeDecoupledJoint(orbit, clock);

  MatXd stm;
  dyn.Propagate(x0, 0.0, 60.0, nullptr, &stm);

  REQUIRE(stm.rows() == x0.size());
  REQUIRE(stm.cols() == x0.size());
  for (int i = 0; i < stm.rows(); ++i)
    for (int j = 0; j < stm.cols(); ++j) REQUIRE(std::isfinite(stm(i, j)));
  // Diagonal must be non-zero (identity-dominated for a short arc).
  for (int i = 0; i < stm.rows(); ++i) REQUIRE(std::abs(stm(i, i)) > 0.0);
}

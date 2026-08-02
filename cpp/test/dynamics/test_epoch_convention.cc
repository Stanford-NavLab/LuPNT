// Pins the dynamics-layer time convention.
//
// `Dynamics::Propagate(x0, t0, tf, ...)` takes RELATIVE simulation time --
// seconds from the global reference epoch `GetLupntEpoch()` -- not an absolute
// epoch. Consumers convert internally (`NBodyDynamics` and
// `FixedPointingDynamics` both compute `t + GetLupntEpoch()`).
//
// Nothing enforced this: both are plain `Real`, so passing an absolute epoch,
// or dropping a consumer's internal offset, compiles fine and silently
// evaluates the ephemeris decades away. Removing the offset from
// FixedPointingDynamics passed the entire 421-case suite, which is why this
// test exists.

#include <lupnt/conversions/time_conversions.h>
#include <lupnt/core/definitions.h>
#include <lupnt/dynamics/attitude_dynamics.h>
#include <lupnt/environment/body.h>
#include <lupnt/states/state.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

TEST_CASE("dynamics.epoch_convention_is_relative") {
  Vec6 rv_vec;
  rv_vec << 7.0e6, 0.0, 0.0, 0.0, 7.5e3, 0.0;
  Cart6 rv(rv_vec, Frame::GCRF);
  Attitude qw0;
  qw0.SetFrame(Frame::GCRF);

  FixedPointingDynamics dyn;
  dyn.SetPrimaryPointing(Axis::X, BodyId::EARTH);
  dyn.SetSecondaryPointing(Axis::Y, BodyId::SUN);

  const Real saved = GetLupntEpoch();

  // Two setups that denote the SAME absolute instant:
  //   (A) reference epoch E, relative time  d
  //   (B) reference epoch E + d, relative time 0
  // Because `t` is relative, both must give the same attitude.
  const Real E = GregorianToTime("2030-01-01T00:00:00");
  const Real d = 3600.0;

  SetLupntEpoch(E);
  Attitude a = Attitude(dyn.Propagate(qw0, Real(0.0), d, rv));

  SetLupntEpoch(E + d);
  Attitude b = Attitude(dyn.Propagate(qw0, Real(0.0), Real(0.0), rv));

  SetLupntEpoch(saved);

  double max_diff = 0.0;
  for (int i = 0; i < 4; ++i) max_diff = std::max(max_diff, std::abs((a(i) - b(i)).val()));
  INFO("max |quaternion difference| = " << max_diff);
  // Same absolute instant reached two ways -> identical attitude. This fails if
  // a consumer's `+ GetLupntEpoch()` offset is dropped or double-applied.
  REQUIRE_THAT(max_diff, WithinAbs(0.0, 1.0e-12));
}

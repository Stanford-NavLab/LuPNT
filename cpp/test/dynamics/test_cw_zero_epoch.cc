// Regression test for the ClohessyWiltshireDynamics::Propagate crash when the
// reference epoch t0 == 0. The integration-constant cache key t0_
// default-initializes to 0, so the old `if (t0 != t0_)` guard was false on a
// first call with t0 == 0, leaving K_ empty and crashing at `Phi * K_`. The fix
// also solves K_ whenever it is empty and updates t0_ so the cache tracks the
// last-solved epoch.

#include <lupnt/lupnt.h>

#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "../utils.cc"

using namespace lupnt;

TEST_CASE("dynamics.clohessy_wiltshire_zero_epoch") {
  const Real a = 1.0;
  const Real n = 1.1e-3;  // ~LEO mean motion [rad/s]
  ClohessyWiltshireDynamics cw(a, n);
  RelCart6 x0(Vec3(10.0, -5.0, 3.0), Vec3(0.1, -0.2, 0.05), Frame::GCRF);

  SECTION("Propagate starting at t0 == 0 returns a finite state (was a crash)") {
    MatXd stm;
    State xf = cw.Propagate(x0, Real(0.0), Real(600.0), nullptr, &stm);
    for (int i = 0; i < 6; ++i) REQUIRE(std::isfinite(xf(i).val()));
    REQUIRE(stm.rows() == 6);
    REQUIRE(stm.cols() == 6);
    REQUIRE(stm.allFinite());
  }

  SECTION("A subsequent call at a different t0 also stays finite (cache updates)") {
    State x_a = cw.Propagate(x0, Real(0.0), Real(300.0), nullptr, nullptr);
    State x_b = cw.Propagate(x0, Real(120.0), Real(420.0), nullptr, nullptr);
    for (int i = 0; i < 6; ++i) {
      REQUIRE(std::isfinite(x_a(i).val()));
      REQUIRE(std::isfinite(x_b(i).val()));
    }
  }
}

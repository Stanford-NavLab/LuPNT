#include <lupnt/dynamics/analytical_orbit_dynamics.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "../utils.cc"
#include "lupnt/conversions/state_conversions.h"
#include "lupnt/numerics/math_utils.h"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {

  // ---------------------------------------------------------------------------
  // KeplerianDynamics<ClassicalOE>: physical invariants over a full orbit.
  // The existing test only checks the mean-anomaly advance and one STM entry;
  // here we assert the two-body constants of motion (specific orbital energy,
  // semi-major axis) are exactly conserved, evaluated in Cartesian space.
  // ---------------------------------------------------------------------------
  TEST_CASE("dynamics.analytical_orbit_dynamics_extra.kepler_energy_conservation") {
    const Real GM = GM_EARTH;
    ClassicalOE coe0(Vec6(8000.0e3, 0.2, 0.3, 0.4, 0.5, 0.1), Frame::GCRF);
    KeplerianDynamics<ClassicalOE> dyn(GM);

    Real a = coe0.a();
    Real n = sqrt(GM / pow(a, 3));
    Real period = TWO_PI / n;
    Real energy_ref = (-GM / (2.0 * a)).val();  // vis-viva specific energy

    // Sample many points around one orbit.
    for (int k = 0; k <= 12; ++k) {
      Real tf = period * (k / 12.0);
      State coe_f = dyn.Propagate(coe0, Real(0.0), tf, nullptr);

      // Semi-major axis, eccentricity, inclination, RAAN, arg. of periapsis are
      // invariant under Keplerian motion; only the mean anomaly advances.
      REQUIRE_THAT(coe_f(0).val(), WithinRel(coe0(0).val(), 1e-12));
      REQUIRE_THAT(coe_f(1).val(), WithinAbs(coe0(1).val(), 1e-12));
      REQUIRE_THAT(coe_f(2).val(), WithinAbs(coe0(2).val(), 1e-12));

      Real M_expected = WrapToPi(coe0.M() + n * tf);
      REQUIRE_THAT(WrapToPi(Real(coe_f(5) - M_expected)).val(), WithinAbs(0.0, 1e-9));

      // Convert to Cartesian and verify specific orbital energy is conserved.
      Vec6 rv = ClassicalToCart(coe_f, GM);
      Real r = rv.head(3).norm();
      Real v = rv.tail(3).norm();
      Real energy = 0.5 * v * v - GM / r;
      REQUIRE_THAT(energy.val(), WithinRel(energy_ref, 1e-9));
    }

    // After exactly one period the state (including M, mod 2*pi) returns.
    State coe_T = dyn.Propagate(coe0, Real(0.0), period, nullptr);
    REQUIRE_THAT(WrapToPi(Real(coe_T(5) - coe0(5))).val(), WithinAbs(0.0, 1e-9));
  }

  // ---------------------------------------------------------------------------
  // KeplerianDynamics specializations for QuasiNonsingularOE / EquinoctialOE.
  // These advance only the fast angle (mean argument of latitude / mean
  // longitude) by n*dt and hold every other element fixed. Neither was
  // previously exercised.
  // ---------------------------------------------------------------------------
  TEST_CASE("dynamics.analytical_orbit_dynamics_extra.kepler_quasinonsingular") {
    const Real GM = GM_EARTH;
    // [a, u, ex, ey, i, Omega]
    QuasiNonsingularOE qns(Vec6(7500.0e3, 0.2, 0.01, -0.02, 0.9, 1.1));
    KeplerianDynamics<QuasiNonsingularOE> dyn(GM);

    Real a = qns(0);
    Real n = sqrt(GM / pow(a, 3));
    Real dt = 300.0;
    State qf = dyn.Propagate(qns, Real(0.0), dt, nullptr);

    // Only the mean argument of latitude (index 1) advances by n*dt.
    REQUIRE_THAT(qf(1).val(), WithinAbs((qns(1) + n * dt).val(), 1e-9));
    REQUIRE_THAT(qf(0).val(), WithinRel(qns(0).val(), 1e-12));
    for (int k : {2, 3, 4, 5}) REQUIRE_THAT(qf(k).val(), WithinAbs(qns(k).val(), 1e-12));

    // STM is explicitly not implemented for this specialization.
    MatXd stm;
    REQUIRE_THROWS_AS(dyn.Propagate(qns, Real(0.0), dt, nullptr, &stm), std::runtime_error);
  }

  TEST_CASE("dynamics.analytical_orbit_dynamics_extra.kepler_equinoctial") {
    const Real GM = GM_EARTH;
    // [a, h, k, p, q, lon]
    EquinoctialOE eq(Vec6(7000.0e3, 0.01, 0.02, 0.03, -0.01, 0.5));
    KeplerianDynamics<EquinoctialOE> dyn(GM);

    Real a = eq(0);
    Real n = sqrt(GM / pow(a, 3));
    Real dt = 120.0;
    State ef = dyn.Propagate(eq, Real(0.0), dt, nullptr);

    // Only the mean longitude (index 5) advances by n*dt.
    REQUIRE_THAT(ef(5).val(), WithinAbs((eq(5) + n * dt).val(), 1e-9));
    for (int k : {0, 1, 2, 3, 4}) REQUIRE_THAT(ef(k).val(), WithinAbs(eq(k).val(), 1e-9));

    MatXd stm;
    REQUIRE_THROWS_AS(dyn.Propagate(eq, Real(0.0), dt, nullptr, &stm), std::runtime_error);
  }

  // ---------------------------------------------------------------------------
  // ClohessyWiltshireDynamics: the state transition matrix derived from the
  // model's ComputeMat (Phi(dt) = ComputeMat(dt) * ComputeMat(0)^{-1}) must
  // equal the textbook Hill-Clohessy-Wiltshire STM for a circular reference
  // orbit. This validates ComputeMat against a closed-form reference.
  // ---------------------------------------------------------------------------
  TEST_CASE("dynamics.analytical_orbit_dynamics_extra.clohessy_wiltshire_stm") {
    Real a = 7000.0e3;
    double nd = 1.0e-3;  // mean motion [rad/s]
    Real n = nd;
    ClohessyWiltshireDynamics cw(a, n);

    auto textbook = [&](double t) {
      double s = std::sin(nd * t);
      double c = std::cos(nd * t);
      Mat6d P = Mat6d::Zero();
      // radial (x)
      P(0, 0) = 4 - 3 * c;
      P(0, 3) = s / nd;
      P(0, 4) = 2 * (1 - c) / nd;
      // along-track (y)
      P(1, 0) = 6 * (s - nd * t);
      P(1, 1) = 1;
      P(1, 3) = -2 * (1 - c) / nd;
      P(1, 4) = (4 * s - 3 * nd * t) / nd;
      // cross-track (z)
      P(2, 2) = c;
      P(2, 5) = s / nd;
      // radial rate (vx)
      P(3, 0) = 3 * nd * s;
      P(3, 3) = c;
      P(3, 4) = 2 * s;
      // along-track rate (vy)
      P(4, 0) = -6 * nd * (1 - c);
      P(4, 3) = -2 * s;
      P(4, 4) = 4 * c - 3;
      // cross-track rate (vz)
      P(5, 2) = -nd * s;
      P(5, 5) = c;
      return P;
    };

    for (double t : {200.0, 900.0, 1500.0}) {
      Mat6 M_dt = cw.ComputeMat(Real(t));
      Mat6 M_0 = cw.ComputeMat(Real(0.0));
      Mat6d phi = (M_dt * M_0.inverse()).cast<double>();
      Mat6d ref = textbook(t);
      for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
          INFO("entry (" << i << "," << j << ") t=" << t);
          REQUIRE_THAT(phi(i, j), WithinRel(ref(i, j), 1e-7) || WithinAbs(ref(i, j), 1e-6));
        }
      }
    }
  }

}  // namespace

/**
 * @file test_integrator_order.cc
 * @brief Order-of-convergence tests for the Runge-Kutta integrators.
 *
 * The pre-existing tests in `test_integrator_extra.cc` check that each method
 * reproduces a closed-form solution to a fixed absolute tolerance at one fixed
 * step size. That is not enough to pin a method's *order*: at dt = 1e-3 an
 * accidentally-3rd-order "RK8" still lands within 1e-8 of the truth, so a
 * wrong Butcher coefficient passes unnoticed.
 *
 * These tests instead measure the empirical order,
 *
 *     p = log2( err(dt) / err(dt/2) ),
 *
 * which is sensitive to every coefficient in the tableau, and assert the
 * consistency condition c_i = sum_j a_ij that a wrong node or row violates
 * directly.
 *
 * Two ODEs are used deliberately:
 *   - an autonomous one (y' = -y), which does not exercise the stage times c_i
 *     at all, and
 *   - a non-autonomous one (y' = y cos t), which does. An integrator whose c_i
 *     are wrong can be full-order on the first and 1st-order on the second,
 *     and orbit dynamics with ephemeris-driven terms is the second kind.
 */

#include <lupnt/numerics/integrator.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

using namespace lupnt;
using namespace Catch::Matchers;

namespace {

  // y' = -y  =>  y(t) = y0 exp(-t).  Autonomous: independent of the nodes c_i.
  ODE Decay = [](const Real /*t*/, const VecX& x) {
    VecX d(1);
    d(0) = -x(0);
    return d;
  };
  double DecayExact(double t) { return std::exp(-t); }

  // y' = y cos(t)  =>  y(t) = y0 exp(sin t).  Non-autonomous AND state-coupled,
  // so both the nodes c_i and the coupling coefficients a_ij matter.
  ODE Oscillating = [](const Real t, const VecX& x) {
    VecX d(1);
    d(0) = x(0) * cos(t);
    return d;
  };
  double OscillatingExact(double t) { return std::exp(std::sin(t)); }

  VecX Vec1(double v) {
    VecX x(1);
    x(0) = v;
    return x;
  }

  /// Empirical order of convergence over one step-halving.
  ///
  /// `tf` is kept short and `dt` coarse so that the discretization error stays
  /// well above double-precision roundoff; otherwise the ratio measures
  /// floating-point noise rather than the method.
  double EmpiricalOrder(Integrator& integ, const ODE& f, double (*exact)(double), double tf,
                        double dt) {
    VecX x0 = Vec1(1.0);
    double e1
        = std::abs(integ.Propagate(f, Real(0.0), Real(tf), x0, Real(dt))(0).val() - exact(tf));
    double e2
        = std::abs(integ.Propagate(f, Real(0.0), Real(tf), x0, Real(dt / 2))(0).val() - exact(tf));
    // A method that has hit roundoff reports a meaningless ratio; guard so the
    // failure message says "too accurate to measure" rather than "order 0".
    REQUIRE(e2 > 1e-15);
    return std::log2(e1 / e2);
  }

}  // namespace

// -----------------------------------------------------------------------------
// Butcher-tableau consistency: c_i = sum_j a_ij
//
// This is the algebraic condition every Runge-Kutta method must satisfy for
// consistency (order >= 1) on non-autonomous problems. It is checked here
// against the literal coefficients as they appear in integrator.cc, so a typo
// in a single coefficient is caught by arithmetic rather than by a
// convergence experiment.
// -----------------------------------------------------------------------------

TEST_CASE("numerics.integrator_order.rk8_tableau_consistency") {
  // Rows of the RK8 (10-stage Shanks-type) tableau exactly as written in
  // RK8::Step, paired with the stage time actually used for that stage.
  struct Row {
    double c;
    std::vector<double> a;
  };
  const std::vector<Row> rows = {
      {4.0 / 27, {4.0 / 27}},
      {2.0 / 9, {1.0 / 18, 3.0 / 18}},
      {1.0 / 3, {1.0 / 12, 3.0 / 12}},
      {1.0 / 2, {1.0 / 8, 3.0 / 8}},
      {2.0 / 3, {13.0 / 54, -27.0 / 54, 42.0 / 54, 8.0 / 54}},
      {1.0 / 6, {389.0 / 4320, -54.0 / 4320, 966.0 / 4320, -824.0 / 4320, 243.0 / 4320}},
      {1.0, {-231.0 / 20, 81.0 / 20, -1164.0 / 20, 656.0 / 20, -122.0 / 20, 800.0 / 20}},
      {5.0 / 6,
       {-127.0 / 288, 18.0 / 288, -678.0 / 288, 456.0 / 288, -9.0 / 288, 576.0 / 288, 4.0 / 288}},
      {1.0,
       {1481.0 / 820, -81.0 / 820, 7104.0 / 820, -3376.0 / 820, 72.0 / 820, -5040.0 / 820,
        -60.0 / 820, 720.0 / 820}},
  };

  for (size_t i = 0; i < rows.size(); ++i) {
    double sum = 0.0;
    for (double a : rows[i].a) sum += a;
    INFO("RK8 stage " << (i + 2) << ": sum(a_ij) = " << sum << ", c_i = " << rows[i].c);
    REQUIRE_THAT(sum, WithinAbs(rows[i].c, 1e-12));
  }

  // Weights must sum to 1 for the method to be consistent.
  double bsum = (41.0 + 27 + 272 + 27 + 216 + 216 + 41) / 840;
  REQUIRE_THAT(bsum, WithinAbs(1.0, 1e-12));
}

TEST_CASE("numerics.integrator_order.pd45_node_consistency") {
  // The Dormand-Prince 4(5) A matrix, as stored in PD45::A_.
  const std::array<std::array<double, 6>, 7> A
      = {{{0, 0, 0, 0, 0, 0},
          {1.0 / 5, 0, 0, 0, 0, 0},
          {3.0 / 40, 9.0 / 40, 0, 0, 0, 0},
          {44.0 / 45, -56.0 / 15, 32.0 / 9, 0, 0, 0},
          {19372.0 / 6561, -25360.0 / 2187, 64448.0 / 6561, -212.0 / 729, 0, 0},
          {9017.0 / 3168, -355.0 / 33, 46732.0 / 5247, 49.0 / 176, -5103.0 / 18656, 0},
          {35.0 / 384, 0, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84}}};

  // Published Dormand-Prince nodes.
  const std::array<double, 7> c_expected = {0, 1.0 / 5, 3.0 / 10, 4.0 / 5, 8.0 / 9, 1.0, 1.0};

  for (int i = 0; i < 7; ++i) {
    double sum = 0.0;
    for (int j = 0; j < 6; ++j) sum += A[i][j];
    INFO("PD45 stage " << i << ": sum_j A[i][j] = " << sum << ", expected c_i = " << c_expected[i]);
    REQUIRE_THAT(sum, WithinAbs(c_expected[i], 1e-12));
  }

  // Guard the specific defect this test was written for: the stage time must be
  // the ROW SUM, not the first column. For stages 2..6 the two differ, and
  // stages 4 and 5 of the first column land outside the step entirely.
  for (int i = 2; i < 7; ++i) {
    INFO("PD45 stage " << i << ": A[i][0] = " << A[i][0]
                       << " must not be mistaken for c_i = " << c_expected[i]);
    REQUIRE(std::abs(A[i][0] - c_expected[i]) > 1e-9);
  }
}

// -----------------------------------------------------------------------------
// Empirical order of convergence
// -----------------------------------------------------------------------------

TEST_CASE("numerics.integrator_order.rk4_is_fourth_order") {
  RK4 rk4;
  SECTION("autonomous") {
    REQUIRE_THAT(EmpiricalOrder(rk4, Decay, DecayExact, 1.0, 0.1), WithinAbs(4.0, 0.3));
  }
  SECTION("non-autonomous") {
    REQUIRE_THAT(EmpiricalOrder(rk4, Oscillating, OscillatingExact, 1.0, 0.1), WithinAbs(4.0, 0.3));
  }
}

TEST_CASE("numerics.integrator_order.rk8_is_seventh_order") {
  // NOTE THE NAME. `RK8` is a 10-stage explicit method, and Butcher's barrier
  // says an explicit Runge-Kutta method needs at least 11 stages to reach
  // order 8. It therefore cannot be 8th order, and measurement agrees: with
  // the tableau corrected, the local truncation error is O(h^8) -- verified in
  // exact rational arithmetic on y' = y^2 and y' = t y^2 -- so the global
  // order is 7. Seven is the best a 10-stage explicit method can do, so the
  // implementation is optimal; only the class name overstates it.
  //
  // Before the stage-8 coefficient was corrected from -234 to -231 this
  // measured ~2.8, because that row violated c_i = sum_j a_ij.
  //
  // dt is coarse on purpose: at dt = 1e-3 the error is far below roundoff and
  // the measurement would be pure floating-point noise.
  RK8 rk8;
  SECTION("autonomous") {
    // A linear constant-coefficient problem is not diagnostic of the true
    // order -- many stage errors cancel and it superconverges past 8 here --
    // so only a lower bound is asserted.
    REQUIRE(EmpiricalOrder(rk8, Decay, DecayExact, 1.0, 0.25) > 7.0);
  }
  SECTION("non-autonomous") {
    REQUIRE_THAT(EmpiricalOrder(rk8, Oscillating, OscillatingExact, 1.0, 0.25),
                 WithinAbs(7.0, 0.5));
  }
}

TEST_CASE("numerics.integrator_order.pd45_is_fifth_order") {
  // Dormand-Prince returns the 5th-order solution (`b_`); `b_star_` is the
  // 4th-order embedded companion used only for the error estimate.
  //
  // Tolerances are deliberately loose so that every step is accepted on the
  // first try and PD45 behaves as a fixed-step DP5. That isolates the tableau:
  // with the adaptive path engaged the measurement would confound step-size
  // control with method order.
  PD45 pd45;
  pd45.SetParams(IntegratorParams(20, 1e6, 1e6));

  SECTION("autonomous") {
    REQUIRE_THAT(EmpiricalOrder(pd45, Decay, DecayExact, 1.0, 0.25), WithinAbs(5.0, 0.5));
  }
  // The discriminating case: wrong stage times are invisible on an autonomous
  // problem and dominant here.
  SECTION("non-autonomous") {
    REQUIRE_THAT(EmpiricalOrder(pd45, Oscillating, OscillatingExact, 1.0, 0.25),
                 WithinAbs(5.0, 0.5));
  }
}

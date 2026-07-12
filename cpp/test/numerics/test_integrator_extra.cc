#include <lupnt/numerics/integrator.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

using namespace lupnt;
using namespace Catch::Matchers;

namespace {

  // --- Reference ODEs with closed-form solutions -------------------------------

  // y' = -y  =>  y(t) = y0 * exp(-(t - t0))
  ODE Decay = [](const Real /*t*/, const VecX& x) {
    VecX d(1);
    d(0) = -x(0);
    return d;
  };

  // y' = c  =>  y(t) = y0 + c*(t - t0)
  ODE Constant = [](const Real /*t*/, const VecX& x) {
    VecX d(1);
    d(0) = 3.0;
    return d;
  };

  const double kOmega = 2.0 * PI;  // harmonic-oscillator angular frequency

  // x'' = -w^2 x  =>  [x, v]' = [v, -w^2 x]; period T = 2*pi/w = 1
  ODE Harmonic = [](const Real /*t*/, const VecX& x) {
    VecX d(2);
    d(0) = x(1);
    d(1) = -kOmega * kOmega * x(0);
    return d;
  };

  VecX Vec1(double v) {
    VecX x(1);
    x(0) = v;
    return x;
  }

  // --- Fixed-step integrators vs analytic solutions ----------------------------

  TEST_CASE("numerics.integrator_extra.fixed_step_closed_form") {
    VecX x0 = Vec1(1.0);
    Real dt = 1e-3;

    SECTION("RK4 exponential decay matches exp(-t)") {
      RK4 rk4;
      VecX xf = rk4.Propagate(Decay, Real(0.0), Real(1.0), x0, dt);
      REQUIRE_THAT(xf(0).val(), WithinAbs(std::exp(-1.0), 1e-8));
    }

    SECTION("RK8 exponential decay matches exp(-t)") {
      RK8 rk8;
      VecX xf = rk8.Propagate(Decay, Real(0.0), Real(1.0), x0, dt);
      REQUIRE_THAT(xf(0).val(), WithinAbs(std::exp(-1.0), 1e-8));
    }

    SECTION("Constant derivative integrates exactly (RK4 is exact for linear-in-t)") {
      RK4 rk4;
      // y(2) = 1 + 3*2 = 7
      VecX xf = rk4.Propagate(Constant, Real(0.0), Real(2.0), x0, Real(0.25));
      REQUIRE_THAT(xf(0).val(), WithinAbs(7.0, 1e-9));
    }
  }

  // --- Adaptive integrators (RKF45, PD45) --------------------------------------

  // Mirrors the existing RK4/RK8 harmonic-oscillator test: after exactly one
  // period the oscillator returns to its initial state. Exercises the adaptive
  // Step() paths of IRKF/RKF45 and PD45.
  TEST_CASE("numerics.integrator_extra.adaptive_harmonic") {
    Real dt = 0.01;

    SECTION("RKF45 returns to initial state after one period") {
      RKF45 integ;
      Real t = 0;
      VecX x(2);
      x(0) = 1;
      x(1) = 0;
      for (int i = 0; i < static_cast<int>(1.0 / dt.val()); i++) {
        x = integ.Step(Harmonic, t, x, dt);
        t += dt;
      }
      REQUIRE_THAT(x(0).val(), WithinAbs(1.0, 1e-4));
      REQUIRE_THAT(x(1).val(), WithinAbs(0.0, 1e-4));
    }

    SECTION("PD45 returns to initial state after one period") {
      PD45 integ;
      Real t = 0;
      VecX x(2);
      x(0) = 1;
      x(1) = 0;
      for (int i = 0; i < static_cast<int>(1.0 / dt.val()); i++) {
        x = integ.Step(Harmonic, t, x, dt);
        t += dt;
      }
      REQUIRE_THAT(x(0).val(), WithinAbs(1.0, 1e-4));
      REQUIRE_THAT(x(1).val(), WithinAbs(0.0, 1e-4));
    }
  }

  TEST_CASE("numerics.integrator_extra.adaptive_decay") {
    VecX x0 = Vec1(1.0);

    SECTION("RKF45 Propagate reaches tf and matches exp(-t)") {
      RKF45 integ;
      IntegratorResult res = integ.PropagateEx(Decay, Real(0.0), Real(1.0), x0, Real(0.05));
      REQUIRE(res.reason == TerminationReason::ReachedTf);
      REQUIRE(res.steps > 0);
      REQUIRE_THAT(res.x(0).val(), WithinAbs(std::exp(-1.0), 1e-5));
      REQUIRE_THAT(res.t.val(), WithinAbs(1.0, 1e-12));
    }

    SECTION("PD45 Step loop matches exp(-t)") {
      PD45 integ;
      Real t = 0;
      Real dt = 0.02;
      VecX x = x0;
      for (int i = 0; i < 50; i++) {
        x = integ.Step(Decay, t, x, dt);
        t += dt;
      }
      REQUIRE_THAT(x(0).val(), WithinAbs(std::exp(-1.0), 1e-5));
    }
  }

  // --- PropagateEx: step count, termination reasons, edge cases ----------------

  TEST_CASE("numerics.integrator_extra.propagate_ex") {
    VecX x0 = Vec1(1.0);

    SECTION("ReachedTf reason with deterministic fixed step count") {
      RK4 rk4;
      // dt = 0.125 is exact in binary, so [0,1] is covered by exactly 8 steps
      // with no floating-point remainder step.
      IntegratorResult res = rk4.PropagateEx(Decay, Real(0.0), Real(1.0), x0, Real(0.125));
      REQUIRE(res.reason == TerminationReason::ReachedTf);
      REQUIRE(res.steps == 8);
      REQUIRE_THAT(res.t.val(), WithinAbs(1.0, 1e-12));
      REQUIRE_THAT(res.x(0).val(), WithinAbs(std::exp(-1.0), 1e-4));
    }

    SECTION("terminate_if predicate stops early with UserCondition") {
      RK4 rk4;
      IntegratorParams p;
      p.terminate_if = [](Real /*t*/, const VecX& x) { return x(0).val() < 0.5; };
      rk4.SetParams(p);
      IntegratorResult res = rk4.PropagateEx(Decay, Real(0.0), Real(10.0), x0, Real(0.01));
      REQUIRE(res.reason == TerminationReason::UserCondition);
      // exp(-t) crosses 0.5 at t = ln(2) ~ 0.6931
      REQUIRE_THAT(res.t.val(), WithinAbs(std::log(2.0), 0.02));
      REQUIRE(res.x(0).val() < 0.5);
    }

    SECTION("SetTerminateIf clears predicate with nullptr") {
      RK4 rk4;
      rk4.SetTerminateIf([](Real, const VecX& x) { return x(0).val() < 0.5; });
      rk4.SetTerminateIf(nullptr);
      IntegratorResult res = rk4.PropagateEx(Decay, Real(0.0), Real(1.0), x0, Real(0.1));
      REQUIRE(res.reason == TerminationReason::ReachedTf);
    }

    SECTION("Zero span returns immediately with zero steps") {
      RK4 rk4;
      IntegratorResult res = rk4.PropagateEx(Decay, Real(5.0), Real(5.0), x0, Real(0.1));
      REQUIRE(res.reason == TerminationReason::ReachedTf);
      REQUIRE(res.steps == 0);
      REQUIRE_THAT(res.x(0).val(), WithinAbs(1.0, 1e-12));
      REQUIRE_THAT(res.t.val(), WithinAbs(5.0, 1e-12));
    }

    SECTION("Backward propagation integrates in reverse time") {
      // y' = -y integrated backward from t0=1 (y=1) to tf=0 gives y(0)=e^{+1}.
      RK4 rk4;
      IntegratorResult res = rk4.PropagateEx(Decay, Real(1.0), Real(0.0), x0, Real(0.001));
      REQUIRE(res.reason == TerminationReason::ReachedTf);
      REQUIRE_THAT(res.t.val(), WithinAbs(0.0, 1e-9));
      REQUIRE_THAT(res.x(0).val(), WithinAbs(std::exp(1.0), 1e-5));
    }
  }

  // --- State-transition matrix via JacobianParallel path -----------------------

  TEST_CASE("numerics.integrator_extra.state_transition_matrix") {
    SECTION("Scalar decay STM equals exp(-dt)") {
      RK4 rk4;
      VecX x0 = Vec1(1.0);
      MatXd J;
      VecX xf = rk4.Propagate(Decay, Real(0.0), Real(1.0), x0, Real(1e-3), &J);
      REQUIRE(J.rows() == 1);
      REQUIRE(J.cols() == 1);
      REQUIRE_THAT(J(0, 0), WithinAbs(std::exp(-1.0), 1e-5));
      REQUIRE_THAT(xf(0).val(), WithinAbs(std::exp(-1.0), 1e-6));
    }

    SECTION("Harmonic oscillator STM over one period is identity") {
      RK4 rk4;
      VecX x0(2);
      x0(0) = 1.0;
      x0(1) = 0.0;
      MatXd J;
      rk4.Propagate(Harmonic, Real(0.0), Real(1.0), x0, Real(1e-3), &J);
      REQUIRE(J.rows() == 2);
      REQUIRE(J.cols() == 2);
      REQUIRE_THAT(J(0, 0), WithinAbs(1.0, 1e-3));
      REQUIRE_THAT(J(1, 1), WithinAbs(1.0, 1e-3));
      REQUIRE_THAT(J(0, 1), WithinAbs(0.0, 1e-3));
      REQUIRE_THAT(J(1, 0), WithinAbs(0.0, 1e-3));
    }
  }

  // --- Parameter validation and error branches ---------------------------------

  TEST_CASE("numerics.integrator_extra.param_validation") {
    SECTION("IntegratorParams rejects non-positive max_iter/tolerances") {
      REQUIRE_THROWS_AS(IntegratorParams(0, 1e-6, 1e-6), std::runtime_error);
      REQUIRE_THROWS_AS(IntegratorParams(20, 0.0, 1e-6), std::runtime_error);
      REQUIRE_THROWS_AS(IntegratorParams(20, 1e-6, -1.0), std::runtime_error);
      // Valid params do not throw.
      REQUIRE_NOTHROW(IntegratorParams(20, 1e-6, 1e-6));
    }

    SECTION("Propagate rejects non-positive step size") {
      RK4 rk4;
      VecX x0 = Vec1(1.0);
      REQUIRE_THROWS_AS(rk4.PropagateEx(Decay, Real(0.0), Real(1.0), x0, Real(-0.1)),
                        std::runtime_error);
      REQUIRE_THROWS_AS(rk4.Propagate(Decay, Real(0.0), Real(1.0), State(x0), Real(0.0)),
                        std::runtime_error);
    }
  }

}  // namespace

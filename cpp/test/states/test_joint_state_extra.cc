#include <lupnt/dynamics/clock_dynamics.h>
#include <lupnt/states/joint_state.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  ClockState3 MakeClock(double b, double d, double dr) {
    ClockState3 c;
    c.b() = b;
    c.d() = d;
    c.dr() = dr;
    return c;
  }
}  // namespace

TEST_CASE("states.joint_state_extra.multi_add_composition") {
  JointState joint;
  joint.Add(MakeClock(1.0, 2.0, 3.0), MakePtr<ClockDynamics>());
  joint.Add(MakeClock(4.0, 5.0, 6.0), MakePtr<ClockDynamics>());

  REQUIRE(joint.GetStateSize() == 6);
  REQUIRE(joint.GetParamSize() == 0);
  REQUIRE(joint.GetStmSize() == 6);

  // GetState concatenates the two sub-states in insertion order.
  State x = joint.GetState();
  REQUIRE(x.size() == 6);
  REQUIRE_THAT(x(0).val(), WithinAbs(1.0, 1e-9));
  REQUIRE_THAT(x(2).val(), WithinAbs(3.0, 1e-9));
  REQUIRE_THAT(x(3).val(), WithinAbs(4.0, 1e-9));
  REQUIRE_THAT(x(5).val(), WithinAbs(6.0, 1e-9));
}

TEST_CASE("states.joint_state_extra.fixed_params_excluded_from_stm") {
  JointState joint;

  VecX pvec(2);
  pvec << 10.0, 20.0;
  ParamState params(pvec, {"est", "fixed"});

  // One estimated, one fixed parameter.
  joint.Add(MakeClock(1.0, 2.0, 3.0), MakePtr<ClockDynamics>(), nullptr, params,
            {ESTIMATED, FIXED});

  REQUIRE(joint.GetParamSize() == 2);
  REQUIRE(joint.GetEstimatedParamSize() == 1);
  REQUIRE(joint.GetFixedParamSize() == 1);
  REQUIRE(joint.GetConsideredParamSize() == 0);

  // Full size = state(3) + params(2) = 5; STM size drops the fixed param: 3 + 1 = 4.
  REQUIRE(joint.GetSize() == 5);
  REQUIRE(joint.GetStmSize() == 4);

  // GetStmState carries the state plus only the estimated/considered params.
  State xs = joint.GetStmState();
  REQUIRE(xs.size() == 4);
  REQUIRE_THAT(xs(0).val(), WithinAbs(1.0, 1e-9));   // clock bias
  REQUIRE_THAT(xs(3).val(), WithinAbs(10.0, 1e-9));  // estimated param, fixed omitted

  // GetState still exposes all parameters.
  ParamState p = joint.GetParams();
  REQUIRE(p.size() == 2);
  REQUIRE_THAT(p(1).val(), WithinAbs(20.0, 1e-9));
}

TEST_CASE("states.joint_state_extra.dynamics_function_shape") {
  JointState joint;
  joint.Add(MakeClock(1.0, 2.0, 3.0), MakePtr<ClockDynamics>());

  FilterDynamicsFunction f = joint.GetDynamicsFunction();
  REQUIRE(static_cast<bool>(f));

  MatXd Phi;
  State xf = f(joint.GetStmState(), 0.0, 30.0, nullptr, &Phi);

  REQUIRE(xf.size() == joint.GetStmSize());
  REQUIRE(Phi.rows() == joint.GetStmSize());
  REQUIRE(Phi.cols() == joint.GetStmSize());
}

TEST_CASE("states.joint_state_extra.process_noise_function") {
  JointState joint;

  VecX pvec(1);
  pvec << 5.0;
  ParamState params(pvec, {"gm"});

  // Estimated parameter with a finite first-order Gauss-Markov time constant.
  double tau = 100.0, q = 0.5;
  joint.Add(MakeClock(1.0, 2.0, 3.0), MakePtr<ClockDynamics>(), nullptr, params, {ESTIMATED},
            {{tau, q}});

  ProcessNoiseFunction pn = joint.GetProcessNoiseFunction();
  REQUIRE(static_cast<bool>(pn));

  double dt = 10.0;
  MatXd Q = pn(joint.GetStmState(), 0.0, dt);

  REQUIRE(Q.rows() == joint.GetStmSize());  // 3 + 1
  REQUIRE(Q.cols() == joint.GetStmSize());
  REQUIRE(Q.isApprox(Q.transpose(), 1e-12));

  // The Gauss-Markov parameter block carries the closed-form steady growth
  // q*tau/2*(1 - exp(-2*dt/tau)) in its diagonal slot (last index).
  double expected = q * tau / 2.0 * (1.0 - std::exp(-2.0 * dt / tau));
  int last = joint.GetStmSize() - 1;
  REQUIRE_THAT(Q(last, last), WithinRel(expected, 1e-9));
  REQUIRE(Q(last, last) > 0.0);
}

#include <lupnt/dynamics/parameter_dynamics.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

// ParameterDynamics models constant parameters: propagation is the identity map
// and its state-transition matrix is the identity, regardless of the time span.
TEST_CASE("dynamics.parameter_dynamics") {
  VecX values(3);
  values << 1.5, -2.0, 42.0;
  ParamState p(values, {"srp", "drag", "bias"});
  ParameterDynamics dyn;

  SECTION("propagation leaves constant parameters unchanged") {
    State xf = dyn.Propagate(p, 0.0, 1234.5);
    VecXd xf_d = xf.cast<double>();
    REQUIRE(xf_d.size() == values.size());
    for (int i = 0; i < values.size(); ++i)
      REQUIRE_THAT(xf_d(i), WithinAbs(values(i).val(), 1e-12));
  }

  SECTION("state-transition matrix is the identity") {
    MatXd stm;
    State xf = dyn.Propagate(p, 100.0, 900.0, nullptr, &stm);
    REQUIRE(stm.rows() == values.size());
    REQUIRE(stm.cols() == values.size());
    MatXd expected = MatXd::Identity(values.size(), values.size());
    for (int i = 0; i < stm.rows(); ++i)
      for (int j = 0; j < stm.cols(); ++j)
        REQUIRE_THAT(stm(i, j), WithinAbs(expected(i, j), 1e-12));
    // And the state itself is still unchanged.
    VecXd xf_d = xf.cast<double>();
    for (int i = 0; i < values.size(); ++i)
      REQUIRE_THAT(xf_d(i), WithinAbs(values(i).val(), 1e-12));
  }

  SECTION("reports the Params state type") { REQUIRE(dyn.GetStateType() == "Params"); }
}

#include <lupnt/agents/satellite.h>
#include <lupnt/dynamics/surface_dynamics.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("agents.satellite") {
  Satellite sat;
  sat.SetName("sat");
  sat.SetTime(10.0);
  sat.SetFrequency(2.0);
  sat.SetDynamics(MakePtr<StaticDynamics>());

  Cart6 state(Vec6(1.0, 2.0, 3.0, 0.1, 0.2, 0.3), Frame::GCRF);
  sat.SetState(state);

  REQUIRE(sat.GetName() == "sat");
  REQUIRE_THAT(sat.GetTime().val(), Catch::Matchers::WithinAbs(10.0, epsilon));
  REQUIRE_THAT(sat.GetFrequency().val(), Catch::Matchers::WithinAbs(2.0, epsilon));
  REQUIRE(sat.GetState().GetFrame() == Frame::GCRF);
  REQUIRE(sat.GetState().isApprox(state, epsilon));
  REQUIRE(sat.GetStateAt(20.0).isApprox(state, epsilon));
  REQUIRE_NOTHROW(sat.Log(10.0));
  REQUIRE_NOTHROW(sat.LogCesium());
}

#include <lupnt/agents/rover.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;

// The refactored `Rover` is a thin, app-hosting agent: its constructor reads
// `name` / `frequency` and attaches the `application:` block, but it integrates
// no dynamics itself -- the hosted app (e.g. `SurfaceRoverNavApp`) drives the
// truth state each epoch via `SetState`, and `GetStateAt` returns it verbatim.
// The 2D unicycle physics that this test used to exercise through the (removed)
// monolithic-Rover propagation path is now covered directly by
// `cpp/test/dynamics/test_surface_dynamics.cc`.
TEST_CASE("agents.rover") {
  Config config = YAML::Load(R"(
name: rover
frequency: 1.0
)");
  Rover rover(config);
  REQUIRE(rover.GetName() == "rover");

  // The hosted application sets the truth state each epoch; GetStateAt returns
  // it verbatim (the query time is ignored -- no dynamics are integrated).
  Vec6 truth;
  truth << 1.0, 2.0, 3.0, 0.1, 0.2, 0.3;
  rover.SetState(Cart6(truth, Frame::MOON_PA));

  Cart6 state = rover.GetStateAt(123.0);
  REQUIRE(state.GetFrame() == Frame::MOON_PA);
  REQUIRE_THAT(state(0).val(), Catch::Matchers::WithinAbs(1.0, 1e-9));
  REQUIRE_THAT(state(2).val(), Catch::Matchers::WithinAbs(3.0, 1e-9));
  REQUIRE_THAT(state(5).val(), Catch::Matchers::WithinAbs(0.3, 1e-9));
}

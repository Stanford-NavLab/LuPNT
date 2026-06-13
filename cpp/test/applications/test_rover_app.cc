#include <lupnt/applications/rover_app.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("applications.rover_app") {
  Config config = YAML::Load(R"(
name: rover_app
dynamics:
  class: SurfaceDynamics2D
filter:
  class: EKF
initial_state:
  sigma_r: 1.0
  sigma_theta: 1.0
trajectory:
  R: 10.0
  v: 1.0
)");
  RoverApp app(config);
  REQUIRE(app.GetName() == "rover_app");
  Config bad_config = YAML::Load("name: bad\nfilter: {class: EKF}\n");
  REQUIRE_THROWS(RoverApp(bad_config));
}

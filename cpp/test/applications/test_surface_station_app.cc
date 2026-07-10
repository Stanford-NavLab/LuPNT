#include <lupnt/applications/lunar_station/surface_station_app.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("applications.surface_station_app") {
  Config config = YAML::Load("name: station_app\nfrequency: 0.5\n");
  SurfaceStationApp app(config);
  REQUIRE(app.GetName() == "station_app");
  REQUIRE_THAT(app.GetFrequency().val(), Catch::Matchers::WithinAbs(0.5, epsilon));
  REQUIRE_NOTHROW(app.Setup());
  REQUIRE_NOTHROW(app.Step(2.0));
}

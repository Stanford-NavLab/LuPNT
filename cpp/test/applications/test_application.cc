#include <lupnt/applications/application.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

namespace {
  class TestApplication : public Application {
  public:
    using Application::Application;
    void Step(Real t) override {
      last_t = t;
      ++count;
    }

    Real last_t = 0.0;
    int count = 0;
  };
}  // namespace

TEST_CASE("applications.application") {
  Config config = YAML::Load("name: app\nfrequency: 4.0\n");
  TestApplication app(config);
  REQUIRE(app.GetName() == "app");
  REQUIRE_THAT(app.GetFrequency().val(), Catch::Matchers::WithinAbs(4.0, epsilon));

  app.SetName("renamed");
  app.SetFrequency(0.0);
  REQUIRE(app.GetName() == "renamed");
  REQUIRE_THAT(app.GetFrequency().val(), Catch::Matchers::WithinAbs(0.0, epsilon));
  REQUIRE_NOTHROW(app.Setup());
  app.Step(12.0);
  REQUIRE(app.count == 1);
  REQUIRE_THAT(app.last_t.val(), Catch::Matchers::WithinAbs(12.0, epsilon));
}

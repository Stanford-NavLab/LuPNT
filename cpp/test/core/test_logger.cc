#include <lupnt/core/logger.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("core.logger") {
  SECTION("Logger accepts all levels without throwing") {
    Logger::SetLogLevel(Logger::DEBUG);

    REQUIRE_NOTHROW(Logger::Debug("debug message", "test", 1.0));
    REQUIRE_NOTHROW(Logger::Info("info message", "test", 1.0));
    REQUIRE_NOTHROW(Logger::Warn("warning message", "test", 1.0));
    REQUIRE_NOTHROW(Logger::Error("error message", "test", 1.0));

    Logger::SetLogLevel(Logger::INFO);
  }

  SECTION("Logger creates a usable progress bar") {
    auto bar = Logger::GetProgressBar(2, "working", "test", 0.0);
    bar->SetMaxUpdateFreq(0.0);

    bar->Update();
    REQUIRE_FALSE(bar->IsDone());

    bar->Update(2);
    REQUIRE(bar->IsDone());
  }
}

#include <lupnt/core/progress_bar.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("core.progress_bar") {
  SECTION("ProgressBar tracks completion and reset") {
    ProgressBar bar(3);
    bar.SetMaxUpdateFreq(0.0);
    bar.SetLeave(false);

    REQUIRE_FALSE(bar.IsDone());

    bar.Update();
    bar.Update(2);
    REQUIRE_FALSE(bar.IsDone());

    bar.Update(3);
    REQUIRE(bar.IsDone());

    bar.Reset();
    REQUIRE_FALSE(bar.IsDone());

    bar.Finish();
    REQUIRE(bar.IsDone());
  }
}

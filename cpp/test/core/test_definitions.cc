#include <lupnt/core/definitions.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("core.definitions") {
  SECTION("global LuPNT epoch can be set and restored") {
    const Real original = GetLupntEpoch();

    SetLupntEpoch(12345.5);
    REQUIRE_THAT(GetLupntEpoch().val(), WithinAbs(12345.5, epsilon));

    SetLupntEpoch(original);
    REQUIRE_THAT(GetLupntEpoch().val(), WithinAbs(original.val(), epsilon));
  }

  SECTION("pointer aliases construct shared objects") {
    auto value = MakePtr<int>(42);

    REQUIRE(value);
    REQUIRE(*value == 42);
  }
}

#include <lupnt/core/constants.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <iostream>
#include <vector>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("core.enums") {
  auto vals = enum_values<BodyId>();
  REQUIRE(vals.size() > 20);
  for (auto naif_id : vals) {
    auto naif_id_str = enum_name(naif_id);
    auto naif_id_val = enum_cast<BodyId>(naif_id_str).value();
    REQUIRE(naif_id == naif_id_val);
  }
}

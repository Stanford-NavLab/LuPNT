#include <lupnt/interfaces/matplot.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("interfaces.matplot") {
  SUCCEED("matplot interface declarations are compile-checked by including the public header");
}

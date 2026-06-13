#include <lupnt/interfaces/python.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("interfaces.python") {
  SUCCEED("Python interface is compile-gated by LUPNT_WITH_PYTHON in C++ tests");
}

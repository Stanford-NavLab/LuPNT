#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <iostream>

// Demonstrate some basic assertions.
TEST_CASE("HelloTest", "[BasicAssertions]") {
  // Expect two strings not to be equal.
  REQUIRE_THAT("hello", !Catch::Matchers::Equals("world"));
  // Expect equality.
  REQUIRE(7 * 6 == 42);

  // print environment variable LUPNT_DATA_PATH
  char *env = std::getenv("LUPNT_DATA_PATH");
  if (env != NULL) {
    std::cout << "LUPNT_DATA_PATH=" << env << '\n';
  } else {
    std::cout << "LUPNT_DATA_PATH is not set\n";
  }
}
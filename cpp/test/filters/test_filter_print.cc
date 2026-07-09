#include <lupnt/numerics/filters/filter_print.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <sstream>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("filters.filter_print") {
  std::ostringstream capture;
  auto* old_buf = std::cout.rdbuf(capture.rdbuf());

  PrintEKFProgressHeaderPVC();
  PrintEKFProgressPVC(120.0, 1.0, 2.0, 3.0);
  PrintEKFProgressHeaderPVC(2);
  PrintEKFProgressPVC(120.0, VecXd::Ones(8), 2);
  PrintEstimationStatistics(VecXd::Ones(4), MatXd::Ones(8, 4), 0.5, 2);

  std::cout.rdbuf(old_buf);
  std::string output = capture.str();
  REQUIRE_THAT(output, ContainsSubstring("Time [min]"));
  REQUIRE_THAT(output, ContainsSubstring("Sat1"));
  REQUIRE_THAT(output, ContainsSubstring("Simulation Statistics"));
}

#include <cspice/SpiceUsr.h>
#include <lupnt/lupnt.h>
#include <lupnt/version.h>
#include <matplot/matplot.h>
#include <omp.h>

#include <Eigen/Core>
#include <autodiff/forward/dual.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <string>

TEST_CASE("LuPNT") {
  using namespace lupnt;

  double expected = 398600.435507;
  REQUIRE_THAT(GM_EARTH, Catch::Matchers::WithinRel(expected, 1e-3)
                             && Catch::Matchers::WithinAbs(expected, 1e-6));
}

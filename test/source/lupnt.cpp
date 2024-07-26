#include <SpiceUsr.h>
#include <doctest/doctest.h>
#include <lupnt/lupnt.h>
#include <lupnt/version.h>
#include <matplot/matplot.h>
#include <omp.h>

#include <Eigen/Core>
#include <autodiff/forward/dual.hpp>
#include <string>

TEST_CASE("LuPNT") {
  using namespace lupnt;

  CHECK(GM_EARTH == doctest::Approx(398600.435507));
}

TEST_CASE("LuPNT version") {
  static_assert(std::string_view(LUPNT_VERSION) == std::string_view("1.0"));
  CHECK(std::string(LUPNT_VERSION) == std::string("1.0"));
}

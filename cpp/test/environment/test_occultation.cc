#include <lupnt/environment/occultation.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("environment.occultation") {
  SECTION("empty body list reports all links visible") {
    VecXd heights(0);
    VecXd min_elevation(0);

    auto vis
        = Occultation::ComputeOccultation(0.0, Vec3::Zero(), Vec3(1.0, 0.0, 0.0), Frame::ICRF,
                                          Frame::ICRF, {}, heights, min_elevation, false, false);

    REQUIRE(vis.size() == 1);
    REQUIRE(vis["all"]);
  }

  SECTION("batched overload validates row counts") {
    MatX3 r1 = MatX3::Zero(2, 3);
    MatX3 r2 = MatX3::Zero(3, 3);
    VecXd heights(0);
    VecXd min_elevation(0);

    REQUIRE_THROWS(Occultation::ComputeOccultation(0.0, r1, r2, Frame::ICRF, Frame::ICRF, {},
                                                   heights, min_elevation, false, false));
  }
}

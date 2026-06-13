#include <lupnt/conversions/frame_converter.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("conversions.frame_converter") {
  SECTION("frame center lookup returns canonical bodies") {
    REQUIRE(GetFrameCenter(Frame::GCRF) == BodyId::EARTH);
    REQUIRE(GetFrameCenter(Frame::ITRF) == BodyId::EARTH);
    REQUIRE(GetFrameCenter(Frame::MOON_CI) == BodyId::MOON);
    REQUIRE(GetFrameCenter(Frame::ICRF) == BodyId::SOLAR_SYSTEM_BARYCENTER);
  }

  SECTION("same-frame conversions return the input unchanged") {
    Vec6 rv(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    Vec3 r(1.0, 2.0, 3.0);

    Vec6 rv_out = ConvertFrame(0.0, rv, Frame::GCRF, Frame::GCRF);
    Vec3 r_out = ConvertFrame(0.0, r, Frame::GCRF, Frame::GCRF);

    for (int i = 0; i < 6; ++i) REQUIRE_THAT(rv_out(i).val(), WithinAbs(rv(i).val(), epsilon));
    for (int i = 0; i < 3; ++i) REQUIRE_THAT(r_out(i).val(), WithinAbs(r(i).val(), epsilon));
  }

  SECTION("matrix overloads preserve row counts") {
    MatX6 rv(2, 6);
    rv.row(0) = Vec6(1, 2, 3, 4, 5, 6).transpose();
    rv.row(1) = Vec6(7, 8, 9, 10, 11, 12).transpose();

    MatX6 out = ConvertFrame(0.0, rv, Frame::GCRF, Frame::GCRF);

    REQUIRE(out.rows() == 2);
    REQUIRE(out.cols() == 6);
    REQUIRE_THAT(out(1, 5).val(), WithinAbs(12.0, epsilon));
  }
}

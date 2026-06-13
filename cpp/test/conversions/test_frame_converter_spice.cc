#include <lupnt/conversions/frame_converter_spice.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("conversions.frame_converter_spice") {
  Vec6 rv;
  rv << 1.0, 2.0, 3.0, 0.1, 0.2, 0.3;
  Vec6 same = spice::ConvertFrameSpice(0.0, rv, Frame::GCRF, Frame::GCRF);
  REQUIRE(same.isApprox(rv, epsilon));

  Vec3 r(1.0, 2.0, 3.0);
  Vec3 r_same = spice::ConvertFrameSpice(0.0, r, Frame::GCRF, Frame::GCRF);
  REQUIRE(r_same.isApprox(r, epsilon));

  VecX times(2);
  times << 0.0, 1.0;
  MatX6 rv_history = spice::ConvertFrameSpice(times, rv, Frame::GCRF, Frame::GCRF);
  REQUIRE(rv_history.rows() == 2);
  REQUIRE(rv_history.row(0).transpose().isApprox(rv, epsilon));

  MatX6 rv_in(1, 6);
  rv_in.row(0) = rv.transpose();
  REQUIRE_THROWS(spice::ConvertFrameSpice(times, rv_in, Frame::GCRF, Frame::GCRF));

  Cart6 state(rv, Frame::GCRF);
  Cart6 state_same = spice::ConvertFrameSpice(0.0, state, Frame::GCRF);
  REQUIRE(state_same.GetFrame() == Frame::GCRF);
  REQUIRE(state_same.isApprox(state, epsilon));
}

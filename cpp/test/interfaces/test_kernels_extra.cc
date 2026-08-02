#include <lupnt/interfaces/kernels.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "../utils.cc"
#include "lupnt/conversions/frame_converter.h"
#include "lupnt/interfaces/spice.h"

using namespace lupnt;
using Catch::Approx;

// Free-function ephemeris accessors in kernels.cc. These wrap spice::GetBodyPosVel
// plus a frame conversion; they are distinct from (and do not touch) the
// segfaulting spice::ConvertTime / et2utc_c paths.

namespace {
  Real Epoch() {
    spice::LoadSpiceKernel();
    return GetLupntEpoch();
  }
  constexpr double kMoonDistMinKm = 350000.0;
  constexpr double kMoonDistMaxKm = 410000.0;
}  // namespace

TEST_CASE("interfaces.kernels_extra.body_pos_vel_earth_moon") {
  Real t = Epoch();
  Vec6 rv = GetBodyPosVel(t, BodyId::EARTH, BodyId::MOON, Frame::GCRF);  // meters
  double dist_km = rv.head(3).norm().val() / 1000.0;
  REQUIRE(dist_km > kMoonDistMinKm);
  REQUIRE(dist_km < kMoonDistMaxKm);

  // Antisymmetry between swapped center/target.
  Vec6 rv_rev = GetBodyPosVel(t, BodyId::MOON, BodyId::EARTH, Frame::GCRF);
  REQUIRE(rv.head(3).isApprox(-rv_rev.head(3), 1e-9));
}

TEST_CASE("interfaces.kernels_extra.two_arg_overload_uses_frame_center") {
  Real t = Epoch();
  // The (target, frame) overload implies the frame's natural center; for GCRF
  // that is Earth, so it must match the explicit Earth-centered call.
  Vec6 rv_2arg = GetBodyPosVel(t, BodyId::MOON, Frame::GCRF);
  Vec6 rv_4arg = GetBodyPosVel(t, BodyId::EARTH, BodyId::MOON, Frame::GCRF);
  REQUIRE(rv_2arg.isApprox(rv_4arg, 1e-9));
}

TEST_CASE("interfaces.kernels_extra.body_pos_matches_pos_vel") {
  Real t = Epoch();
  Vec3 r = GetBodyPos(t, BodyId::EARTH, BodyId::MOON, Frame::GCRF);
  Vec6 rv = GetBodyPosVel(t, BodyId::EARTH, BodyId::MOON, Frame::GCRF);
  REQUIRE(r.isApprox(rv.head(3), 1e-9));
}

TEST_CASE("interfaces.kernels_extra.unit_system_scaling") {
  Real t = Epoch();
  Vec6 rv_si = GetBodyPosVel(t, BodyId::EARTH, BodyId::MOON, Frame::GCRF);
  Vec6 rv_km = GetBodyPosVel(t, BodyId::EARTH, BodyId::MOON, Frame::GCRF, KM_S_KG_UNITS);
  REQUIRE(rv_km.head(3).isApprox(rv_si.head(3) / 1000.0, 1e-9));
  REQUIRE(rv_km.tail(3).isApprox(rv_si.tail(3) / 1000.0, 1e-9));
}

TEST_CASE("interfaces.kernels_extra.vectorized_overload") {
  Real t = Epoch();
  VecX times(3);
  times << t, t + 600.0, t + 1200.0;
  MatX6 rvs = GetBodyPosVel(times, BodyId::EARTH, BodyId::MOON, Frame::GCRF);

  REQUIRE(rvs.rows() == 3);
  REQUIRE(rvs.cols() == 6);
  for (int i = 0; i < 3; i++) {
    Vec6 rv = GetBodyPosVel(times(i), BodyId::EARTH, BodyId::MOON, Frame::GCRF);
    REQUIRE(rvs.row(i).transpose().isApprox(rv, 1e-9));
  }
}

TEST_CASE("interfaces.kernels_extra.tt_tdb_difference_not_implemented") {
  // GetTtTdbDifference is a declared-but-unimplemented stub; it must throw.
  REQUIRE_THROWS(GetTtTdbDifference(0.0));
}

// The vectorized GetLunarMantleData must agree row-for-row with the scalar
// overload it delegates to. It previously divided the scalar result by M_KM a
// second time, so every vectorized libration angle came out 1000x too small.
// The bug was latent -- nothing in-tree called the vector overload -- but it is
// public API in kernels.h.
TEST_CASE("interfaces.kernels_extra.lunar_mantle_vector_matches_scalar") {
  Real t = Epoch();

  VecX ts(3);
  ts << t, t + 3600.0, t + 7200.0;

  MatX6 rows = GetLunarMantleData(ts);
  for (int i = 0; i < ts.size(); ++i) {
    Vec6 one = GetLunarMantleData(Real(ts(i)));
    INFO("row " << i << " vector=" << rows(i, 0).val() << " scalar=" << one(0).val());
    for (int j = 0; j < 6; ++j) {
      REQUIRE(rows(i, j).val() == Approx(one(j).val()).epsilon(1e-12));
    }
  }

  // The libration angles are radians of order unity, not milliradians: a
  // second division by M_KM would put them near 1e-3 of this.
  REQUIRE(std::abs(rows(0, 0).val()) > 1e-3);
}

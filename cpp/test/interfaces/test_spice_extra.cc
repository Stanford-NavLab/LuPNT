#include <lupnt/interfaces/spice.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "../utils.cc"

using namespace lupnt;
using Catch::Approx;

// NOTE: this file deliberately exercises only the SPICE runtime helpers that
// route through str2et_c / spkezr_c / spkpos_c / sxform_c / pxform_c. It avoids
// the unitim_c / et2utc_c code paths (ConvertTime, StringToTai, TAItoStringUTC,
// TDBtoStringUTC), which segfault in this build (a pre-existing issue also
// noted in test_spice.cc / test_spice_interface.cc).

namespace {
  // A fixed, well-defined epoch (matches the value asserted in
  // test_spice_interface.cc's disabled StringToTdb case).
  Real Epoch() {
    spice::LoadSpiceKernel();
    return spice::StringToTdb("2023-04-15 00:00:00 TDB");
  }
  // Bounds on the Earth-Moon distance [km] (perigee ~356500, apogee ~406700).
  constexpr double kMoonDistMinKm = 350000.0;
  constexpr double kMoonDistMaxKm = 410000.0;
}  // namespace

TEST_CASE("interfaces.spice_extra.string_to_tdb") {
  Real et = Epoch();
  // 2023-04-15 00:00:00 TDB == 734788800 s past J2000 (TDB).
  REQUIRE(et.val() == Approx(734788800.0).margin(1.0));
}

TEST_CASE("interfaces.spice_extra.body_pos_vel_earth_moon") {
  Real et = Epoch();
  Vec6 rv = spice::GetBodyPosVel(et, BodyId::EARTH, BodyId::MOON);
  double dist = rv.head(3).norm().val();  // km
  REQUIRE(dist > kMoonDistMinKm);
  REQUIRE(dist < kMoonDistMaxKm);

  // Antisymmetry: r(Earth->Moon) == -r(Moon->Earth).
  Vec6 rv_rev = spice::GetBodyPosVel(et, BodyId::MOON, BodyId::EARTH);
  REQUIRE(rv.head(3).isApprox(-rv_rev.head(3), 1e-9));
  REQUIRE(rv.tail(3).isApprox(-rv_rev.tail(3), 1e-9));

  // A body relative to itself is the zero state.
  Vec6 rv_self = spice::GetBodyPosVel(et, BodyId::EARTH, BodyId::EARTH);
  REQUIRE(rv_self.norm().val() == Approx(0.0).margin(1e-6));
}

TEST_CASE("interfaces.spice_extra.body_pos_vel_vectorized") {
  Real et = Epoch();
  VecX times(3);
  times << et, et + 3600.0, et + 7200.0;
  MatX6 rvs = spice::GetBodyPosVel(times, BodyId::EARTH, BodyId::MOON);

  REQUIRE(rvs.rows() == 3);
  REQUIRE(rvs.cols() == 6);

  // Each row must agree with the per-epoch scalar overload.
  for (int i = 0; i < 3; i++) {
    Vec6 rv = spice::GetBodyPosVel(times(i), BodyId::EARTH, BodyId::MOON);
    REQUIRE(rvs.row(i).transpose().isApprox(rv, 1e-9));
  }
}

TEST_CASE("interfaces.spice_extra.body_pos_spice_consistency") {
  Real et = Epoch();

  // spkpos_c (GetBodyPosSpice) vs spkezr_c (GetBodyPosVelSpice): position agrees.
  Vec3d r = spice::GetBodyPosSpice(et, BodyId::EARTH, BodyId::MOON);
  Vec6d rv = spice::GetBodyPosVelSpice(et, BodyId::EARTH, BodyId::MOON);
  REQUIRE(r.isApprox(rv.head(3), 1e-9));

  double dist = r.norm();
  REQUIRE(dist > kMoonDistMinKm);
  REQUIRE(dist < kMoonDistMaxKm);

  // Direct SPICE lookup vs the Chebyshev-cached GetBodyPosVel: agree to ~meters.
  Vec6 rv_cheby = spice::GetBodyPosVel(et, BodyId::EARTH, BodyId::MOON);
  REQUIRE((rv.head(3) - rv_cheby.head(3).cast<double>()).norm() < 1.0);  // < 1 km
}

TEST_CASE("interfaces.spice_extra.frame_conversion_matrix_orthonormal") {
  Real et = Epoch();
  Mat6d M = spice::GetFrameConversionMat(et, "J2000", "ITRF93");
  REQUIRE(M.rows() == 6);
  REQUIRE(M.cols() == 6);

  // Top-left 3x3 is a proper rotation: R R^T = I, det = +1.
  Mat3d R = M.topLeftCorner(3, 3);
  REQUIRE((R * R.transpose()).isApprox(Mat3d::Identity(), 1e-9));
  REQUIRE(R.determinant() == Approx(1.0).margin(1e-9));

  // Forward then inverse frame map recovers identity (6x6, includes the
  // angular-velocity coupling block).
  Mat6d Minv = spice::GetFrameConversionMat(et, "ITRF93", "J2000");
  REQUIRE((M * Minv).isApprox(Mat6d::Identity(), 1e-9));
}

TEST_CASE("interfaces.spice_extra.planet_orientation_earth") {
  Real et = Epoch();
  // Earth's IAU rotation pole sits essentially at the J2000 celestial pole:
  // alpha0 ~ 0, delta0 ~ 90 deg.
  Vec3d o = spice::GetPlanetOrientation(BodyId::EARTH, et);
  REQUIRE(std::isfinite(o(0)));
  REQUIRE(std::isfinite(o(1)));
  REQUIRE(std::isfinite(o(2)));
  REQUIRE(o(1) == Approx(PI / 2.0).margin(1e-2));  // delta0 ~ 90 deg (89.87 at this epoch)

  // Moon orientation angles are finite and the prime-meridian angle is a
  // well-defined value.
  Vec3d om = spice::GetPlanetOrientation(BodyId::MOON, et);
  REQUIRE(std::isfinite(om(0)));
  REQUIRE(std::isfinite(om(1)));
  REQUIRE(std::isfinite(om(2)));
}

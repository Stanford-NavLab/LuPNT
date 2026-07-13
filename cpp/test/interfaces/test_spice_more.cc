#include <lupnt/interfaces/spice.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "../utils.cc"

using namespace lupnt;
using Catch::Approx;

// Deepens spice.cc coverage beyond test_spice.cc / test_spice_extra.cc, staying
// on the str2et_c / spkezr_c / spkpos_c / sxform_c / pxform_c / bod*_c runtime
// paths and avoiding the unitim_c / et2utc_c paths (ConvertTime / StringToTai /
// TAItoStringUTC / TDBtoStringUTC), which segfault in this build.

namespace {
  Real Epoch() {
    spice::LoadSpiceKernel();
    return spice::StringToTdb("2023-04-15 00:00:00 TDB");
  }
  constexpr double kMoonDistMinKm = 350000.0;
  constexpr double kMoonDistMaxKm = 410000.0;
  constexpr double kAuKm = 149597870.7;
}  // namespace

TEST_CASE("interfaces.spice_more.naif_id_name_roundtrip") {
  spice::LoadSpiceKernel();

  // Name <-> ID round trips for several bodies.
  REQUIRE(spice::GetNaifId("SUN") == static_cast<int>(BodyId::SUN));
  REQUIRE(spice::GetNaifId("MARS") == static_cast<int>(BodyId::MARS));
  REQUIRE(spice::GetNaifName(static_cast<int>(BodyId::SUN)) == "SUN");
  REQUIRE(spice::GetNaifName(static_cast<int>(BodyId::EARTH)) == "EARTH");

  REQUIRE(spice::HasNaifBody("SUN"));
  REQUIRE(spice::HasNaifBody("MARS"));

  // GetNaifName returns the numeric ID as a string when no name is defined
  // (it does not throw).
  REQUIRE(spice::GetNaifName(987654321) == "987654321");

  // GetNaifId throws on an unresolvable body name.
  REQUIRE_THROWS(spice::GetNaifId("__LUPNT_UNKNOWN_BODY__"));
}

TEST_CASE("interfaces.spice_more.body_pos_vel_spice_sun_au") {
  Real et = Epoch();

  // Direct spkezr lookup of the Sun relative to Earth: ~1 AU.
  Vec6d rv = spice::GetBodyPosVelSpice(et, BodyId::EARTH, BodyId::SUN);
  double dist = rv.head(3).norm();
  REQUIRE(dist > 0.98 * kAuKm);
  REQUIRE(dist < 1.02 * kAuKm);

  // Earth's orbital speed about the Sun is ~30 km/s.
  double speed = rv.tail(3).norm();
  REQUIRE(speed > 25.0);
  REQUIRE(speed < 35.0);

  // A body relative to itself is the zero state.
  Vec6d rv_self = spice::GetBodyPosVelSpice(et, BodyId::EARTH, BodyId::EARTH);
  REQUIRE(rv_self.norm() == Approx(0.0).margin(1e-6));
}

TEST_CASE("interfaces.spice_more.body_pos_frame_and_aberration") {
  Real et = Epoch();

  // Same Earth->Moon vector in J2000 and in the Earth body-fixed frame
  // (IAU_EARTH): the magnitude is rotation-invariant.
  Vec3d r_j2000 = spice::GetBodyPosSpice(et, BodyId::EARTH, BodyId::MOON, "J2000", "NONE");
  Vec3d r_iau = spice::GetBodyPosSpice(et, BodyId::EARTH, BodyId::MOON, "IAU_EARTH", "NONE");
  REQUIRE(r_iau.norm() == Approx(r_j2000.norm()).epsilon(1e-9));
  REQUIRE(r_j2000.norm() > kMoonDistMinKm);
  REQUIRE(r_j2000.norm() < kMoonDistMaxKm);

  // Light-time ("LT") aberration correction shifts the apparent Moon position
  // by roughly a light-time of lunar motion (~1 km) relative to "NONE".
  Vec3d r_lt = spice::GetBodyPosSpice(et, BodyId::EARTH, BodyId::MOON, "J2000", "LT");
  double shift = (r_lt - r_j2000).norm();
  REQUIRE(shift > 0.0);
  REQUIRE(shift < 100.0);  // km
}

TEST_CASE("interfaces.spice_more.frame_conversion_mat_moon") {
  Real et = Epoch();

  // J2000 <-> IAU_MOON 6x6 state rotation: top-left 3x3 is a proper rotation.
  Mat6d M = spice::GetFrameConversionMat(et, "J2000", "IAU_MOON");
  REQUIRE(M.rows() == 6);
  REQUIRE(M.cols() == 6);
  Mat3d R = M.topLeftCorner(3, 3);
  REQUIRE((R * R.transpose()).isApprox(Mat3d::Identity(), 1e-9));
  REQUIRE(R.determinant() == Approx(1.0).margin(1e-9));

  // Forward then inverse recovers identity (6x6, incl. angular-velocity block).
  Mat6d Minv = spice::GetFrameConversionMat(et, "IAU_MOON", "J2000");
  REQUIRE((M * Minv).isApprox(Mat6d::Identity(), 1e-9));
}

TEST_CASE("interfaces.spice_more.planet_orientation_bodies") {
  Real et = Epoch();

  // (alpha0, delta0, W) for Mars / Jupiter: finite, with the pole declination
  // in the valid [-pi/2, pi/2] range and a well-defined prime-meridian angle.
  for (BodyId id : {BodyId::MARS, BodyId::JUPITER}) {
    Vec3d o = spice::GetPlanetOrientation(id, et);
    REQUIRE(std::isfinite(o(0)));
    REQUIRE(std::isfinite(o(1)));
    REQUIRE(std::isfinite(o(2)));
    REQUIRE(o(1) >= -PI / 2.0 - 1e-6);
    REQUIRE(o(1) <= PI / 2.0 + 1e-6);
  }

  // Mars' rotation pole declination is well away from the celestial pole
  // (~52.9 deg), i.e. clearly distinct from Earth's ~90 deg -- a coarse,
  // convention-independent sanity check that the body is resolved distinctly.
  Vec3d mars = spice::GetPlanetOrientation(BodyId::MARS, et);
  REQUIRE(std::abs(mars(1)) < PI / 2.0 - 0.2);
}

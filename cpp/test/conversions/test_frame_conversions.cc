#include <cspice/SpiceUsr.h>
#include <lupnt/conversions/frame_conversions.h>
#include <lupnt/conversions/time_conversions.h>
#include <lupnt/interfaces/spice.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <string>
#include <vector>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("conversions.frame_conversions") {
  SECTION("GCRF and EME rotations round trip position and velocity") {
    Vec6 rv_gcrf(7000.0e3, -1200.0e3, 3000.0e3, 1200.0, 7300.0, -900.0);

    Vec6 rv_eme = GcrfToEme(rv_gcrf);
    Vec6 recovered = EmeToGcrf(rv_eme);

    for (int i = 0; i < 6; ++i)
      REQUIRE_THAT(recovered(i).val(), WithinAbs(rv_gcrf(i).val(), 1.0e-8));
  }

  SECTION("frame-bias rotation is orthonormal") {
    Mat3d R = RotGcrfToEme();
    Mat3d identity = R * R.transpose();

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        REQUIRE_THAT(identity(i, j), WithinAbs(i == j ? 1.0 : 0.0, 1.0e-12));
  }
}

// ---------------------------------------------------------------------
// SPICE-fitted high-accuracy Earth/lunar orientation calibration
// ---------------------------------------------------------------------
//
// Validates InitFrameConversionFromSpice (see frame_conversions.h/.cc):
// after fitting, LuPNT's native (non-SPICE) GcrfToItrf / MoonCiToPa should
// closely reproduce the rotation that SPICE's high-accuracy orientation
// kernels produce -- the same `sxform_c("J2000", "ITRF93"/"IAU_MOON", et,
// xform)` call that spice::ConvertFrameSpice / spice::GetFrameConversionMat
// use internally to do GCRF<->ITRF and MOON_CI<->MOON_PA conversions.
//
// NOTE: this test deliberately does *not* invoke spice::ConvertFrameSpice /
// spice::GetFrameConversionMat with frame_in != frame_out. This test samples
// the same SPICE frame pairs directly with TDB epochs so the fitted native path
// and the SPICE reference use the same epoch scale.
namespace {
  Mat6d SampleSpiceXform(double t_tdb, const std::string& from, const std::string& to) {
    double xform[6][6];
    sxform_c(from.c_str(), to.c_str(), t_tdb, xform);
    Mat6d M;
    for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++) M(i, j) = xform[i][j];
    return M;
  }
}  // namespace

TEST_CASE("conversions.frame_conversions.spice_fit") {
  spice::LoadSpiceKernel();
  ClearFrameConversionFit();

  const double t_start = 0.0;          // seconds since J2000, TDB
  const double t_end = 3.0 * 86400.0;  // 3-day calibration window
  const double t_mid = 0.5 * (t_start + t_end);
  const double t_far_before = t_start - 30.0 * 86400.0;
  const double t_far_after = t_end + 30.0 * 86400.0;
  const double rot_tol = 1.0e-6;  // rotation-matrix-element tolerance

  // Unit-magnitude "rotation-only" probe vectors: applying a 3x3 rotation
  // to a unit vector means the absolute difference between the native and
  // SPICE-sampled results directly reflects the rotation-matrix-element
  // error (no magnitude amplification from large position/velocity values).
  Vec3 r_gcrf(1.0, 0.0, 0.0);
  Vec3 r_mi(0.0, 1.0, 0.0);

  std::vector<double> probe_epochs = {t_start, t_start + 0.37 * (t_end - t_start), t_mid, t_end};

  SECTION("fitted native conversions closely match SPICE inside the fitted window") {
    InitFrameConversionFromSpice(Real(t_start), Real(t_end));

    for (double t : probe_epochs) {
      INFO("epoch t_tdb = " << t);

      // Earth: GCRF -> ITRF
      REQUIRE(HasFittedEarthOrientation(Real(t)));
      Vec3 r_itrf_native = GcrfToItrf(Real(t), r_gcrf);
      Mat6d xform_e = SampleSpiceXform(t, "J2000", "ITRF93");
      Vec3d r_itrf_spice = xform_e.topLeftCorner<3, 3>() * r_gcrf.cast<double>();
      for (int i = 0; i < 3; i++)
        REQUIRE_THAT(r_itrf_native(i).val(), WithinAbs(r_itrf_spice(i), rot_tol));

      // Moon: MOON_CI -> MOON_PA
      REQUIRE(HasFittedLunarOrientation(Real(t)));
      Vec3 r_pa_native = MoonCiToPa(Real(t), r_mi);
      Mat6d xform_m = SampleSpiceXform(t, "J2000", "IAU_MOON");
      Vec3d r_pa_spice = xform_m.topLeftCorner<3, 3>() * r_mi.cast<double>();
      for (int i = 0; i < 3; i++)
        REQUIRE_THAT(r_pa_native(i).val(), WithinAbs(r_pa_spice(i), rot_tol));

      // Round-trip: ItrfToGcrf(GcrfToItrf(r)) must recover r exactly.
      // With angle-based fitting, R_fit is EXACTLY orthonormal (guaranteed
      // by the ZXZ composition), so R_fit^T * R_fit = I to machine precision
      // and the round-trip error is limited by floating-point rounding only.
      Vec3 r_gcrf_rt = ItrfToGcrf(Real(t), r_itrf_native);
      Vec3 r_mi_rt = MoonPaToCi(Real(t), r_pa_native);
      for (int i = 0; i < 3; i++) {
        REQUIRE_THAT(r_gcrf_rt(i).val(), WithinAbs(r_gcrf(i).val(), 1.0e-12));
        REQUIRE_THAT(r_mi_rt(i).val(), WithinAbs(r_mi(i).val(), 1.0e-12));
      }

      // Explicit orthonormality: reconstruct R by applying it to the 3
      // standard basis vectors, then verify R^T*R = I to near machine
      // precision.  This checks the key property that angle-based fitting
      // guarantees -- it would have FAILED at ~4e-8 under the old
      // matrix-element-fitting approach.
      {
        const double orth_tol = 1.0e-14;
        Vec3 e0(1.0, 0.0, 0.0), e1(0.0, 1.0, 0.0), e2(0.0, 0.0, 1.0);
        // Earth columns
        Vec3 Re0 = GcrfToItrf(Real(t), e0);
        Vec3 Re1 = GcrfToItrf(Real(t), e1);
        Vec3 Re2 = GcrfToItrf(Real(t), e2);
        REQUIRE_THAT(Re0.squaredNorm().val(), WithinAbs(1.0, orth_tol));
        REQUIRE_THAT(Re1.squaredNorm().val(), WithinAbs(1.0, orth_tol));
        REQUIRE_THAT(Re2.squaredNorm().val(), WithinAbs(1.0, orth_tol));
        REQUIRE_THAT(Re0.dot(Re1).val(), WithinAbs(0.0, orth_tol));
        REQUIRE_THAT(Re0.dot(Re2).val(), WithinAbs(0.0, orth_tol));
        REQUIRE_THAT(Re1.dot(Re2).val(), WithinAbs(0.0, orth_tol));
        // Moon columns
        Vec3 Rm0 = MoonCiToPa(Real(t), e0);
        Vec3 Rm1 = MoonCiToPa(Real(t), e1);
        Vec3 Rm2 = MoonCiToPa(Real(t), e2);
        REQUIRE_THAT(Rm0.squaredNorm().val(), WithinAbs(1.0, orth_tol));
        REQUIRE_THAT(Rm1.squaredNorm().val(), WithinAbs(1.0, orth_tol));
        REQUIRE_THAT(Rm2.squaredNorm().val(), WithinAbs(1.0, orth_tol));
        REQUIRE_THAT(Rm0.dot(Rm1).val(), WithinAbs(0.0, orth_tol));
        REQUIRE_THAT(Rm0.dot(Rm2).val(), WithinAbs(0.0, orth_tol));
        REQUIRE_THAT(Rm1.dot(Rm2).val(), WithinAbs(0.0, orth_tol));
      }
    }
  }

  SECTION("HasFittedEarthOrientation / HasFittedLunarOrientation report window coverage") {
    // No fit installed yet -- nothing should report as covered.
    REQUIRE_FALSE(HasFittedEarthOrientation(Real(t_mid)));
    REQUIRE_FALSE(HasFittedLunarOrientation(Real(t_mid)));

    InitFrameConversionFromSpice(Real(t_start), Real(t_end));

    REQUIRE(HasFittedEarthOrientation(Real(t_start)));
    REQUIRE(HasFittedEarthOrientation(Real(t_mid)));
    REQUIRE(HasFittedEarthOrientation(Real(t_end)));
    REQUIRE(HasFittedLunarOrientation(Real(t_start)));
    REQUIRE(HasFittedLunarOrientation(Real(t_mid)));
    REQUIRE(HasFittedLunarOrientation(Real(t_end)));

    REQUIRE_FALSE(HasFittedEarthOrientation(Real(t_far_before)));
    REQUIRE_FALSE(HasFittedEarthOrientation(Real(t_far_after)));
    REQUIRE_FALSE(HasFittedLunarOrientation(Real(t_far_before)));
    REQUIRE_FALSE(HasFittedLunarOrientation(Real(t_far_after)));
  }

  SECTION("ComputeEopFromSpice round-trips the SPICE J2000->ITRF93 rotation") {
    for (double t : probe_epochs) {
      INFO("epoch t_tdb = " << t);
      Real t_tdb(t);

      SpiceEopParams p = ComputeEopFromSpice(t_tdb);

      Real t_tt = ConvertTime(t_tdb, Time::TDB, Time::TT);
      Real t_utc = ConvertTime(t_tdb, Time::TDB, Time::UTC);
      Real t_ut1 = t_utc + p.ut1_utc;
      Real theta_era = EarthRotationAngle(t_ut1);
      Real sp = -47e-6 * RAD_ARCSEC * (t_tt / DAYS_CENTURY);

      Mat3 R_po = RotX(-p.y_pole) * RotY(-p.x_pole) * RotZ(sp);
      Mat3 R_s = RotZ(theta_era);
      Mat3 R_pn = RotPrecessionNutation(t_tdb);
      Mat3 R_reconstructed = R_po * R_s * R_pn;

      Mat6d xform = SampleSpiceXform(t, "J2000", "ITRF93");
      Mat3d R_spice = xform.topLeftCorner<3, 3>();

      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          REQUIRE_THAT(R_reconstructed(i, j).val(), WithinAbs(R_spice(i, j), 1.0e-9));
    }
  }

  SECTION("GcrfToItrf with SPICE-fitted EOP vs. the bundled IERS EOP table agree closely") {
    Vec3 r_itrf_table = GcrfToItrf(Real(t_mid), r_gcrf);
    REQUIRE_FALSE(HasFittedEarthOrientation(Real(t_mid)));

    InitFrameConversionFromSpice(Real(t_start), Real(t_end));
    REQUIRE(HasFittedEarthOrientation(Real(t_mid)));
    Vec3 r_itrf_fit = GcrfToItrf(Real(t_mid), r_gcrf);

    for (int i = 0; i < 3; i++)
      REQUIRE_THAT(r_itrf_fit(i).val(), WithinAbs(r_itrf_table(i).val(), 1.0e-5));
  }

  SECTION("ClearFrameConversionFit reverts native conversions to the original analytic results") {
    // Capture the original (pre-fit) analytic results.
    Vec3 r_itrf_before = GcrfToItrf(Real(t_mid), r_gcrf);
    Vec3 r_pa_before = MoonCiToPa(Real(t_mid), r_mi);
    REQUIRE_FALSE(HasFittedEarthOrientation(Real(t_mid)));
    REQUIRE_FALSE(HasFittedLunarOrientation(Real(t_mid)));

    InitFrameConversionFromSpice(Real(t_start), Real(t_end));
    REQUIRE(HasFittedEarthOrientation(Real(t_mid)));
    REQUIRE(HasFittedLunarOrientation(Real(t_mid)));

    ClearFrameConversionFit();
    REQUIRE_FALSE(HasFittedEarthOrientation(Real(t_mid)));
    REQUIRE_FALSE(HasFittedLunarOrientation(Real(t_mid)));

    Vec3 r_itrf_after = GcrfToItrf(Real(t_mid), r_gcrf);
    Vec3 r_pa_after = MoonCiToPa(Real(t_mid), r_mi);
    for (int i = 0; i < 3; i++) {
      REQUIRE_THAT(r_itrf_after(i).val(), WithinAbs(r_itrf_before(i).val(), 1.0e-12));
      REQUIRE_THAT(r_pa_after(i).val(), WithinAbs(r_pa_before(i).val(), 1.0e-12));
    }
  }

  ClearFrameConversionFit();
}

TEST_CASE("conversions.frame_conversions.planet_frame_predicates") {
  // The predicates recognise the generic-IAU planet frames (Mercury..Neptune),
  // and deliberately exclude Earth (GCRF/ITRF) and the Moon, which have their own
  // dedicated frame machinery.
  SECTION("planet body-fixed frames") {
    REQUIRE(IsPlanetFixedFrame(Frame::MARS_FIXED));
    REQUIRE(IsPlanetFixedFrame(Frame::JUPITER_FIXED));
    REQUIRE(IsPlanetFixedFrame(Frame::NEPTUNE_FIXED));
    REQUIRE_FALSE(IsPlanetFixedFrame(Frame::MARS_CI));
    REQUIRE_FALSE(IsPlanetFixedFrame(Frame::ITRF));
    REQUIRE_FALSE(IsPlanetFixedFrame(Frame::GCRF));
  }
  SECTION("planet inertial frames") {
    REQUIRE(IsPlanetCiFrame(Frame::MARS_CI));
    REQUIRE(IsPlanetCiFrame(Frame::VENUS_CI));
    REQUIRE_FALSE(IsPlanetCiFrame(Frame::MARS_FIXED));
    REQUIRE_FALSE(IsPlanetCiFrame(Frame::GCRF));
  }
  SECTION("IsPlanetFrame is the union of fixed and inertial") {
    REQUIRE(IsPlanetFrame(Frame::MARS_FIXED));
    REQUIRE(IsPlanetFrame(Frame::MARS_CI));
    REQUIRE_FALSE(IsPlanetFrame(Frame::GCRF));
    REQUIRE_FALSE(IsPlanetFrame(Frame::MOON_CI));
  }
}

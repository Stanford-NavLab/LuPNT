#include <lupnt/conversions/frame_converter.h>
#include <lupnt/conversions/time_conversions.h>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

// ============================================================================
// Cross-validation against GMAT
//
// The reference values below were generated once via GMAT (R2022a) and
// checked into `cpp/test/gmat/data/gmat_reference.json` by
// `cpp/test/gmat/gen_gmat_reference.py`. This test does **not** require
// GMAT to run -- see `cpp/test/gmat/README.md` for how the fixture was
// produced and how to regenerate it.
//
// `r_gcrf`/`v_gcrf` below are the same state expressed in GMAT's "ICRF"
// coordinate system (LuPNT/Orekit's GCRF), `r_eme2000`/`v_eme2000` are GMAT's
// "EarthMJ2000Eq" axes (LuPNT/Orekit's EME2000/J2000), and `r_itrf`/`v_itrf`
// are GMAT's "BodyFixed" axes (LuPNT/Orekit's ITRF).
// ============================================================================

namespace {
  Vec6 ReadVec6(const nlohmann::json& r, const nlohmann::json& v) {
    Vec6 rv;
    for (int i = 0; i < 3; i++) {
      rv(i) = r.at(i).get<double>();
      rv(i + 3) = v.at(i).get<double>();
    }
    return rv;
  }
}  // namespace

TEST_CASE("conversions.frame_conversions_gmat_reference") {
  nlohmann::json data = LoadTestJson("gmat/data/gmat_reference.json");

  for (const auto& tc : data["frames"]) {
    std::string epoch = tc["epoch_utc"].get<std::string>();
    Real t_tdb = ConvertTime(GregorianToTime(epoch), Time::UTC, Time::TDB);

    Vec6 rv_gcrf = ReadVec6(tc["r_gcrf"], tc["v_gcrf"]);
    Vec6 rv_eme2000_ref = ReadVec6(tc["r_eme2000"], tc["v_eme2000"]);
    Vec6 rv_itrf_ref = ReadVec6(tc["r_itrf"], tc["v_itrf"]);

    DYNAMIC_SECTION("epoch " << epoch << ", r_gcrf = " << rv_gcrf.head<3>().transpose()) {
      // GCRF <-> EME2000 (J2000) is a fixed frame-bias rotation that does
      // not depend on Earth Orientation Parameters. However, GMAT's "ICRF"
      // -> "EarthMJ2000Eq" rotation differs from LuPNT's GCRF -> EME2000
      // frame-bias matrix (which matches Orekit's to sub-mm, see
      // test_frame_conversions_orekit.cc) by a small angle. Across LEO to
      // trans-lunar radii and 2010-2025 epochs, the observed difference
      // scales with |r| as ~0.8-1.6e-7 * |r| (i.e. a ~1e-7 rad rotation
      // difference, slightly epoch-dependent), so the position tolerance
      // scales as 3e-7 * |r| (with a 5 m floor), ~2x the observed maximum
      // and still ~7 orders of magnitude tighter than a gross error (wrong
      // rotation, wrong sign), which would produce an error comparable to
      // |r| itself.
      Vec6 rv_eme2000 = ConvertFrame(t_tdb, rv_gcrf, Frame::GCRF, Frame::EME);
      double r_tol_eme = std::max(5.0, 3e-7 * rv_gcrf.head<3>().norm().val());
      for (int i = 0; i < 3; i++) {
        RequireNear(rv_eme2000(i), rv_eme2000_ref(i), r_tol_eme);     // [m]
        RequireNear(rv_eme2000(i + 3), rv_eme2000_ref(i + 3), 2e-3);  // [m/s]
      }

      // Round trip back to GCRF should recover the input exactly (to within
      // numerical precision).
      Vec6 rv_gcrf_rt = ConvertFrame(t_tdb, rv_eme2000, Frame::EME, Frame::GCRF);
      for (int i = 0; i < 6; i++) RequireNear(rv_gcrf_rt(i), rv_gcrf(i), 1e-6);

      // GCRF <-> ITRF is the full Earth-orientation (precession, nutation,
      // sidereal rotation, polar motion) transform, so it depends on the
      // Earth Orientation Parameters (EOP) table in use. LuPNT's bundled EOP
      // table and GMAT's bundled EOP/SPICE kernels correspond to different
      // IERS bulletin "vintages" and can disagree in UT1-UTC (and polar
      // motion) by ~tens of ms. A small error `dtheta` in the
      // Earth-rotation angle rotates the position vector by `dtheta` about
      // the polar axis, so the resulting Cartesian error scales with the
      // *radius*, not a fixed distance: `|delta r| ~ dtheta * |r|`. The
      // tolerance below uses the same formula (and the same ~tens-of-ms
      // order of magnitude) as the Orekit ITRF comparison in
      // test_frame_conversions_orekit.cc.
      Vec6 rv_itrf = ConvertFrame(t_tdb, rv_gcrf, Frame::GCRF, Frame::ITRF);
      double r_tol_itrf = std::max(100.0, 2e-5 * rv_gcrf.head<3>().norm().val());
      for (int i = 0; i < 3; i++) {
        RequireNear(rv_itrf(i), rv_itrf_ref(i), r_tol_itrf);   // [m]
        RequireNear(rv_itrf(i + 3), rv_itrf_ref(i + 3), 1.0);  // [m/s]
      }
    }
  }
}

// ----------------------------------------------------------------------------
// Moon-centered frames
//
// GMAT R2022a ships DE405/DE421/DE424 (no DE440); the fixture was generated
// with DE421, the closest to LuPNT's DE440. The dominant difference in both
// comparisons below is therefore the DE421 vs DE440 lunar-ephemeris offset
// (observed ~55-85 m in the 2024-2025 epochs), not a frame-modeling error.
// ----------------------------------------------------------------------------
TEST_CASE("conversions.moon_frames_gmat_reference") {
  nlohmann::json data = LoadTestJson("gmat/data/gmat_reference.json");

  for (const auto& tc : data["moon_frames"]) {
    std::string epoch = tc["epoch_utc"].get<std::string>();
    Real t_tdb = ConvertTime(GregorianToTime(epoch), Time::UTC, Time::TDB);

    Vec6 rv_gcrf = ReadVec6(tc["r_gcrf"], tc["v_gcrf"]);
    Vec6 rv_moon_ci_ref = ReadVec6(tc["r_moon_ci"], tc["v_moon_ci"]);
    Vec6 rv_moon_fixed_ref = ReadVec6(tc["r_moon_fixed"], tc["v_moon_fixed"]);

    DYNAMIC_SECTION("epoch " << epoch << ", r_moon_ci = " << rv_moon_ci_ref.head<3>().transpose()) {
      // GCRF -> MOON_CI is a pure translation by the Earth->Moon ephemeris
      // vector: LuPNT evaluates DE440, GMAT evaluated DE421, and the two
      // ephemerides differ by ~55-85 m (position) / ~2e-4 m/s (velocity)
      // at these epochs. 200 m / 1e-3 m/s gives ~2.5-5x margin while still
      // catching any real error (wrong center, wrong frame, km-vs-m),
      // which would show up at 1e3-1e8 m scale.
      Vec6 rv_moon_ci = ConvertFrame(t_tdb, rv_gcrf, Frame::GCRF, Frame::MOON_CI);
      for (int i = 0; i < 3; i++) {
        RequireNear(rv_moon_ci(i), rv_moon_ci_ref(i), 200.0);         // [m]
        RequireNear(rv_moon_ci(i + 3), rv_moon_ci_ref(i + 3), 1e-3);  // [m/s]
      }

      // Round trip back to GCRF should recover the input exactly.
      Vec6 rv_gcrf_rt = ConvertFrame(t_tdb, rv_moon_ci, Frame::MOON_CI, Frame::GCRF);
      for (int i = 0; i < 6; i++) RequireNear(rv_gcrf_rt(i), rv_gcrf(i), 1e-6);

      // GMAT's Luna "BodyFixed" axes use a DE-era lunar *principal axes*
      // frame (GMAT loads a SPICE Luna frame kernel), so it is compared
      // against LuPNT's MOON_PA (DE440 principal axes). The observed
      // difference (~55-90 m, ~1e-3 m/s) is dominated by the same DE421 vs
      // DE440 translation as the MOON_CI comparison above, i.e. the PA
      // *orientations* agree at the sub-arcsecond level. Note: comparing
      // against LuPNT's MOON_ME instead fails by ~10-40x more
      // (~0.6-3.4 km), confirming GMAT's Luna BodyFixed is the
      // principal-axes frame, complementary to the Orekit comparison
      // (which validates MOON_ME against the IAU model).
      Vec6 rv_moon_pa = ConvertFrame(t_tdb, rv_gcrf, Frame::GCRF, Frame::MOON_PA);
      for (int i = 0; i < 3; i++) {
        RequireNear(rv_moon_pa(i), rv_moon_fixed_ref(i), 300.0);         // [m]
        RequireNear(rv_moon_pa(i + 3), rv_moon_fixed_ref(i + 3), 5e-3);  // [m/s]
      }
    }
  }
}

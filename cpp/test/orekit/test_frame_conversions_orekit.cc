#include <lupnt/conversions/frame_converter.h>
#include <lupnt/conversions/time_conversions.h>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

// ============================================================================
// Cross-validation against Orekit
//
// The reference values below were generated once via Orekit and checked
// into `cpp/test/orekit/data/orekit_reference.json` by
// `cpp/test/orekit/gen_orekit_reference.py`. This test does **not**
// require Orekit/Java to run -- see `cpp/test/orekit/README.md` for how
// the fixture was produced and how to regenerate it.
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

TEST_CASE("conversions.frame_conversions_orekit_reference") {
  nlohmann::json data = LoadTestJson("orekit/data/orekit_reference.json");

  for (const auto& tc : data["frames"]) {
    std::string epoch = tc["epoch_utc"].get<std::string>();
    Real t_tdb = ConvertTime(GregorianToTime(epoch), Time::UTC, Time::TDB);

    Vec6 rv_gcrf = ReadVec6(tc["r_gcrf"], tc["v_gcrf"]);
    Vec6 rv_eme2000_ref = ReadVec6(tc["r_eme2000"], tc["v_eme2000"]);
    Vec6 rv_itrf_ref = ReadVec6(tc["r_itrf"], tc["v_itrf"]);

    DYNAMIC_SECTION("epoch " << epoch << ", r_gcrf = " << rv_gcrf.head<3>().transpose()) {
      // GCRF <-> EME2000 (J2000) is a fixed frame-bias rotation that does
      // *not* depend on Earth Orientation Parameters, so LuPNT and Orekit
      // should agree to well below a millimeter.
      Vec6 rv_eme2000 = ConvertFrame(t_tdb, rv_gcrf, Frame::GCRF, Frame::EME);
      for (int i = 0; i < 3; i++) {
        RequireNear(rv_eme2000(i), rv_eme2000_ref(i), 1e-3);          // [m]
        RequireNear(rv_eme2000(i + 3), rv_eme2000_ref(i + 3), 1e-6);  // [m/s]
      }

      // Round trip back to GCRF should recover the input exactly (to within
      // numerical precision).
      Vec6 rv_gcrf_rt = ConvertFrame(t_tdb, rv_eme2000, Frame::EME, Frame::GCRF);
      for (int i = 0; i < 6; i++) RequireNear(rv_gcrf_rt(i), rv_gcrf(i), 1e-6);

      // GCRF <-> ITRF is the full Earth-orientation (precession, nutation,
      // sidereal rotation, polar motion) transform, so it depends on the
      // Earth Orientation Parameters (EOP) table in use. LuPNT's bundled EOP
      // table and Orekit's `orekit-data` table correspond to different IERS
      // bulletin "vintages" and can disagree in UT1-UTC (and polar motion)
      // by ~tens of ms. A small error `dtheta` in the Earth-rotation angle
      // rotates the position vector by `dtheta` about the polar axis, so the
      // resulting Cartesian error scales with the *radius*, not a fixed
      // distance: `|delta r| ~ dtheta * |r|`. For a LEO case
      // (|r| ~ 7.3e6 m) this is only ~11 m, but for a GEO-altitude case
      // (|r| ~ 3.1e7 m) the same `dtheta` gives ~400 m. The tolerance below
      // is therefore scaled with |r_gcrf| (with a 100 m floor for small
      // |r|), calibrated so that a UT1-UTC/polar-motion difference of order
      // tens of ms still passes for any orbit radius, while a gross error
      // (wrong rotation, wrong sign, missing polar motion, etc. -- which
      // would produce an error comparable to |r| itself) is still caught.
      Vec6 rv_itrf = ConvertFrame(t_tdb, rv_gcrf, Frame::GCRF, Frame::ITRF);
      double r_tol = std::max(100.0, 2e-5 * rv_gcrf.head<3>().norm().val());
      for (int i = 0; i < 3; i++) {
        RequireNear(rv_itrf(i), rv_itrf_ref(i), r_tol);        // [m]
        RequireNear(rv_itrf(i + 3), rv_itrf_ref(i + 3), 1.0);  // [m/s]
      }
    }
  }
}

// ----------------------------------------------------------------------------
// Moon-centered frames
// ----------------------------------------------------------------------------
TEST_CASE("conversions.moon_frames_orekit_reference") {
  nlohmann::json data = LoadTestJson("orekit/data/orekit_reference.json");

  for (const auto& tc : data["moon_frames"]) {
    std::string epoch = tc["epoch_utc"].get<std::string>();
    Real t_tdb = ConvertTime(GregorianToTime(epoch), Time::UTC, Time::TDB);

    Vec6 rv_gcrf = ReadVec6(tc["r_gcrf"], tc["v_gcrf"]);
    Vec6 rv_moon_ci_ref = ReadVec6(tc["r_moon_ci"], tc["v_moon_ci"]);
    Vec6 rv_moon_fixed_ref = ReadVec6(tc["r_moon_fixed_iau"], tc["v_moon_fixed_iau"]);

    DYNAMIC_SECTION("epoch " << epoch << ", r_moon_ci = " << rv_moon_ci_ref.head<3>().transpose()) {
      // GCRF -> MOON_CI is a pure translation by the Earth->Moon ephemeris
      // vector (both LuPNT and Orekit evaluate DE440-series ephemerides),
      // so this is effectively a cross-check of the lunar ephemeris
      // evaluation through completely independent code paths (LuPNT's
      // Chebyshev-coefficient reader vs Orekit's JPL-DE loader, on
      // separately distributed copies of the DE440 data). Observed
      // agreement is ~0.2-0.3 m in position (sub-nm/s-level relative) --
      // consistent with small differences between DE440 distributions --
      // so 1 m / 1e-5 m/s catches any real indexing/units/time-argument
      // error (which would show up at km scale) with ~3x margin.
      Vec6 rv_moon_ci = ConvertFrame(t_tdb, rv_gcrf, Frame::GCRF, Frame::MOON_CI);
      for (int i = 0; i < 3; i++) {
        RequireNear(rv_moon_ci(i), rv_moon_ci_ref(i), 1.0);           // [m]
        RequireNear(rv_moon_ci(i + 3), rv_moon_ci_ref(i + 3), 1e-5);  // [m/s]
      }

      // Round trip back to GCRF should recover the input exactly.
      Vec6 rv_gcrf_rt = ConvertFrame(t_tdb, rv_moon_ci, Frame::MOON_CI, Frame::GCRF);
      for (int i = 0; i < 6; i++) RequireNear(rv_gcrf_rt(i), rv_gcrf(i), 1e-6);

      // GCRF -> MOON_ME (Moon body-fixed, mean-Earth/polar axes). The
      // reference was computed with Orekit's body-oriented Moon frame, which
      // uses the analytical IAU pole/prime-meridian model, whereas LuPNT
      // derives MOON_PA from the DE440 principal-axes kernel and applies the
      // fixed PA->ME rotation. The IAU model approximates the DE mean-Earth
      // orientation to ~1e-4 rad (~150 m at the lunar surface, per the
      // NAIF/IAU WG reports); the Cartesian error of a small rotation
      // difference scales with the lunar-centric radius
      // (|delta r| ~ dtheta * |r_moon_ci|). Observed differences here are
      // ~3e-5 * |r| (57-200 m at 2000-6900 km radii), so the 3e-4 * |r|
      // tolerance has ~10x margin while still being ~30x tighter than a
      // wrong-axes/wrong-angle bug (error ~ |r|). Note: comparing against
      // LuPNT's MOON_PA instead fails by ~15x more (~880-3100 m),
      // confirming the IAU model tracks the *mean-Earth* axes, not the
      // principal axes. The velocity difference is dominated by the same
      // angular offset applied to the ~km/s moon-relative velocities
      // (observed up to ~0.024 m/s).
      Vec6 rv_moon_me = ConvertFrame(t_tdb, rv_gcrf, Frame::GCRF, Frame::MOON_ME);
      double r_tol_me = 3e-4 * rv_moon_ci_ref.head<3>().norm().val();
      for (int i = 0; i < 3; i++) {
        RequireNear(rv_moon_me(i), rv_moon_fixed_ref(i), r_tol_me);      // [m]
        RequireNear(rv_moon_me(i + 3), rv_moon_fixed_ref(i + 3), 0.05);  // [m/s]
      }

      // MOON_CI <-> MOON_PA <-> MOON_ME round trip (LuPNT-internal
      // consistency, anchored on the Orekit-verified MOON_CI state).
      Vec6 rv_moon_pa = ConvertFrame(t_tdb, rv_moon_ci, Frame::MOON_CI, Frame::MOON_PA);
      Vec6 rv_moon_ci_rt = ConvertFrame(t_tdb, rv_moon_pa, Frame::MOON_PA, Frame::MOON_CI);
      for (int i = 0; i < 6; i++) RequireNear(rv_moon_ci_rt(i), rv_moon_ci(i), 1e-6);
    }
  }
}

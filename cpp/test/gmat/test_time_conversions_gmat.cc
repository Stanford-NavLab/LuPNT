#include <lupnt/conversions/time_conversions.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

// ============================================================================
// Cross-validation against GMAT
//
// The reference values below were generated once via GMAT (R2026a) and
// checked into `cpp/test/gmat/data/gmat_reference.json` by
// `cpp/test/gmat/gen_gmat_reference.py`. This test does **not** require
// GMAT to run -- see `cpp/test/gmat/README.md` for how the fixture was
// produced and how to regenerate it.
//
// Unlike the Orekit reference, GMAT does not expose a UT1 or GPS time scale
// as a spacecraft parameter, so only TAI-UTC, TT-TAI, and TDB-TAI are
// cross-checked here.
// ============================================================================

TEST_CASE("conversions.time_conversions_gmat_reference") {
  // GMAT and Orekit both use an analytic (Fairhead-Bretagnon style) TT<->TDB
  // series. LuPNT's default is now the integrated DE440t model, which differs
  // from those series by ~20-30 us -- LuPNT is the more accurate of the two, so
  // this is not a LuPNT error. To keep this a like-for-like check of the
  // *algorithm*, compare against LuPNT's analytic series here.
  const bool autofit_was = GetTtTdbAutoFit();
  SetTtTdbAutoFit(false);
  ClearTtMinusTdbFit();

  nlohmann::json data = LoadTestJson("gmat/data/gmat_reference.json");

  for (const auto& tc : data["time_scales"]) {
    std::string epoch = tc["epoch_utc"].get<std::string>();
    Real t_utc = GregorianToTime(epoch);

    DYNAMIC_SECTION("epoch " << epoch) {
      Real t_tai = ConvertTime(t_utc, Time::UTC, Time::TAI);
      Real t_tt = ConvertTime(t_utc, Time::UTC, Time::TT);
      Real t_tdb = ConvertTime(t_utc, Time::UTC, Time::TDB);

      // TAI-UTC (leap seconds) and TT-TAI (exact 32.184 s) are constant
      // offsets defined by the time-scale standards themselves, so
      // algorithmically LuPNT and GMAT should agree exactly. GMAT reports
      // these as `<X>ModJulian` fields (days), so `gen_gmat_reference.py`
      // recovers the offset as `(X_mjd - TAI_mjd) * 86400`. Subtracting two
      // ~30000-day Modified Julian Dates and multiplying by 86400 amplifies
      // the ULP of the underlying double (~1e-7 s), so GMAT's reported
      // offsets differ from the exact constants by up to ~2e-7 s. This is a
      // floating-point representation artifact (on GMAT's side, mirroring
      // the same kind of artifact documented in
      // test_time_conversions_orekit.cc), not an algorithmic discrepancy --
      // 1e-6 s (1 microsecond) comfortably absorbs it while still catching a
      // real (e.g. wrong leap-second count, wrong sign) error.
      RequireNear(t_tai - t_utc, Real(tc["tai_minus_utc"].get<double>()), 1e-6);
      RequireNear(t_tt - t_tai, Real(tc["tt_minus_tai"].get<double>()), 1e-6);

      // TDB-TAI has a small (~ms-amplitude) periodic term that depends on
      // the analytical expression used for the relativistic clock
      // correction. LuPNT's and GMAT's implementations use slightly
      // different truncated series for this term, so they can disagree by a
      // few microseconds (in addition to the ~1e-7 s representation
      // artifact above) -- still many orders of magnitude smaller than the
      // ~1.7 ms amplitude of the periodic term itself, so 1e-5 s (10
      // microseconds) is generous while still catching a gross error (wrong
      // sign, missing term, etc.).
      RequireNear(t_tdb - t_tai, Real(tc["tdb_minus_tai"].get<double>()), 1e-5);
    }
  }

  SetTtTdbAutoFit(autofit_was);
}

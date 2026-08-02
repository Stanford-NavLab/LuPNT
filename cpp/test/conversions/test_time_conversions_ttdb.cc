// Evaluation of the high-fidelity TT<->TDB path.
//
//   1. spice::ConvertTime now reads TT-TDB from the DE440t time ephemeris
//      (NAIF body 1000000001), not CSPICE unitim_c's truncated analytic
//      series, so it differs from the native analytic TDBToTt by tens of us.
//   2. InitTtMinusTdbFit() builds a piecewise-Chebyshev fit of that same
//      DE440t TT-TDB, so native TDBToTt()/TtToTdb() reproduce the SPICE result
//      with NO SPICE call (safe for OpenMP-parallel regions).
//
// Note on precision: LuPNT epochs are doubles holding seconds from J2000, so
// one ULP is |t|*2^-52 (~0.29 us at 2024, ~0.5 ns near J2000). The sub-us
// TT-TDB offset therefore only survives in an *absolute* TT epoch near J2000.
// Tight (< 1 ns) checks here are done close to J2000; away from it, absolute
// TT epochs from the fit and from SPICE still agree to the double-precision
// floor.

#include <lupnt/conversions/time_conversions.h>
#include <lupnt/interfaces/kernels.h>
#include <lupnt/interfaces/spice.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  const double day = SECS_DAY;
  // TT - TDB [s] at TDB epoch t, from SPICE (DE440t), evaluated as a small
  // offset (subtracting the large epoch back off).
  double SpiceTtMinusTdb(Real t) { return (spice::ConvertTime(t, Time::TDB, Time::TT) - t).val(); }
  // TT - TDB [s] from the native TDBToTt (fitted fast path when a fit is active).
  double NativeTtMinusTdb(Real t) { return (TDBToTt(t) - t).val(); }
}  // namespace

TEST_CASE("conversions.tt_tdb_high_fidelity") {
  SECTION("spice::ConvertTime TDB<->TT round-trips") {
    for (double d : {0.0, 5.0, 20.0, -15.0}) {
      Real t_tdb = d * day;  // near J2000 so the round-trip floor is < 1 ns
      Real t_tt = spice::ConvertTime(t_tdb, Time::TDB, Time::TT);
      Real t_tdb_rt = spice::ConvertTime(t_tt, Time::TT, Time::TDB);
      REQUIRE_THAT(t_tdb_rt.val(), WithinAbs(t_tdb.val(), 1.0e-9));
    }
  }

  SECTION("DE440t path differs from the native analytic series by tens of us") {
    // The DE440t integrated TT-TDB and the analytic series are both valid
    // TT-TDB models, so they must agree at the ~ms level yet differ
    // meaningfully (> 1 us) -- that difference is the fidelity gain.
    //
    // The native default is now the DE440t Chebyshev fit (auto-built on demand),
    // so the analytic series must be requested explicitly for this comparison;
    // otherwise both sides are DE440t and agree to ~1e-14 s.
    const bool autofit_was = GetTtTdbAutoFit();
    SetTtTdbAutoFit(false);
    ClearTtMinusTdbFit();
    double max_diff = 0.0;
    for (int i = 0; i <= 40; i++) {
      Real t_tdb = (i * 9.0) * day;  // 0 .. 360 days from J2000
      max_diff = std::max(max_diff, std::abs(SpiceTtMinusTdb(t_tdb) - NativeTtMinusTdb(t_tdb)));
    }
    INFO("max |DE440t - analytic| = " << max_diff << " s");
    SetTtTdbAutoFit(autofit_was);
    REQUIRE(max_diff > 1.0e-6);  // meaningfully different (> 1 us)
    REQUIRE(max_diff < 1.0e-3);  // but both are sane TT-TDB (< 1 ms apart)
  }

  SECTION("native Chebyshev fit reproduces the DE440t SPICE TT-TDB to < 1 ns") {
    // Fit near J2000 so the offset-level comparison is not ULP-limited.
    const Real t_start = -20.0 * day;
    const Real t_end = 20.0 * day;
    InitTtMinusTdbFit(t_start, t_end);  // defaults: 16-day segments, 13 coeffs

    double max_err = 0.0;
    for (int i = 0; i <= 400; i++) {
      Real t_tdb = t_start + (t_end - t_start) * (i / 400.0);
      REQUIRE(HasFittedTtMinusTdb(t_tdb));
      max_err = std::max(max_err, std::abs(NativeTtMinusTdb(t_tdb) - SpiceTtMinusTdb(t_tdb)));
    }
    INFO("max |fitted - DE440t SPICE| = " << max_err << " s");
    REQUIRE(max_err < 1.0e-9);

    ClearTtMinusTdbFit();
    REQUIRE_FALSE(HasFittedTtMinusTdb(0.0));
  }

  SECTION("fitted native TDBToTt matches spice::ConvertTime across a long window") {
    // End-to-end over a 2024-2025 mission window: the two absolute TT epochs
    // agree to the double-precision floor (~1 ULP ~ 0.3 us at these epochs).
    const Real t0 = GregorianToTime(2024, 1, 1, 0, 0, 0.0);
    const Real t_start = t0 - 30.0 * day;
    const Real t_end = t0 + 400.0 * day;
    InitTtMinusTdbFit(t_start, t_end);

    double max_err = 0.0, ulp = 0.0;
    for (int i = 0; i <= 200; i++) {
      Real t_tdb = t_start + (t_end - t_start) * (i / 200.0);
      REQUIRE(HasFittedTtMinusTdb(t_tdb));
      Real fitted = TDBToTt(t_tdb);
      Real spice_ref = spice::ConvertTime(t_tdb, Time::TDB, Time::TT);
      max_err = std::max(max_err, std::abs((fitted - spice_ref).val()));
      ulp = std::max(ulp, std::abs(fitted.val()) * std::pow(2.0, -52));
    }
    INFO("max |fitted - DE440t SPICE| = " << max_err << " s, ulp ~ " << ulp << " s");
    REQUIRE(max_err <= ulp);  // agree to the last representable bit
    ClearTtMinusTdbFit();
  }

  SECTION("native fitted TDB->TT->TDB round-trips to < 1 ns") {
    const Real t_start = -20.0 * day;
    const Real t_end = 20.0 * day;
    InitTtMinusTdbFit(t_start, t_end);

    for (double d : {0.0, 5.0, 12.0, -18.0}) {
      Real t_tdb = d * day;
      Real t_tt = TDBToTt(t_tdb);
      Real t_tdb_rt = TtToTdb(t_tt);
      REQUIRE_THAT(t_tdb_rt.val(), WithinAbs(t_tdb.val(), 1.0e-9));
    }
    ClearTtMinusTdbFit();
  }

  SECTION("SetTtTdbModel(ANALYTIC) returns the series even with auto-fit on") {
    // Regression: the model selector must win over the auto-fit fast path.
    // Previously TtMinusTdb consulted the fit first, so with auto-fit on (the
    // default) an explicit ANALYTIC request was silently served the DE440t fit.
    const bool autofit_was = GetTtTdbAutoFit();
    const TtTdbModel model_was = GetTtTdbModel();
    SetTtTdbAutoFit(true);  // the condition under which the bug appeared
    ClearTtMinusTdbFit();

    SetTtTdbModel(TtTdbModel::ANALYTIC);
    const double analytic = NativeTtMinusTdb(50.0 * day);
    SetTtTdbModel(TtTdbModel::FITTED);
    const double fitted = NativeTtMinusTdb(50.0 * day);

    SetTtTdbModel(model_was);
    SetTtTdbAutoFit(autofit_was);

    // The two models are genuinely different: ANALYTIC must NOT collapse onto
    // the DE440t fit. Their gap is the tens-of-us fidelity difference.
    INFO("analytic = " << analytic << ", fitted = " << fitted);
    REQUIRE(std::abs(analytic - fitted) > 1.0e-6);
  }

  SECTION("outside the fit window, TDBToTt falls back to the analytic series") {
    InitTtMinusTdbFit(Real(0.0), 100.0 * day);

    Real t_outside = 500.0 * day;
    REQUIRE_FALSE(HasFittedTtMinusTdb(t_outside));

    Real fitted_call = TDBToTt(t_outside);  // fit active but doesn't cover t
    ClearTtMinusTdbFit();
    Real analytic_call = TDBToTt(t_outside);  // no fit -> analytic
    REQUIRE_THAT(fitted_call.val(), WithinAbs(analytic_call.val(), 1.0e-15));
  }
}

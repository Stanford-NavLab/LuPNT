#include <lupnt/conversions/frame_conversions.h>
#include <lupnt/conversions/time_conversions.h>
#include <lupnt/core/constants.h>
#include <lupnt/core/file.h>
#include <lupnt/interfaces/eop.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("data.eop") {
  EopFileData* file_data = GetEopFileData();
  REQUIRE(file_data != nullptr);
  REQUIRE(file_data->mjds_utc.size() > 10);
  REQUIRE(file_data->mjds_utc(1) > file_data->mjds_utc(0));

  Real mjd = file_data->mjds_utc(file_data->mjds_utc.size() / 2);
  EopData data = GetEopData(mjd);
  REQUIRE(std::isfinite(data.x_pole.val()));
  REQUIRE(std::isfinite(data.y_pole.val()));
  REQUIRE(std::isfinite(data.ut1_utc.val()));
  REQUIRE(std::isfinite(data.lod.val()));
  REQUIRE_THAT(GetUt1UtcDifference(mjd).val(), WithinAbs(data.ut1_utc.val(), epsilon));
}

// ----------------------------------------------------------------------------
// Coverage reporting and the clamp-outside-the-table behavior. The clamp is the
// dominant Earth-orientation error for simulations run at epochs past the end of
// the loaded file, so it is pinned here rather than left implicit.
// ----------------------------------------------------------------------------
TEST_CASE("data.eop_coverage") {
  EopFileData* file_data = GetEopFileData();
  EopCoverage cov = GetEopCoverage();

  REQUIRE_THAT(cov.mjd_first, WithinAbs(file_data->mjds_utc(0), epsilon));
  REQUIRE_THAT(cov.mjd_last,
               WithinAbs(file_data->mjds_utc(file_data->mjds_utc.size() - 1), epsilon));
  REQUIRE(cov.mjd_last > cov.mjd_first);

  // Past the end of the table: the last row is held constant, not extrapolated.
  int last = file_data->mjds_utc.size() - 1;
  EopData beyond = GetEopData(Real(cov.mjd_last + 365.0));
  REQUIRE_THAT(beyond.x_pole.val(), WithinAbs(file_data->x(last) * RAD_ARCSEC, epsilon));
  REQUIRE_THAT(beyond.y_pole.val(), WithinAbs(file_data->y(last) * RAD_ARCSEC, epsilon));
  REQUIRE_THAT(beyond.ut1_utc.val(), WithinAbs(file_data->ut1_utc(last), epsilon));

  // Before the start: likewise held at the first row.
  EopData before = GetEopData(Real(cov.mjd_first - 365.0));
  REQUIRE_THAT(before.x_pole.val(), WithinAbs(file_data->x(0) * RAD_ARCSEC, epsilon));
  REQUIRE_THAT(before.ut1_utc.val(), WithinAbs(file_data->ut1_utc(0), epsilon));

  // The sigma_* fields are the C04 formal errors, not CIP offsets or TAI-UTC.
  REQUIRE_THAT(before.sigma_x_pole.val(), WithinAbs(file_data->xErr(0) * RAD_ARCSEC, epsilon));
  REQUIRE_THAT(before.sigma_y_pole.val(), WithinAbs(file_data->yErr(0) * RAD_ARCSEC, epsilon));
  REQUIRE_THAT(before.sigma_ut1_utc.val(), WithinAbs(file_data->ut1_utc_err(0), epsilon));
}

// ----------------------------------------------------------------------------
// The IERS finals ("Bulletin A") reader, exercised offline against the finals.all snapshot
// bundled with the repo's orekit data. Restores the bundled C04 table afterwards so the rest of
// the suite (calibrated against that vintage) is unaffected.
// ----------------------------------------------------------------------------
TEST_CASE("data.eop_finals") {
  std::filesystem::path finals = std::filesystem::path(GetDataPath()) / ".." / "orekit-data-main"
                                 / "Earth-Orientation-Parameters" / "IAU-1980" / "finals.all";
  if (!std::filesystem::exists(finals)) {
    WARN("orekit finals.all snapshot not present; skipping finals reader test");
    return;
  }

  LoadEopFinalsFileData(finals, true);
  EopFileData* fd = GetEopFileData();
  EopCoverage cov = GetEopCoverage();

  REQUIRE(fd->source == EopSource::Finals);
  REQUIRE(fd->mjds_utc.size() > 10000);
  REQUIRE(fd->is_prediction.size() == fd->mjds_utc.size());
  // The series begins 1973-01-02 (MJD 41684) and is strictly daily.
  REQUIRE_THAT(cov.mjd_first, WithinAbs(41684.0, epsilon));
  REQUIRE_THAT(fd->mjds_utc(1) - fd->mjds_utc(0), WithinAbs(1.0, epsilon));
  REQUIRE(fd->years(0) == 1973);
  REQUIRE(fd->months(0) == 1);
  REQUIRE(fd->days(0) == 2);

  // First row, checked field by field against the raw file (units: x/y in arcsec, UT1-UTC in s,
  // LOD converted from ms, dpsi/deps converted from mas).
  // 73 1 2 41684.00 I  0.120733 0.009786  0.136966 0.015902  I 0.8084178 0.0002710  0.0000 0.1916
  REQUIRE_THAT(fd->x(0), WithinAbs(0.120733, epsilon));
  REQUIRE_THAT(fd->y(0), WithinAbs(0.136966, epsilon));
  REQUIRE_THAT(fd->ut1_utc(0), WithinAbs(0.8084178, epsilon));
  REQUIRE_THAT(fd->xErr(0), WithinAbs(0.009786, epsilon));
  REQUIRE_THAT(fd->yErr(0), WithinAbs(0.015902, epsilon));
  REQUIRE_THAT(fd->ut1_utc_err(0), WithinAbs(0.0002710, epsilon));
  REQUIRE_THAT(fd->lod_err(0), WithinAbs(0.1916e-3, epsilon));
  REQUIRE(fd->is_prediction(0) == 0);

  // The table must contain both measured and predicted rows, in that order, and the predictions
  // must extend past the measured span -- that reach is the whole point of using finals.
  REQUIRE(cov.mjd_last_measured > cov.mjd_first);
  REQUIRE(cov.mjd_last > cov.mjd_last_measured);
  REQUIRE(fd->is_prediction(0) == 0);
  REQUIRE(fd->is_prediction(fd->is_prediction.size() - 1) == 1);
  for (Eigen::Index i = 1; i < fd->mjds_utc.size(); ++i) {
    REQUIRE(fd->mjds_utc(i) > fd->mjds_utc(i - 1));                         // strictly increasing
    if (fd->is_prediction(i - 1) == 1) REQUIRE(fd->is_prediction(i) == 1);  // no measured-after-P
    REQUIRE((fd->mjds_utc(i) <= cov.mjd_last_measured) == (fd->is_prediction(i) == 0));
  }

  // Values must stay physical across the whole table, predictions included: polar motion within
  // a second of arc, |UT1-UTC| under a second (leap seconds keep it there), LOD sub-millisecond.
  for (Eigen::Index i = 0; i < fd->mjds_utc.size(); ++i) {
    REQUIRE(std::abs(fd->x(i)) < 1.0);
    REQUIRE(std::abs(fd->y(i)) < 1.0);
    REQUIRE(std::abs(fd->ut1_utc(i)) < 1.0);
    REQUIRE(std::abs(fd->lod(i)) < 0.01);
  }

  // Interpolating in the predicted span works and returns finite values.
  EopData predicted = GetEopData(Real(0.5 * (cov.mjd_last_measured + cov.mjd_last)));
  REQUIRE(std::isfinite(predicted.x_pole.val()));
  REQUIRE(std::isfinite(predicted.ut1_utc.val()));
  REQUIRE(std::isfinite(predicted.lod.val()));

  // Finals must agree with the independent C04 series where both are final. Compare at an epoch
  // well inside both (MJD 59000, 2020-05-31); tolerances are a few mas / tens of microseconds,
  // which is the genuine C04-vs-Bulletin-A solution difference.
  Real mjd_common = 59000.0;
  EopData f = GetEopData(mjd_common);
  LoadEopFileData(GetFilePath(EOP_FILENAME), true);
  EopData c = GetEopData(mjd_common);
  REQUIRE_THAT((f.x_pole - c.x_pole).val() / RAD_ARCSEC, WithinAbs(0.0, 5e-3));  // < 5 mas
  REQUIRE_THAT((f.y_pole - c.y_pole).val() / RAD_ARCSEC, WithinAbs(0.0, 5e-3));
  REQUIRE_THAT((f.ut1_utc - c.ut1_utc).val(), WithinAbs(0.0, 1e-4));  // < 0.1 ms
  REQUIRE_THAT((f.lod - c.lod).val(), WithinAbs(0.0, 1e-4));

  // ...and the restored table is C04 again.
  REQUIRE(GetEopFileData()->source == EopSource::C04);
}

// ----------------------------------------------------------------------------
// SetEopSource: the selection must be independent of call order. This is the defect it exists to
// fix -- a bare Load*(path, force=false) issued after anything has touched Earth orientation is
// dropped, so a run can believe it selected one product while carrying another.
// ----------------------------------------------------------------------------
TEST_CASE("data.eop_source_selection") {
  std::filesystem::path finals = std::filesystem::path(GetDataPath()) / ".." / "orekit-data-main"
                                 / "Earth-Orientation-Parameters" / "IAU-1980" / "finals.all";
  if (!std::filesystem::exists(finals)) {
    WARN("orekit finals.all snapshot not present; skipping source selection test");
    return;
  }

  // Default is C04.
  LoadEopFileData(GetFilePath(EOP_FILENAME), true);
  REQUIRE(GetEopSource() == EopSource::C04);
  REQUIRE(GetEopCoverage().source == EopSource::C04);

  // Selecting with a table already loaded switches immediately rather than deferring (which
  // would amount to ignoring the call).
  SetEopSource(EopSource::Finals, finals);
  REQUIRE(GetEopSource() == EopSource::Finals);
  REQUIRE(GetEopCoverage().source == EopSource::Finals);
  REQUIRE(GetEopFileData()->source == EopSource::Finals);

  // Selecting the source already in force is a no-op, not a reload.
  double mjd_last = GetEopCoverage().mjd_last;
  SetEopSource(EopSource::Finals, finals);
  REQUIRE(GetEopCoverage().mjd_last == mjd_last);

  // ...and back.
  SetEopSource(EopSource::C04);
  REQUIRE(GetEopCoverage().source == EopSource::C04);

  // Selecting a source whose file is absent fails at the point of the call, not later inside a
  // frame conversion.
  REQUIRE_THROWS(SetEopSource(EopSource::Finals, finals.parent_path() / "no_such_finals_file"));
  // ...and a failed selection leaves the previous choice intact.
  REQUIRE(GetEopSource() == EopSource::C04);
  REQUIRE(GetEopCoverage().source == EopSource::C04);

  // A non-forcing load with a table already present is dropped -- the behavior SetEopSource
  // exists to route around. Pinned here so it cannot regress into a silent switch.
  LoadEopFinalsFileData(finals, false);
  REQUIRE(GetEopCoverage().source == EopSource::C04);

  LoadEopFileData(GetFilePath(EOP_FILENAME), true);
}

// ----------------------------------------------------------------------------
// Loading the latest EOP series from the IERS data center. This test does not
// require network access: if the download fails (offline environment), it
// just verifies that the bundled EOP data remains usable.
// ----------------------------------------------------------------------------
TEST_CASE("data.eop_latest_iers") {
  double mjd_last_before = GetEopFileData()->mjds_utc(GetEopFileData()->mjds_utc.size() - 1);

  bool downloaded = LoadLatestEopFromIers(true);

  EopFileData* file_data = GetEopFileData();
  REQUIRE(file_data != nullptr);
  REQUIRE(file_data->mjds_utc.size() > 10);

  if (downloaded) {
    double mjd_last_after = file_data->mjds_utc(file_data->mjds_utc.size() - 1);
    REQUIRE(mjd_last_after > mjd_last_before);

    // 2025-01-01 in MJD, well within the latest IERS series.
    Real mjd_recent = 60676.0;
    EopData data = GetEopData(mjd_recent);
    REQUIRE(std::isfinite(data.x_pole.val()));
    REQUIRE(std::isfinite(data.y_pole.val()));
    REQUIRE(std::isfinite(data.ut1_utc.val()));

    // Restore the bundled EOP data so other tests (e.g. ITRF cross-validation,
    // which is calibrated against the bundled EOP vintage) are unaffected.
    LoadEopFileData(GetFilePath(EOP_FILENAME), true);
  } else {
    Real mjd = file_data->mjds_utc(file_data->mjds_utc.size() / 2);
    EopData data = GetEopData(mjd);
    REQUIRE(std::isfinite(data.x_pole.val()));
    REQUIRE(std::isfinite(data.ut1_utc.val()));
  }
}

// ----------------------------------------------------------------------------
// EOP perturbation injector (Step 2a). Pins the small-angle convention by central difference
// rather than asserting it -- the signs of eps are exactly the kind of thing that is easy to get
// backwards and impossible to notice downstream.
//
//   delta r_gcrf = R_itrf->gcrf * (eps x r_itrf),
//   eps = (-dy_pole, -dx_pole, OMEGA_ERA * dut1)   [rad, in ITRF]
//
// Note the UT1 step: LuPNT epochs are seconds since J2000 (~1.9e9), so a 1e-6 s step is at
// double-precision resolution and the difference is pure roundoff. 1e-2 s is needed.
// ----------------------------------------------------------------------------
TEST_CASE("data.eop_perturbation") {
  LoadEopFileData(GetFilePath(EOP_FILENAME), true);

  // No perturbation by default, and the fast-path flag is off for a C04 table.
  ClearEopPerturbation();
  REQUIRE(EopHasCelestialPoleOffsets() == false);
  EopPerturbation none = GetEopPerturbation(Real(59000.0));
  REQUIRE(none.dut1.val() == 0.0);

  Real mjd = 59000.0;
  EopData nominal = GetEopData(mjd);

  // A constant perturbation must land additively on every field it names.
  EopPerturbation p;
  p.dx_pole = 1e-8;
  p.dy_pole = -2e-8;
  p.dut1 = 3e-3;
  p.dlod = 4e-6;
  SetEopPerturbation(p);
  REQUIRE(EopHasCelestialPoleOffsets() == true);
  EopData perturbed = GetEopData(mjd);
  REQUIRE_THAT((perturbed.x_pole - nominal.x_pole).val(), WithinAbs(1e-8, 1e-14));
  REQUIRE_THAT((perturbed.y_pole - nominal.y_pole).val(), WithinAbs(-2e-8, 1e-14));
  REQUIRE_THAT((perturbed.ut1_utc - nominal.ut1_utc).val(), WithinAbs(3e-3, 1e-9));
  REQUIRE_THAT((perturbed.lod - nominal.lod).val(), WithinAbs(4e-6, 1e-12));

  // A time-varying perturbation is evaluated at the epoch requested, not at a fixed one.
  ClearEopPerturbation();
  EopData nominal_10 = GetEopData(Real(59010.0));
  SetEopPerturbationFunction([](Real m) {
    EopPerturbation q;
    q.dut1 = (m - 59000.0) * 1e-3;
    return q;
  });
  REQUIRE_THAT((GetEopData(Real(59010.0)).ut1_utc - nominal_10.ut1_utc).val(),
               WithinAbs(1e-2, 1e-9));
  REQUIRE_THAT((GetEopData(Real(59000.0)).ut1_utc - nominal.ut1_utc).val(), WithinAbs(0.0, 1e-12));
  ClearEopPerturbation();
  REQUIRE(EopHasCelestialPoleOffsets() == false);
  REQUIRE(GetEopPerturbation(Real(59010.0)).dut1.val() == 0.0);

  // ---- the small-angle convention, by central difference on GCRF<->ITRF ----
  Real t_tdb = ConvertTime(MjdToTime(mjd), Time::UTC, Time::TDB);
  Vec3 r_itrf(4197160.8, 815845.4, 4716876.3);

  auto rotated = [&](const EopPerturbation& q) {
    SetEopPerturbation(q);
    Vec3 out = ItrfToGcrf(t_tdb, r_itrf);
    ClearEopPerturbation();
    return out;
  };
  Vec3 r_nominal = rotated(EopPerturbation{});
  // No public accessor for the composite rotation, so build it from the public conversion:
  // ITRF->GCRF is a pure rotation (both geocentric), so its columns are the images of the axes.
  Mat3 R;
  R.col(0) = ItrfToGcrf(t_tdb, Vec3(1, 0, 0));
  R.col(1) = ItrfToGcrf(t_tdb, Vec3(0, 1, 0));
  R.col(2) = ItrfToGcrf(t_tdb, Vec3(0, 0, 1));

  struct Case {
    const char* name;
    double step;
    Vec3 axis;  // eps per unit of the parameter, in ITRF
    EopPerturbation (*make)(double);
  };
  const double OMEGA_ERA = TWO_PI * 1.00273781191135448 / SECS_DAY;
  std::vector<Case> cases{
      {"dx_pole", 1e-8, Vec3(0, -1, 0),
       [](double h) {
         EopPerturbation q;
         q.dx_pole = h;
         return q;
       }},
      {"dy_pole", 1e-8, Vec3(-1, 0, 0),
       [](double h) {
         EopPerturbation q;
         q.dy_pole = h;
         return q;
       }},
      {"dut1", 1e-2, Vec3(0, 0, 1),
       [](double h) {
         EopPerturbation q;
         q.dut1 = h;
         return q;
       }},
  };
  for (const Case& c : cases) {
    Vec3 num = (rotated(c.make(c.step)) - rotated(c.make(-c.step))) / (2.0 * c.step);
    Vec3 eps = c.axis;
    if (std::string(c.name) == "dut1") eps *= OMEGA_ERA;
    Vec3 ana = R * eps.cross(r_itrf);
    double rel = (num - ana).norm().val() / ana.norm().val();
    INFO("parameter " << c.name << " rel.err " << rel);
    REQUIRE(rel < 1e-4);
  }

  // The unperturbed rotation must be exactly restored.
  REQUIRE_THAT((rotated(EopPerturbation{}) - r_nominal).norm().val(), WithinAbs(0.0, 1e-9));
  ClearEopPerturbation();
}

#include <lupnt/core/constants.h>
#include <lupnt/core/file.h>
#include <lupnt/data/eop.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

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

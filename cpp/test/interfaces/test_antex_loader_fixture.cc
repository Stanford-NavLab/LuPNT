#include <lupnt/lupnt.h>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  void RequireNearVec3d(const Vec3d& a, const Vec3d& b, double tol) {
    for (int i = 0; i < 3; i++) REQUIRE_THAT(a(i), WithinAbs(b(i), tol));
  }

  // Small, trimmed ANTEX file committed under `cpp/test/fixtures/gnss/` (full
  // header + the GPS G01 antenna block valid from 2024-12-17 onward, plus one
  // extra GPS block) so the loader gets CI coverage without the ~51 MB product.
  // `LUPNT_TEST_FIXTURES_DIR` is defined by `cpp/test/CMakeLists.txt`.
  std::filesystem::path FixtureFile() {
    return std::filesystem::path(LUPNT_TEST_FIXTURES_DIR) / "gnss" / "igs20_trimmed.atx";
  }
}  // namespace

const double epsilon = 1e-6;

TEST_CASE("interfaces.antex_fixture") {
  const std::filesystem::path antex_file = FixtureFile();

  // GPS PRN 01 (BLOCK IIIA, SVN G080) is valid from 2024-12-17 onward (open-ended);
  // this matches the epoch of the SP3/RINEX-nav test data (2026-01-14).
  Real t_tai = GregorianToTime(2026, 1, 14, 0, 0, 0);

  // SP3 (center-of-mass) ECEF position sample for G01 at the same epoch.
  const Vec3d r_sat_ecef(21175.826427e3, 13096.416599e3, 9341.310042e3);

  SECTION("Default-constructed loader has no satellites loaded") {
    AntexLoader antex;
    REQUIRE_FALSE(antex.HasSatellite(GnssConst::GPS, 1));
  }

  SECTION("Single-file constructor parses satellite antenna entries") {
    AntexLoader antex(antex_file);

    REQUIRE(antex.HasSatellite(GnssConst::GPS, 1));
    REQUIRE_FALSE(antex.HasSatellite(GnssConst::GPS, 99));

    std::vector<std::string> codes = antex.GetAvailableFreqCodes(GnssConst::GPS, 1, t_tai);
    REQUIRE_FALSE(codes.empty());
    REQUIRE(std::find(codes.begin(), codes.end(), "G01") != codes.end());

    // GnssFreq overload and raw-ANTEX-frequency-code overload of GetPco agree
    Vec3d pco_l1 = antex.GetPco(GnssConst::GPS, 1, GnssFreq::L1, t_tai);
    Vec3d pco_g01 = antex.GetPco("G", 1, "G01", t_tai);
    RequireNearVec3d(pco_l1, pco_g01, epsilon);

    // PCO magnitude is on the order of meters (NEU antenna frame, ANTEX stores mm)
    REQUIRE(pco_l1.norm() > 0.1);
    REQUIRE(pco_l1.norm() < 100.0);

    // Querying an unknown satellite / frequency throws
    REQUIRE_THROWS(antex.GetPco(GnssConst::GPS, 99, GnssFreq::L1, t_tai));
    REQUIRE_THROWS(antex.GetPco(GnssConst::GPS, 1, GnssFreq::E1, t_tai));
  }

  SECTION("LoadFile merges additional satellite entries") {
    AntexLoader antex;
    REQUIRE_FALSE(antex.HasSatellite(GnssConst::GPS, 1));
    antex.LoadFile(antex_file);
    REQUIRE(antex.HasSatellite(GnssConst::GPS, 1));
  }

  SECTION("ComputeIjkToEcefRotation points kvec toward nadir") {
    Mat3d Cijk = AntexLoader::ComputeIjkToEcefRotation(t_tai, r_sat_ecef);

    // kvec (3rd column) = -normalize(r_sat_ecef)
    RequireNearVec3d(Cijk.col(2), Vec3d(-r_sat_ecef.normalized()), epsilon);

    REQUIRE_THAT(Cijk.col(0).norm(), WithinAbs(Cijk.col(1).cross(Cijk.col(2)).norm(), epsilon));
    RequireNearVec3d(Cijk.col(0), Vec3d(Cijk.col(1).cross(Cijk.col(2))), epsilon);
  }

  SECTION("ApplyPcoCorrectionEcef offsets the SP3 position by Cijk * pco_neu") {
    AntexLoader antex(antex_file);
    Vec3d pco_neu_m = antex.GetPco(GnssConst::GPS, 1, GnssFreq::L1, t_tai);

    Vec3d corrected = AntexLoader::ApplyPcoCorrectionEcef(t_tai, r_sat_ecef, pco_neu_m);
    Mat3d Cijk = AntexLoader::ComputeIjkToEcefRotation(t_tai, r_sat_ecef);
    Vec3d expected = r_sat_ecef + Cijk * pco_neu_m;
    RequireNearVec3d(corrected, expected, epsilon);

    REQUIRE((corrected - r_sat_ecef).norm() > 0.0);
  }

  SECTION("GnssLetter / SatId identifier helpers") {
    REQUIRE(AntexLoader::GnssLetter(GnssConst::GPS) == "G");
    REQUIRE(AntexLoader::GnssLetter(GnssConst::GLONASS) == "R");
    REQUIRE(AntexLoader::GnssLetter(GnssConst::GALILEO) == "E");
    REQUIRE(AntexLoader::GnssLetter(GnssConst::BEIDOU) == "C");
    REQUIRE(AntexLoader::GnssLetter(GnssConst::QZSS) == "J");

    REQUIRE(AntexLoader::SatId(GnssConst::GPS, 1) == "G01");
    REQUIRE(AntexLoader::SatId(GnssConst::GPS, 32) == "G32");
    REQUIRE(AntexLoader::SatId(GnssConst::GALILEO, 11) == "E11");
    REQUIRE(AntexLoader::SatId(GnssConst::BEIDOU, 6) == "C06");
  }
}

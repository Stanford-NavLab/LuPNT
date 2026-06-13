#include <lupnt/lupnt.h>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <filesystem>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  // See `test_sp3_loader.cc` for why this differs from `GetOutputDir("gnss_files")`.
  std::filesystem::path GnssFilesDir() {
    return GetDataPath().parent_path().parent_path() / "output" / "gnss_files";
  }

  // `GnssConstellation` consumes precomputed ECI ephemerides -- it does not
  // propagate orbits itself -- so a simple analytic equatorial circular orbit
  // is sufficient to exercise its geometry / link-budget / measurement APIs.
  MatXd CircularOrbitEci(double radius_m, double phase0_rad, const VecXd& t_tai) {
    double omega = std::sqrt(GM_EARTH / (radius_m * radius_m * radius_m));
    MatXd rv(t_tai.size(), 6);
    for (int k = 0; k < t_tai.size(); k++) {
      double th = phase0_rad + omega * t_tai(k);
      rv(k, 0) = radius_m * std::cos(th);
      rv(k, 1) = radius_m * std::sin(th);
      rv(k, 2) = 0.0;
      rv(k, 3) = -radius_m * omega * std::sin(th);
      rv(k, 4) = radius_m * omega * std::cos(th);
      rv(k, 5) = 0.0;
    }
    return rv;
  }
}  // namespace

const double epsilon = 1e-6;

TEST_CASE("agents.gnss_constellation") {
  const double t0 = GregorianToTime(2026, 1, 14, 0, 0, 0).val();
  const double gps_radius_m = 26560e3;

  VecXd t_tai(5);
  for (int k = 0; k < 5; k++) t_tai(k) = t0 + k * 3600.0;  // hourly samples over 4 hours

  const std::vector<int> prns = {1, 2};
  const std::vector<MatXd> rv_eci_list
      = {CircularOrbitEci(gps_radius_m, 0.0, t_tai), CircularOrbitEci(gps_radius_m, PI, t_tai)};

  // ---- Constructors --------------------------------------------------------

  SECTION("Default / GnssConst / Config constructors") {
    GnssConstellation c_default;
    REQUIRE(c_default.GetNumSatellites() == 0);
    REQUIRE(c_default.GetGnssConst() == GnssConst::GPS);
    REQUIRE(c_default.GetPrns().empty());

    GnssConstellation c_galileo(GnssConst::GALILEO);
    REQUIRE(c_galileo.GetGnssConst() == GnssConst::GALILEO);

    Config config = YAML::Load("gnss_const: GPS\nname: TestGps\n");
    GnssConstellation c_config(config);
    REQUIRE(c_config.GetGnssConst() == GnssConst::GPS);

    c_config.SetSatelliteStates(prns, t_tai, rv_eci_list);
    REQUIRE(c_config.GetNumSatellites() == 2);
  }

  // ---- SetSatelliteStates / basic queries -----------------------------------

  GnssConstellation constellation(GnssConst::GPS);
  constellation.SetSatelliteStates(prns, t_tai, rv_eci_list);

  SECTION("SetSatelliteStates populates PRNs / GetNumSatellites / GetPrns") {
    REQUIRE(constellation.GetNumSatellites() == 2);
    REQUIRE(constellation.GetPrns() == prns);
    REQUIRE(constellation.GetGnssConst() == GnssConst::GPS);
  }

  SECTION("SetFaultPrns / IsFaultPrn") {
    REQUIRE_FALSE(constellation.IsFaultPrn(1));
    REQUIRE_FALSE(constellation.IsFaultPrn(2));

    constellation.SetFaultPrns({2});
    REQUIRE_FALSE(constellation.IsFaultPrn(1));
    REQUIRE(constellation.IsFaultPrn(2));
  }

  SECTION("GetSatelliteStateEci interpolates the precomputed ephemeris") {
    const double interp_epsilon = 1e-3;

    // At a sample epoch, the interpolated state matches the input table exactly
    Vec6 s0 = constellation.GetSatelliteStateEci(1, t_tai(0));
    Vec6 s0_ref = rv_eci_list[0].row(0).transpose().cast<Real>();
    RequireNear(s0, s0_ref, interp_epsilon);

    // The interpolated state at the bracketing epoch should equal the average
    // of the two samples (GnssConstellation uses linear interpolation).
    Real t_mid = Real(0.5 * (t_tai(0) + t_tai(1)));
    Vec6 s_mid = constellation.GetSatelliteStateEci(1, t_mid);
    Vec6 s_avg = ((rv_eci_list[0].row(0) + rv_eci_list[0].row(1)) / 2.0).transpose().cast<Real>();
    RequireNear(s_mid, s_avg, interp_epsilon);

    REQUIRE_THROWS(constellation.GetSatelliteStateEci(99, t_tai(0)));
  }

  SECTION("ComputeRange / ComputeRangeRate against a receiver state") {
    Vec3 r_rx(384400e3, 0.0, 0.0);  // ~lunar distance, on the +x axis
    Vec6 rv_rx;
    rv_rx << r_rx, Vec3(0.0, 1000.0, 0.0);

    Real t = Real(t_tai(0));
    Vec6 rv_tx = constellation.GetSatelliteStateEci(1, t);

    Real range = constellation.ComputeRange(1, r_rx, t);
    RequireNear(range, Real((r_rx - rv_tx.head(3)).norm()), epsilon);
    REQUIRE(range.val() > 0.0);

    Real range_rate = constellation.ComputeRangeRate(1, rv_rx, t);
    Vec3 dr = rv_rx.head(3) - rv_tx.head(3);
    Vec3 dv = rv_rx.tail(3) - rv_tx.tail(3);
    RequireNear(range_rate, Real(dr.dot(dv) / dr.norm()), epsilon);
  }

  // ---- Ephemeris HDF5 round trip --------------------------------------------

  SECTION("SaveEphemeris / LoadEphemeris round trip via HDF5") {
    std::filesystem::path h5_path
        = std::filesystem::temp_directory_path() / "lupnt_test_gnss_constellation_ephemeris.h5";

    constellation.SaveEphemeris(h5_path);
    REQUIRE(std::filesystem::exists(h5_path));

    GnssConstellation loaded(GnssConst::GPS);
    loaded.LoadEphemeris(h5_path);

    REQUIRE(loaded.GetNumSatellites() == constellation.GetNumSatellites());
    REQUIRE(loaded.GetPrns() == constellation.GetPrns());
    for (int prn : prns) {
      Vec6 s_orig = constellation.GetSatelliteStateEci(prn, t_tai(2));
      Vec6 s_load = loaded.GetSatelliteStateEci(prn, t_tai(2));
      RequireNear(s_orig, s_load, epsilon);
    }

    std::filesystem::remove(h5_path);
  }

  // ---- Attitude (static helper) ---------------------------------------------

  SECTION("ComputeAttitude (static) delegates to GnssAttitude::Compute") {
    Vec3 r_sat(gps_radius_m, 0.0, 0.0);
    Vec3 r_sun(0.0, AU, 0.0);

    Vec3 ex1, ey1, ez1, ex2, ey2, ez2;
    GnssConstellation::ComputeAttitude(r_sat, r_sun, ex1, ey1, ez1);
    GnssAttitude::Compute(r_sat, r_sun, ex2, ey2, ez2);

    RequireNear(ex1, ex2, epsilon);
    RequireNear(ey1, ey2, epsilon);
    RequireNear(ez1, ez2, epsilon);
  }

  // ---- Transmitters / measurement noise -------------------------------------

  SECTION("SetupTransmitters loads antenna patterns & transmit power for each PRN") {
    constellation.SetupTransmitters();
    REQUIRE_FALSE(constellation.GetFreqList().empty());

    // GPS PRN 1 (Block IIF) transmits on L1/L2/L5
    const auto& freqs = constellation.GetFreqList();
    REQUIRE(std::find(freqs.begin(), freqs.end(), GnssFreq::L1) != freqs.end());
    REQUIRE(constellation.HasTransmitterInfo(1, GnssFreq::L1));
    REQUIRE(std::isfinite(constellation.GetTransmitPowerDbw(1, GnssFreq::L1).val()));
  }

  SECTION("measurement-noise sigmas") {
    Real cn0_dbhz = Real(45.0);
    Real sigma_range = constellation.ComputeSigmaRange(cn0_dbhz, GnssFreq::L1);
    Real sigma_rate = constellation.ComputeSigmaRangeRate(cn0_dbhz, GnssFreq::L1);
    Real sigma_phase = constellation.ComputeSigmaCarrierPhase(cn0_dbhz, GnssFreq::L1);
    REQUIRE(sigma_range.val() > 0.0);
    REQUIRE(sigma_rate.val() > 0.0);
    REQUIRE(sigma_phase.val() > 0.0);

    // Higher C/N0 should produce smaller (or equal) measurement-noise std devs
    Real cn0_better = Real(55.0);
    REQUIRE(constellation.ComputeSigmaRange(cn0_better, GnssFreq::L1).val() <= sigma_range.val());
    REQUIRE(constellation.ComputeSigmaRangeRate(cn0_better, GnssFreq::L1).val()
            <= sigma_rate.val());
    REQUIRE(constellation.ComputeSigmaCarrierPhase(cn0_better, GnssFreq::L1).val()
            <= sigma_phase.val());
  }

  // ---- Constellation interface ----------------------------------------------

  SECTION("Setup / Step") {
    // Setup() lazily calls SetupTransmitters() the first time
    GnssConstellation c2(GnssConst::GPS);
    c2.SetSatelliteStates(prns, t_tai, rv_eci_list);
    REQUIRE(c2.GetFreqList().empty());
    c2.Setup();
    REQUIRE_FALSE(c2.GetFreqList().empty());

    // Step() is a (precomputed-ephemeris) no-op; just verify it doesn't throw
    c2.Step(Real(t_tai(0)));
  }
}

TEST_CASE("agents.gnss_constellation.setup_from_files") {
  const std::filesystem::path sp3_file
      = GnssFilesDir() / "sp3" / "COD0MGXFIN_20260140000_01D_05M_ORB.SP3";
  const std::filesystem::path antex_file = GetDataPath() / "gnss" / "igs20.atx";

  Sp3Loader sp3_probe(sp3_file);
  auto [t_min, t_max] = sp3_probe.GetTimeSpan("G01");
  double span = t_max - t_min;

  VecXd t_tai(3);
  t_tai << t_min + 0.1 * span, t_min + 0.5 * span, t_min + 0.9 * span;

  GnssConstellation constellation(GnssConst::GPS);
  constellation.SetupSatelliteStatesFromFiles({sp3_file}, antex_file, t_tai, GnssFreq::L1, {1});

  REQUIRE(constellation.GetNumSatellites() == 1);
  REQUIRE(constellation.GetPrns() == std::vector<int>{1});

  Vec6 rv_eci = constellation.GetSatelliteStateEci(1, t_tai(1));
  Real r = rv_eci.head(3).norm();
  REQUIRE(r.val() > 2.0e7);  // GPS-orbit-altitude geocentric radius
  REQUIRE(r.val() < 3.0e7);

  // The PCO correction should shift the antenna phase center away from the SP3
  // (center-of-mass) position by approximately the PCO magnitude (a few meters).
  Vec6 rv_sp3_ecef = sp3_probe.GetPosVel("G01", Real(t_tai(1)));
  Real t_tdb = ConvertTime(Real(t_tai(1)), Time::TAI, Time::TDB);
  Vec6 rv_sp3_eci = ConvertFrame(t_tdb, rv_sp3_ecef, Frame::ECEF, Frame::ECI, false);
  Real delta = (rv_eci.head(3) - rv_sp3_eci.head(3)).norm();

  AntexLoader antex_probe(antex_file);
  Vec3d pco = antex_probe.GetPco(GnssConst::GPS, 1, GnssFreq::L1, Real(t_tai(1)));
  REQUIRE_THAT(delta.val(), WithinAbs(pco.norm(), 1.0));
}

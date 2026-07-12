// Example 5: GNSS Measurement Simulation for a Lunar Satellite
// ------------------------------------------------------------------------------
// C++ counterpart of python/examples/ex5_gnss_measurement_sim.ipynb.
//
// Simulates which Earth-GNSS signals a lunar-orbiting receiver can track. It
// mirrors the notebook flow but reports numerically instead of plotting:
//
//   1. define + propagate a lunar ELFO receiver (one orbital period, N-body);
//   2. load GPS and Galileo precise ephemerides from SP3 (via
//      `GnssConstellation::SetupSatelliteStatesFromFiles`, which uses
//      `Sp3Loader`/`AntexLoader` and applies the antenna phase-center offset);
//   3. build a `GNSSMeasurements` model (GPS L1 + Galileo E1, Earth + Moon
//      occultation, moongpsr Earth-pointing receiver antenna);
//   4. `Precompute` visibility + link-budget (C/N0) at every receiver epoch;
//   5. print per-epoch visible-satellite counts and C/N0 statistics.
//
// The basic observable is pseudorange; this example focuses on the geometry and
// link budget (which links close, and how strong they are) that later filters
// such as Example 6 consume.
//
// Data note: the SP3 products are downloaded from NASA CDDIS (Earthdata Login;
// ~/.netrc or EARTHDATA_USERNAME/PASSWORD, see docs/pages/sp3_download) and the
// bundled IGS20 ANTEX file is used for the phase-center offset. Absent the
// data/credentials the program prints an informative note and exits cleanly,
// so it always compiles and runs.

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "lupnt/lupnt.h"

using namespace lupnt;

int main() {
  std::cout << std::fixed;

  // ==== 1. Epoch ============================================================
  const Real t0_tdb = GregorianToTime(2026, 1, 1, 0, 0, 0.0);
  const Real t0_tai = ConvertTime(t0_tdb, Time::TDB, Time::TAI);
  std::cout << "Epoch (TDB): " << TimeToGregorianString(t0_tdb) << "\n";

  // ==== 2. ELFO receiver ====================================================
  const Real a = 6541.4e3;
  const Real e = 0.6;
  const Real i = 56.2 * RAD;
  const Real raan = 0.0 * RAD;
  const Real w = 90.0 * RAD;
  const Real M0 = 0.0 * RAD;
  Vec6 coe0_op(a, e, i, raan, w, M0);
  const Real period = GetOrbitalPeriod(a, GM_MOON);
  std::cout << "Orbital period: " << std::setprecision(1) << period / SECS_HOUR << " h\n";

  const Vec6 rv0_op = ClassicalToCart(coe0_op, GM_MOON);
  const Vec6 rv0_mci = ConvertFrame(t0_tdb, rv0_op, Frame::MOON_OP, Frame::MOON_CI);

  // ==== 3. N-body dynamics ==================================================
  NBodyDynamics dyn;
  dyn.SetFrame(Frame::MOON_CI);
  dyn.SetUnits(SI_UNITS);
  dyn.AddBody(Body::Moon(2, 2));
  dyn.AddBody(Body::Earth());
  dyn.AddBody(Body::Sun());
  dyn.SetSrpCoefficient(1.8, 0.02, 10.0);
  dyn.SetIntegrator(IntegratorType::RKF45);
  IntegratorParams params;
  params.reltol = 1e-10;
  params.abstol = 1e-6;
  dyn.SetIntegratorParams(params);

  // ==== 4. Propagate one orbital period =====================================
  const double dt_rx = 300.0;
  const int n_rx = static_cast<int>(period.val() / dt_rx) + 1;
  VecX t_tdb_arr(n_rx), t_tai_arr(n_rx);
  for (int k = 0; k < n_rx; ++k) {
    t_tdb_arr(k) = t0_tdb + k * dt_rx;
    t_tai_arr(k) = ConvertTime(t_tdb_arr(k), Time::TDB, Time::TAI);
  }
  const MatX6 rv_mci = dyn.Propagate(rv0_mci, t_tdb_arr);
  const MatX6 rv_eci = ConvertFrame(t_tdb_arr, rv_mci, Frame::MOON_CI, Frame::ECI);
  {
    double alt_max = -1e300, alt_min = 1e300;
    for (int k = 0; k < n_rx; ++k) {
      const double alt = (rv_mci.row(k).head(3).norm().val() - R_MOON) / 1e3;
      alt_max = std::max(alt_max, alt);
      alt_min = std::min(alt_min, alt);
    }
    std::cout << "Propagated " << n_rx << " epochs (" << dt_rx << "-s steps): apoapsis "
              << std::setprecision(0) << alt_max << " km, periapsis " << alt_min
              << " km above the surface\n";
  }

  // 8-element receiver states [r, v, clock_bias=0, clock_drift=0].
  std::vector<State> rx_states;
  std::vector<Real> receive_times;
  rx_states.reserve(n_rx);
  receive_times.reserve(n_rx);
  for (int k = 0; k < n_rx; ++k) {
    State s(8);
    s.setZero();
    s.head(6) = rv_eci.row(k).transpose();
    rx_states.push_back(s);
    receive_times.push_back(t_tai_arr(k));
  }

  // ==== 5-8. Constellations + measurement model (data-driven) ===============
  try {
    // ---- Download SP3 covering the arc (+/- 1 day margin) ------------------
    const double sim_t_s = period.val() + 3600.0;
    const int n_days = static_cast<int>(std::ceil(sim_t_s / 86400.0));
    const double base_day_tai
        = ConvertTime(GregorianToTime(2026, 1, 1, 12, 0, 0), Time::TDB, Time::TAI).val();
    std::cout << "\nDownloading SP3 (COD MGEX final) from CDDIS ...\n";
    std::vector<std::filesystem::path> sp3_paths;
    for (int d = -1; d <= n_days; ++d)
      sp3_paths.push_back(Sp3Loader::DownloadFileForEpoch(base_day_tai + d * 86400.0, Time::TAI));
    const std::filesystem::path antex_path = GetDataPath() / "gnss" / "igs20.atx";

    // ---- Transmitter sample grid (one step before t0 for light-time) ------
    const double dt_sat = 300.0;
    const double margin_s = dt_sat;
    const int n_sat = static_cast<int>((sim_t_s + margin_s) / dt_sat) + 1;
    VecXd t_sat_tai(n_sat);
    for (int k = 0; k < n_sat; ++k) t_sat_tai(k) = (t0_tai.val() - margin_s) + k * dt_sat;

    // ---- Build GPS + Galileo constellations from the precise files --------
    auto gps_const = MakePtr<GnssConstellation>(GnssConst::GPS);
    gps_const->SetupSatelliteStatesFromFiles(sp3_paths, antex_path, t_sat_tai, GnssFreq::L1);
    gps_const->SetupTransmitters();
    std::cout << "GnssConstellation (GPS):     " << gps_const->GetNumSatellites()
              << " satellites\n";

    auto gal_const = MakePtr<GnssConstellation>(GnssConst::GALILEO);
    gal_const->SetupSatelliteStatesFromFiles(sp3_paths, antex_path, t_sat_tai, GnssFreq::E1);
    gal_const->SetupTransmitters();
    std::cout << "GnssConstellation (Galileo): " << gal_const->GetNumSatellites()
              << " satellites\n";

    // ---- Configure GNSSMeasurements (GPS L1 + Galileo E1) -----------------
    GnssMeasurementOptions opts;
    opts.frame = Frame::ECI;
    opts.receive_time_scale = Time::TAI;
    opts.ephemeris_time_scale = Time::TAI;
    opts.solve_light_time = true;
    opts.apply_transmitter_relativity = true;
    opts.apply_shapiro_delay = false;
    opts.apply_visibility = true;
    opts.apply_cn0_threshold = true;
    opts.cn0_acquisition_threshold_dbhz = 22.0;
    opts.cn0_tracking_threshold_dbhz = 20.0;
    opts.apply_ionosphere_plasma_delay = false;

    GnssReceiverParams rx_params;
    rx_params.Bp = 1.0;
    rx_params.Bn = 0.7;
    rx_params.T = 20e-3;

    GnssOccludingBody earth_body;
    earth_body.radius_m = R_EARTH;
    earth_body.position_m = Vec3::Zero();

    GnssOccludingBody moon_body;
    moon_body.radius_m = R_MOON;
    moon_body.position_provider = [](Real t_tai) {
      const Real t_tdb = ConvertTime(t_tai, Time::TAI, Time::TDB);
      return GetBodyPos(t_tdb, BodyId::EARTH, BodyId::MOON, Frame::ECI);
    };

    GNSSMeasurements gnss_meas(gps_const);
    gnss_meas.SetFrequency(GnssFreq::L1);
    gnss_meas.AddConstellation(gal_const, GnssFreq::E1);
    gnss_meas.SetOptions(opts);
    gnss_meas.SetReceiverParams(rx_params);
    gnss_meas.SetReceiverAntenna(Antenna("moongpsr"));  // Earth-pointing, +14 dBi
    gnss_meas.SetOccludingBodies({earth_body, moon_body});
    gnss_meas.SetSunPositionProvider([](Real t_tai) {
      const Real t_tdb = ConvertTime(t_tai, Time::TAI, Time::TDB);
      return GetBodyPos(t_tdb, BodyId::EARTH, BodyId::SUN, Frame::ECI);
    });
    gnss_meas.SetBoresightTargetProvider([](Real) { return Vec3::Zero(); });  // Earth at ECI origin

    // ==== Precompute visibility + C/N0 over the arc =========================
    std::cout << "Running Precompute over " << n_rx << " epochs ...\n";
    const std::vector<GNSSMeasurementsEpoch>& epochs
        = gnss_meas.Precompute(receive_times, rx_states, /*compute_jacobians=*/false);

    // ==== Extract results ===================================================
    int n_gps_links = 0, n_gal_links = 0, total_links = 0;
    int peak_vis = 0, min_vis = 1 << 30;
    long sum_vis = 0;
    std::vector<double> cn0_all;
    for (const GNSSMeasurementsEpoch& ep : epochs) {
      const int nvis = static_cast<int>(ep.channels.size());
      sum_vis += nvis;
      peak_vis = std::max(peak_vis, nvis);
      min_vis = std::min(min_vis, nvis);
      for (const GnssChannel& ch : ep.channels) {
        ++total_links;
        if (ch.gnss_const == GnssConst::GPS)
          ++n_gps_links;
        else if (ch.gnss_const == GnssConst::GALILEO)
          ++n_gal_links;
        const double cn0 = ch.cn0_dbhz.val();
        if (std::isfinite(cn0)) cn0_all.push_back(cn0);
      }
    }

    std::cout << "\n=== Visibility summary (C/N0 >= 22 dB-Hz to acquire) ===\n";
    std::cout << "  Total visible links : " << total_links << "  (GPS " << n_gps_links
              << ", Galileo " << n_gal_links << ")\n";
    std::cout << std::setprecision(1);
    std::cout << "  Visible sats / epoch: mean " << static_cast<double>(sum_vis) / n_rx << ", peak "
              << peak_vis << ", min " << min_vis << "\n";
    if (!cn0_all.empty()) {
      const auto [mn, mx] = std::minmax_element(cn0_all.begin(), cn0_all.end());
      double mean = 0.0;
      for (double c : cn0_all) mean += c;
      mean /= static_cast<double>(cn0_all.size());
      std::cout << "  C/N0 [dB-Hz]        : min " << *mn << ", mean " << mean << ", max " << *mx
                << "\n";
    }

    // Per-epoch visible-count time history (subsampled to keep the table short).
    std::cout << "\n  Visible-satellite time history (every ~1 h):\n";
    std::cout << "    " << std::setw(8) << "t[h]" << std::setw(8) << "GPS" << std::setw(10)
              << "Galileo" << std::setw(8) << "total" << "\n";
    const int stride = std::max(1, static_cast<int>(3600.0 / dt_rx));
    for (int k = 0; k < static_cast<int>(epochs.size()); k += stride) {
      int ng = 0, ne = 0;
      for (const GnssChannel& ch : epochs[k].channels) {
        if (ch.gnss_const == GnssConst::GPS)
          ++ng;
        else if (ch.gnss_const == GnssConst::GALILEO)
          ++ne;
      }
      std::cout << "    " << std::setw(8) << std::setprecision(2) << (k * dt_rx / 3600.0)
                << std::setw(8) << ng << std::setw(10) << ne << std::setw(8) << (ng + ne) << "\n";
    }

  } catch (const std::exception& ex) {
    std::cout << "\n[skipped] SP3/ANTEX data unavailable: " << ex.what() << "\n"
              << "This example needs date-specific SP3 products from NASA CDDIS (Earthdata\n"
              << "Login) and the bundled IGS20 ANTEX file. See docs/pages/sp3_download.\n"
              << "The program compiled and ran; it skipped the data-driven measurement\n"
              << "simulation because the inputs are absent.\n";
    return 0;
  }

  return 0;
}

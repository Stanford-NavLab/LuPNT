// Example 3: GNSS Interface
// ------------------------------------------------------------------------------
// C++ counterpart of python/examples/ex3_gnss_interface.ipynb.
//
// Exercises LuPNT's GNSS file interfaces the way the notebook does, but reports
// numerically instead of plotting:
//
//   Part 1  Satellite transmit-antenna gain patterns. Load a few GPS block /
//           Galileo / QZSS antenna tables with `Antenna` and print their peak
//           gain (the main lobe illuminates Earth; the sidelobes are what a
//           lunar receiver actually tracks).
//   Part 2  Antenna phase-center offsets (PCO) from an IGS20 ANTEX file, via the
//           native `AntexLoader`.
//   Part 3  Broadcast-vs-precise ephemeris comparison. Download the daily COD
//           MGEX final SP3 (`Sp3Loader`) and multi-GNSS BRDC broadcast
//           (`RinexNavLoader`) products for a short window, build a
//           `GnssConstellation` from the precise files, then for a handful of
//           common GPS/Galileo PRNs compare the broadcast orbit against the
//           precise orbit (shifted center-of-mass -> antenna phase center with
//           `AntexLoader::ApplyPcoCorrectionEcef`) in the radial/along-track/
//           cross-track (RTN) frame, and the broadcast vs precise satellite
//           clock (after removing the per-constellation time-reference bias).
//           An RTN + clock error summary table is printed.
//
// Data note: Parts 1-2 use bundled data (data/LuPNT_data/antenna,
// data/LuPNT_data/gnss/igs20.atx). Part 3 downloads SP3/BRDC from NASA CDDIS,
// which needs an Earthdata Login (~/.netrc or EARTHDATA_USERNAME/PASSWORD; see
// docs/pages/sp3_download). Absent the data/credentials the program prints an
// informative note and exits cleanly, so it always compiles and runs.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "lupnt/lupnt.h"

using namespace lupnt;

namespace {

  // ---- Small statistics helpers over std::vector<double> --------------------
  double Mean(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    double s = 0.0;
    for (double x : v) s += x;
    return s / static_cast<double>(v.size());
  }
  double Std(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    const double m = Mean(v);
    double s = 0.0;
    for (double x : v) s += (x - m) * (x - m);
    return std::sqrt(s / static_cast<double>(v.size()));
  }
  double Rms(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    double s = 0.0;
    for (double x : v) s += x * x;
    return std::sqrt(s / static_cast<double>(v.size()));
  }
  double Percentile95Abs(std::vector<double> v) {
    if (v.empty()) return 0.0;
    for (double& x : v) x = std::abs(x);
    std::sort(v.begin(), v.end());
    const double idx = 0.95 * static_cast<double>(v.size() - 1);
    const size_t lo = static_cast<size_t>(std::floor(idx));
    const size_t hi = static_cast<size_t>(std::ceil(idx));
    if (lo == hi) return v[lo];
    const double frac = idx - static_cast<double>(lo);
    return v[lo] * (1.0 - frac) + v[hi] * frac;
  }
  double Median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
  }

  // Decompose an ECI position difference into radial/along-track/cross-track.
  void RtnDecompose(const Vec3d& r_ref, const Vec3d& v_ref, const Vec3d& dr, double& R, double& T,
                    double& N) {
    const Vec3d r_hat = r_ref.normalized();
    const Vec3d h = r_ref.cross(v_ref);
    const Vec3d n_hat = h.normalized();
    const Vec3d t_hat = n_hat.cross(r_hat);
    R = dr.dot(r_hat);
    T = dr.dot(t_hat);
    N = dr.dot(n_hat);
  }

  struct SatSeries {
    std::string tag;
    std::vector<double> R, T, N;  // RTN orbit differences [m]
    std::vector<double> dclk_ns;  // broadcast - precise clock difference [ns]
  };

}  // namespace

int main() {
  std::cout << std::fixed;

  // ==== Part 1: transmit-antenna gain patterns ==============================
  std::cout << "Part 1: satellite transmit-antenna gain patterns (peak gain)\n";
  try {
    const std::vector<std::pair<std::string, std::string>> antennas
        = {{"GPS IIA (L1)", "Block-IIA_ACE"},     {"GPS IIR (L1)", "Block-IIR_ACE"},
           {"GPS IIR-M (L1)", "Block-IIR-M_ACE"}, {"GPS IIF (L1)", "Block-IIF_ACE"},
           {"Galileo E1", "Galileo_E1"},          {"Galileo E5a", "Galileo_E5a"},
           {"QZSS PRN02 L1", "QZSS_02_L1"},       {"QZSS PRN03 L1", "QZSS_03_L1"}};
    std::cout << "  " << std::left << std::setw(16) << "Signal" << std::right << std::setw(14)
              << "peak gain[dB]" << "\n";
    std::cout << "  " << std::string(30, '-') << "\n";
    for (const auto& [label, name] : antennas) {
      Antenna ant(name);
      MatXd gain = ant.GetGainMatrix();
      double peak = -1e300;
      for (int i = 0; i < gain.rows(); ++i)
        for (int j = 0; j < gain.cols(); ++j)
          if (std::isfinite(gain(i, j))) peak = std::max(peak, gain(i, j));
      std::cout << "  " << std::left << std::setw(16) << label << std::right << std::setw(14)
                << std::setprecision(1) << peak << "\n";
    }
  } catch (const std::exception& e) {
    std::cout << "  [skipped] antenna patterns unavailable: " << e.what() << "\n";
  }

  // ==== Part 2 + 3: PCO + broadcast-vs-precise comparison ===================
  // A single try/catch guards the data-dependent flow: bundled ANTEX plus the
  // CDDIS SP3/BRDC downloads. Any missing file / failed download exits cleanly.
  try {
    // ---- Comparison window (mirrors the notebook) --------------------------
    const int start_year = 2026, start_month = 1, start_day = 1;
    const int sim_days = 3;

    auto day_epoch_tai = [](int y, int m, int d) {
      const Real t_tdb = GregorianToTime(y, m, d, 12, 0, 0);
      return ConvertTime(t_tdb, Time::TDB, Time::TAI);
    };

    // Mid-day TAI epoch of the window's start date; other days are +/- whole
    // days from here (correct across month/year boundaries, unlike day + d).
    const double base_day_tai = day_epoch_tai(start_year, start_month, start_day).val();

    // ---- Part 2: ANTEX phase-center offsets -------------------------------
    const std::filesystem::path antex_path = GetDataPath() / "gnss" / "igs20.atx";
    AntexLoader antex(antex_path);
    const Real t_pco_tai = base_day_tai;

    std::cout << "\nPart 2: GPS L1 antenna phase-center offsets (IGS20 ANTEX)\n";
    std::cout << "  " << std::setw(5) << "PRN" << std::setw(10) << "N[mm]" << std::setw(10)
              << "E[mm]" << std::setw(10) << "U[mm]" << "\n";
    for (int prn = 1; prn <= 5; ++prn) {
      if (!antex.HasPco(GnssConst::GPS, prn, GnssFreq::L1, t_pco_tai)) continue;
      const Vec3d pco = antex.GetPco(GnssConst::GPS, prn, GnssFreq::L1, t_pco_tai);
      std::cout << "  G" << std::setfill('0') << std::setw(2) << prn << std::setfill(' ')
                << std::setw(10) << std::setprecision(1) << pco(0) * 1e3 << std::setw(10)
                << pco(1) * 1e3 << std::setw(10) << pco(2) * 1e3 << "\n";
    }

    // ---- Part 3: download the SP3 + BRDC products --------------------------
    std::cout << "\nPart 3: downloading SP3 + BRDC from CDDIS (" << sim_days
              << "-day window starting " << start_year << "-" << std::setfill('0') << std::setw(2)
              << start_month << "-" << std::setw(2) << start_day << std::setfill(' ') << ") ...\n";
    std::vector<std::filesystem::path> sp3_paths, brdc_paths;
    // One daily product per day, plus a +/-1 day margin for the interpolators.
    for (int d = -1; d <= sim_days; ++d) {
      const Real te = base_day_tai + d * 86400.0;
      sp3_paths.push_back(Sp3Loader::DownloadFileForEpoch(te, Time::TAI));
      brdc_paths.push_back(RinexNavLoader::DownloadFileForEpoch(te, Time::TAI));
    }

    Sp3Loader sp3(sp3_paths);
    RinexNavLoader brdc(brdc_paths);

    // ---- Build a GnssConstellation from the precise files (interface demo) -
    // The query grid: 15-min cadence over the window (mirrors the notebook).
    const Real t_start_tai = ConvertTime(
        GregorianToTime(start_year, start_month, start_day, 0, 0, 0), Time::TDB, Time::TAI);
    const double dt_q = 15.0 * 60.0;
    const int n_q = sim_days * 24 * 4 + 1;
    VecXd t_q_tai(n_q);
    for (int k = 0; k < n_q; ++k) t_q_tai(k) = t_start_tai.val() + k * dt_q;

    GnssConstellation gps_const(GnssConst::GPS);
    gps_const.SetupSatelliteStatesFromFiles(sp3_paths, antex_path, t_q_tai, GnssFreq::L1);
    std::cout << "  Built GnssConstellation (GPS): " << gps_const.GetNumSatellites()
              << " satellites; SP3 file(s) hold " << sp3.GetSatellites().size() << " satellites\n";

    // ---- Select up to 5 common GPS + Galileo PRNs --------------------------
    auto common_prns = [&](char sys, GnssConst gc, GnssFreq freq) {
      std::vector<int> out;
      for (int prn = 1; prn <= 40 && out.size() < 5; ++prn) {
        char id[8];
        std::snprintf(id, sizeof(id), "%c%02d", sys, prn);
        if (sp3.HasSatellite(id) && brdc.HasSatellite(id) && antex.HasPco(gc, prn, freq, t_pco_tai))
          out.push_back(prn);
      }
      return out;
    };
    const std::vector<int> gps_prns = common_prns('G', GnssConst::GPS, GnssFreq::L1);
    const std::vector<int> gal_prns = common_prns('E', GnssConst::GALILEO, GnssFreq::E1);

    // ---- Per-PRN RTN + clock differences -----------------------------------
    auto compare = [&](char sys, const std::vector<int>& prns, GnssConst gc, GnssFreq freq) {
      std::vector<SatSeries> out;
      for (int prn : prns) {
        char id[8];
        std::snprintf(id, sizeof(id), "%c%02d", sys, prn);
        SatSeries s;
        s.tag = id;
        for (int k = 0; k < n_q; ++k) {
          const Real t_tai = t_q_tai(k);
          const Real t_tdb = ConvertTime(t_tai, Time::TAI, Time::TDB);
          Vec6 rv_sp3_ecef, rv_brdc_ecef;
          Real clk_sp3, clk_brdc;
          try {
            sp3.GetPosVelClock(id, t_tai, rv_sp3_ecef, clk_sp3);
            brdc.GetPosVelClock(id, t_tai, rv_brdc_ecef, clk_brdc);
          } catch (const std::exception&) {
            continue;  // epoch outside coverage
          }
          // Shift SP3 center-of-mass -> antenna phase center (ECEF), keep vel.
          const Vec3d r_com_ecef = rv_sp3_ecef.head(3).cast<double>();
          const Vec3d pco = antex.GetPco(gc, prn, freq, t_tai);
          const Vec3d r_apc_ecef = AntexLoader::ApplyPcoCorrectionEcef(t_tai, r_com_ecef, pco);
          Vec6 rv_apc_ecef = rv_sp3_ecef;
          rv_apc_ecef.head(3) = r_apc_ecef.cast<Real>();
          // ECEF -> ECI for both.
          const Vec6 rv_apc_eci = ConvertFrame(t_tdb, rv_apc_ecef, Frame::ECEF, Frame::ECI);
          const Vec6 rv_brdc_eci = ConvertFrame(t_tdb, rv_brdc_ecef, Frame::ECEF, Frame::ECI);
          const Vec3d r_brdc = rv_brdc_eci.head(3).cast<double>();
          const Vec3d v_brdc = rv_brdc_eci.tail(3).cast<double>();
          const Vec3d dr = r_brdc - rv_apc_eci.head(3).cast<double>();
          double R, T, N;
          RtnDecompose(r_brdc, v_brdc, dr, R, T, N);
          s.R.push_back(R);
          s.T.push_back(T);
          s.N.push_back(N);
          // Store raw (brdc, sp3) clocks in ns; per-constellation bias removed below.
          s.dclk_ns.push_back((clk_brdc.val() - clk_sp3.val()) * 1e9);
        }
        out.push_back(std::move(s));
      }
      return out;
    };

    std::vector<SatSeries> gps = compare('G', gps_prns, GnssConst::GPS, GnssFreq::L1);
    std::vector<SatSeries> gal = compare('E', gal_prns, GnssConst::GALILEO, GnssFreq::E1);

    // Remove the systematic broadcast-vs-precise clock bias per constellation
    // (median over every satellite and epoch), per Montenbruck & Steigenberger.
    auto debias = [](std::vector<SatSeries>& sats) {
      std::vector<double> pooled;
      for (const SatSeries& s : sats)
        for (double x : s.dclk_ns) pooled.push_back(x);
      const double bias = Median(pooled);
      for (SatSeries& s : sats)
        for (double& x : s.dclk_ns) x -= bias;
      return bias;
    };
    const double gps_bias = debias(gps);
    const double gal_bias = debias(gal);

    // ---- Print the summary tables ------------------------------------------
    auto print_table = [](const std::string& title, const std::vector<SatSeries>& sats,
                          double clock_bias_ns) {
      std::cout << "\n"
                << title << "  (per-constellation clock bias removed: " << std::showpos
                << std::setprecision(3) << clock_bias_ns << " ns)" << std::noshowpos << "\n";
      std::cout << "  " << std::setw(5) << "PRN" << std::setw(11) << "R std[m]" << std::setw(11)
                << "T std[m]" << std::setw(11) << "N std[m]" << std::setw(13) << "clk std[ns]"
                << "\n";
      std::cout << "  " << std::string(49, '-') << "\n";
      std::vector<double> pool_R, pool_T, pool_N, pool_c;
      for (const SatSeries& s : sats) {
        std::cout << "  " << std::setw(5) << s.tag << std::setw(11) << std::setprecision(2)
                  << Std(s.R) << std::setw(11) << Std(s.T) << std::setw(11) << Std(s.N)
                  << std::setw(13) << Std(s.dclk_ns) << "\n";
        pool_R.insert(pool_R.end(), s.R.begin(), s.R.end());
        pool_T.insert(pool_T.end(), s.T.begin(), s.T.end());
        pool_N.insert(pool_N.end(), s.N.begin(), s.N.end());
        pool_c.insert(pool_c.end(), s.dclk_ns.begin(), s.dclk_ns.end());
      }
      std::cout << "  pooled (" << pool_R.size() << " samples): "
                << "RMS  R=" << std::setprecision(2) << Rms(pool_R) << "  T=" << Rms(pool_T)
                << "  N=" << Rms(pool_N) << " m   clk=" << Rms(pool_c) << " ns\n";
      std::cout << "                     95th R=" << Percentile95Abs(pool_R)
                << "  T=" << Percentile95Abs(pool_T) << "  N=" << Percentile95Abs(pool_N)
                << " m   clk=" << Percentile95Abs(pool_c) << " ns\n";
    };

    std::cout << "\nBRDC - (SP3 + PCO) orbit & clock differences (" << n_q << " epochs, 15-min):";
    print_table("GPS", gps, gps_bias);
    print_table("Galileo", gal, gal_bias);

  } catch (const std::exception& e) {
    std::cout << "\n[Part 3 skipped] SP3/BRDC/ANTEX data unavailable: " << e.what() << "\n"
              << "This example needs date-specific SP3/BRDC products from NASA CDDIS\n"
              << "(Earthdata Login) and the bundled IGS20 ANTEX file. See\n"
              << "docs/pages/sp3_download. The program compiled and ran; it skipped the\n"
              << "data-driven comparison because the inputs are absent.\n";
    return 0;
  }

  return 0;
}

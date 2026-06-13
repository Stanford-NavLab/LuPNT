/**
 * @file rinex_nav_loader.cc
 * @author Stanford NAV LAB
 * @brief RINEX navigation (broadcast ephemeris / BRDC) file loader
 * @version 0.1
 * @date 2025-06-07
 *
 * @copyright Copyright (c) 2025
 */
#include "lupnt/interfaces/rinex_nav_loader.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "lupnt/conversions/time_conversions.h"
#include "lupnt/core/error.h"
#include "lupnt/core/file.h"

namespace lupnt {

  namespace {
    constexpr double kGmEarth = 3.986005e14;  // [m^3/s^2] WGS-84 value used in broadcast nav
    constexpr double kOmegaDotEarth = 7.2921151467e-5;  // [rad/s] Earth rotation rate
    constexpr double kSecWeek = 604800.0;
    constexpr double kRelCorrFactor = -4.442807633e-10;  // [s / sqrt(m)] relativistic correction

    std::string Trim(const std::string& s) {
      size_t b = s.find_first_not_of(" \t\r\n");
      if (b == std::string::npos) return "";
      size_t e = s.find_last_not_of(" \t\r\n");
      return s.substr(b, e - b + 1);
    }

    std::string SatId(char sys, int prn) {
      std::ostringstream ss;
      ss << sys << std::setfill('0') << std::setw(2) << prn;
      return ss.str();
    }

    /// @brief TAI seconds at the GPS time epoch (1980-01-06 00:00:00 GPS).
    double GpsEpochTai() {
      static const double t0
          = ConvertTime(GregorianToTime(1980, 1, 6, 0, 0, 0.0), Time::GPS, Time::TAI).val();
      return t0;
    }

    /// @brief Convert TAI seconds to (GPS week, seconds-into-week).
    void TaiToGpsWeeks(double t_tai, double& gps_week, double& sec_week) {
      double delta = t_tai - GpsEpochTai();
      double delta_days = std::floor(delta / 86400.0);
      gps_week = std::floor(delta_days / 7.0);
      sec_week = delta - gps_week * 7.0 * 86400.0;
    }
  }  // namespace

  // Constructors ***************************************************************

  RinexNavLoader::RinexNavLoader(const std::filesystem::path& filepath) { LoadFile(filepath); }

  RinexNavLoader::RinexNavLoader(const std::vector<std::filesystem::path>& filepaths) {
    for (const auto& filepath : filepaths) LoadFile(filepath);
  }

  // Loading ********************************************************************

  void RinexNavLoader::LoadFile(const std::filesystem::path& filepath) {
    ParseFile(filepath);

    sats_.clear();
    for (const auto& [sat, _] : nav_) sats_.push_back(sat);
    std::sort(sats_.begin(), sats_.end());
  }

  void RinexNavLoader::ParseFile(const std::filesystem::path& filepath) {
    LUPNT_CHECK(std::filesystem::exists(filepath), "RINEX nav file not found: " + filepath.string(),
                "RinexNavLoader");

    std::ifstream file = OpenFile<std::ifstream>(filepath);

    bool in_header = true;
    std::string line;

    while (std::getline(file, line)) {
      if (in_header) {
        if (line.find("END OF HEADER") != std::string::npos) in_header = false;
        continue;
      }
      if (line.empty()) continue;

      char c0 = line[0];
      if (c0 != 'G' && c0 != 'E' && c0 != 'C' && c0 != 'J' && c0 != 'R') continue;

      std::istringstream hdr(line);
      std::string sat_token;
      int year, month, day, hh, mm;
      if (!(hdr >> sat_token >> year >> month >> day >> hh >> mm)) continue;
      if (sat_token.size() < 2) continue;

      char sys = sat_token[0];
      int prn = 0;
      try {
        prn = std::stoi(sat_token.substr(1));
      } catch (...) {
        continue;
      }

      if (sys == 'R') {
        // GLONASS uses tabulated orbital state vectors (not Keplerian
        // elements); not supported here -- skip the 3 continuation lines so
        // subsequent record parsing stays aligned.
        for (int k = 0; k < 3; k++) {
          if (!std::getline(file, line)) break;
        }
        continue;
      }
      if (sys != 'G' && sys != 'E' && sys != 'C' && sys != 'J') continue;
      if (line.size() < 80) continue;

      int ss;
      double af0, af1, af2;
      try {
        ss = std::stoi(line.substr(20, 2));
        af0 = std::stod(line.substr(23, 19));
        af1 = std::stod(line.substr(42, 19));
        af2 = std::stod(line.substr(61, 19));
      } catch (...) {
        continue;
      }

      double epoch_tai
          = ConvertTime(GregorianToTime(year, month, day, hh, mm, static_cast<double>(ss)),
                        Time::GPS, Time::TAI)
                .val();

      // Read the 7 continuation lines, each holding up to 4 fixed-width
      // values (column widths follow RINEX 3.x "D" exponent notation; the
      // first value on each line is 1 character wider iff the (stripped)
      // line starts with '-'). Mirrors `BRDCLoader.parse_brdc`'s value loop.
      std::vector<double> vals;
      vals.reserve(28);
      for (int k = 0; k < 7; k++) {
        if (!std::getline(file, line)) break;
        std::string l = Trim(line);
        if (l.empty()) continue;
        size_t start = 0;
        for (int kk = 0; kk < 4; kk++) {
          size_t lennum = (kk == 0) ? ((!l.empty() && l[0] == '-') ? 19 : 18) : 19;
          if (start >= l.size()) break;
          size_t len = std::min(lennum, l.size() - start);
          std::string val = Trim(l.substr(start, len));
          if (!val.empty()) {
            try {
              vals.push_back(std::stod(val));
            } catch (...) {
            }
          }
          start += lennum;
        }
      }

      // Positional mapping of `vals` to `ORBITS[sys][5:]` fields (this field
      // ordering is shared by GPS / Galileo / BeiDou / QZSS broadcast nav
      // messages in RINEX 3.x -- only the first element's name differs:
      // IODE / IODnav / AODE / IODE):
      //   0: IODE/IODnav/AODE   1: Crs        2: Delta_n   3: M0
      //   4: Cuc                5: ecc        6: Cus       7: sqrtA
      //   8: Toe                9: Cic       10: Omega0   11: Cis
      //  12: i0                13: Crc       14: omega    15: Omega_dot
      //  16: IDOT              17: Codes_L2/DataSrc/Spare1
      //  18: Week              ... (SV_accuracy/health/group-delay/T_trans/...)
      auto get = [&](size_t idx) -> double { return idx < vals.size() ? vals[idx] : 0.0; };

      NavMessage msg;
      msg.epoch_tai = epoch_tai;
      msg.af0 = af0;
      msg.af1 = af1;
      msg.af2 = af2;
      msg.crs = get(1);
      msg.delta_n = get(2);
      msg.m0 = get(3);
      msg.cuc = get(4);
      msg.ecc = get(5);
      msg.cus = get(6);
      msg.sqrt_a = get(7);
      msg.toe = get(8);
      msg.cic = get(9);
      msg.omega0 = get(10);
      msg.cis = get(11);
      msg.i0 = get(12);
      msg.crc = get(13);
      msg.omega = get(14);
      msg.omega_dot = get(15);
      msg.idot = get(16);
      msg.week = get(18);

      nav_[SatId(sys, prn)].push_back(msg);
    }
  }

  // Queries ********************************************************************

  bool RinexNavLoader::HasSatellite(const std::string& sat_id) const {
    return nav_.find(sat_id) != nav_.end();
  }

  int RinexNavLoader::FindClosestMessage(const std::string& sat_id, double t_tai) const {
    auto it = nav_.find(sat_id);
    LUPNT_CHECK(it != nav_.end() && !it->second.empty(),
                "Satellite '" + sat_id + "' not found in navigation data", "RinexNavLoader");
    const auto& msgs = it->second;

    int best = 0;
    double best_diff = std::abs(msgs[0].epoch_tai - t_tai);
    for (size_t i = 1; i < msgs.size(); i++) {
      double diff = std::abs(msgs[i].epoch_tai - t_tai);
      if (diff < best_diff) {
        best_diff = diff;
        best = static_cast<int>(i);
      }
    }
    return best;
  }

  void RinexNavLoader::GetPosVelClock(const std::string& sat_id, Real t_tai, Vec6& rv_ecef,
                                      Real& clock_corr_s) const {
    LUPNT_CHECK(sat_id.size() >= 1, "Invalid satellite identifier", "RinexNavLoader");
    char sys = sat_id[0];
    LUPNT_CHECK(sys == 'G' || sys == 'E' || sys == 'C' || sys == 'J',
                "Satellite system '" + std::string(1, sys)
                    + "' is not supported by RinexNavLoader (only G/E/C/J Keplerian-element "
                      "broadcast ephemerides are implemented)",
                "RinexNavLoader");

    double t = t_tai.val();
    int idx = FindClosestMessage(sat_id, t);
    const NavMessage& nav = nav_.at(sat_id)[idx];

    double gps_week, sec_week;
    TaiToGpsWeeks(t, gps_week, sec_week);

    // Time from ephemeris reference epoch (handle week rollover)
    double t_k = sec_week - nav.toe;
    if (t_k > 302400.0) {
      t_k -= kSecWeek;
    } else if (t_k < -302400.0) {
      t_k += kSecWeek;
    }

    // Mean anomaly & Kepler's equation (3 fixed-point iterations, matching the Python port)
    double n0 = std::sqrt(kGmEarth) / std::pow(nav.sqrt_a, 3);
    double n = n0 + nav.delta_n;
    double M = nav.m0 + n * t_k;

    double E = M;
    for (int iter = 0; iter < 3; iter++) {
      E = E + (M - (E - nav.ecc * std::sin(E))) / (1.0 - nav.ecc * std::cos(E));
    }
    double nu = 2.0 * std::atan(std::sqrt((1.0 + nav.ecc) / (1.0 - nav.ecc)) * std::tan(E / 2.0));

    // Clock correction: polynomial + relativistic terms.
    // (The Galileo-specific GST/GPST "GAGP" system-time term `t_corr_sys` is
    // intentionally omitted -- see the file-level note in `rinex_nav_loader.h`;
    // it does not affect the broadcast position/velocity computed below.)
    double t_corr_poly = nav.af0 + nav.af1 * t_k + nav.af2 * t_k * t_k;
    double t_corr_rel = kRelCorrFactor * nav.ecc * nav.sqrt_a * std::sin(E);
    double t_clk_corr = t_corr_poly + t_corr_rel;

    // Orbital-plane position with second-harmonic perturbation corrections
    double A = nav.sqrt_a * nav.sqrt_a;
    double phi = nu + nav.omega;
    double sin2phi = std::sin(2.0 * phi);
    double cos2phi = std::cos(2.0 * phi);

    double delta_u = nav.cus * sin2phi + nav.cuc * cos2phi;
    double delta_r = nav.crs * sin2phi + nav.crc * cos2phi;
    double delta_i = nav.cis * sin2phi + nav.cic * cos2phi;

    double u_k = phi + delta_u;
    double r_k = A * (1.0 - nav.ecc * std::cos(E)) + delta_r;
    double i_k = nav.i0 + delta_i + nav.idot * t_k;
    double Omega_k = nav.omega0 + (nav.omega_dot - kOmegaDotEarth) * t_k - kOmegaDotEarth * nav.toe;

    double cos_uk = std::cos(u_k), sin_uk = std::sin(u_k);
    double cos_Ok = std::cos(Omega_k), sin_Ok = std::sin(Omega_k);
    double cos_ik = std::cos(i_k), sin_ik = std::sin(i_k);

    double x_k_hat = r_k * cos_uk;
    double y_k_hat = r_k * sin_uk;

    double x_k = x_k_hat * cos_Ok - y_k_hat * cos_ik * sin_Ok;
    double y_k = x_k_hat * sin_Ok + y_k_hat * cos_ik * cos_Ok;
    double z_k = y_k_hat * sin_ik;

    // Velocities (analytic time derivatives)
    double Ekdot = n / (1.0 - nav.ecc * std::cos(E));
    double nudot = Ekdot * std::sqrt(1.0 - nav.ecc * nav.ecc) / (1.0 - nav.ecc * std::cos(E));
    double didot_dt = nav.idot + 2.0 * nudot * (nav.cis * cos2phi - nav.cic * sin2phi);
    double udot = nudot + 2.0 * nudot * (nav.cus * cos2phi - nav.cuc * sin2phi);
    double rdot
        = nav.ecc * A * Ekdot * std::sin(E) + 2.0 * nudot * (nav.crs * cos2phi - nav.crc * sin2phi);
    double Omega_k_dot = nav.omega_dot - kOmegaDotEarth;

    double xdot_hat = rdot * cos_uk - r_k * udot * sin_uk;
    double ydot_hat = rdot * sin_uk + r_k * udot * cos_uk;

    double xdot_k = -x_k_hat * Omega_k_dot * sin_Ok + xdot_hat * cos_Ok - ydot_hat * sin_Ok * cos_ik
                    - y_k_hat * (Omega_k_dot * cos_Ok * cos_ik - didot_dt * sin_Ok * sin_ik);
    double ydot_k = x_k_hat * Omega_k_dot * cos_Ok + xdot_hat * sin_Ok + ydot_hat * cos_Ok * cos_ik
                    - y_k_hat * (Omega_k_dot * sin_Ok * cos_ik + didot_dt * cos_Ok * sin_ik);
    double zdot_k = y_k_hat * didot_dt * cos_ik + ydot_hat * sin_ik;

    rv_ecef << x_k, y_k, z_k, xdot_k, ydot_k, zdot_k;
    clock_corr_s = t_clk_corr;
  }

  Vec6 RinexNavLoader::GetPosVel(const std::string& sat_id, Real t_tai) const {
    Vec6 rv_ecef;
    Real clock_corr_s;
    GetPosVelClock(sat_id, t_tai, rv_ecef, clock_corr_s);
    return rv_ecef;
  }

}  // namespace lupnt

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

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <system_error>

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

    // --- BRDC download helpers (mirror Sp3Loader::DownloadFileForEpoch) --------

    std::string ShellQuote(const std::string& s) {
      std::string out = "'";
      for (char c : s) {
        if (c == '\'') {
          out += "'\\''";
        } else {
          out += c;
        }
      }
      out += "'";
      return out;
    }

    bool IsTruthyEnv(const char* value) {
      if (value == nullptr) return false;
      std::string s(value);
      std::transform(s.begin(), s.end(), s.begin(),
                     [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
      return !(s.empty() || s == "0" || s == "false" || s == "off" || s == "no");
    }

    bool IsBrdcDownloadDisabled() { return IsTruthyEnv(std::getenv("LUPNT_SKIP_BRDC_DOWNLOAD")); }

    bool RunShellCommand(const std::string& cmd) { return std::system(cmd.c_str()) == 0; }

    bool LooksLikeHtml(const std::filesystem::path& filepath) {
      std::ifstream file(filepath, std::ios::binary);
      if (!file.is_open()) return false;

      std::string head(256, '\0');
      file.read(head.data(), static_cast<std::streamsize>(head.size()));
      head.resize(static_cast<size_t>(file.gcount()));

      auto first = std::find_if_not(head.begin(), head.end(),
                                    [](unsigned char c) { return std::isspace(c) != 0; });
      std::string trimmed(first, head.end());
      std::transform(trimmed.begin(), trimmed.end(), trimmed.begin(),
                     [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
      return trimmed.rfind("<!doctype html", 0) == 0 || trimmed.rfind("<html", 0) == 0;
    }

    bool DownloadToFileEarthdata(const std::string& url, const std::filesystem::path& dest_path) {
      if (IsBrdcDownloadDisabled()) return false;

      std::filesystem::path tmp_path = dest_path;
      tmp_path += ".part";
      std::error_code ec;
      std::filesystem::remove(tmp_path, ec);
      std::filesystem::create_directories(dest_path.parent_path(), ec);

      const bool has_env_auth = std::getenv("EARTHDATA_USERNAME") != nullptr
                                && std::getenv("EARTHDATA_PASSWORD") != nullptr;
      std::string auth_arg;
      if (has_env_auth) auth_arg = "-u \"$EARTHDATA_USERNAME:$EARTHDATA_PASSWORD\" ";

      // See the identical comment in sp3_loader.cc's DownloadToFileEarthdata: NASA's URS OAuth
      // redirect chain hands off session state via a cookie, so without -c/-b persisting it
      // across redirects, the final hop loops back to an unauthenticated request and restarts
      // the OAuth dance until curl's 50-redirect limit is hit -- even with correct credentials.
      const char* home = std::getenv("HOME");
      std::string cookie_jar = home != nullptr
                                   ? std::string(home) + "/.urs_cookies"
                                   : (dest_path.parent_path() / ".urs_cookies").string();

      std::string cmd = fmt::format(
          "curl -fsSL --netrc-optional -c {0} -b {0} {1}--connect-timeout 10 --max-time 600 -o "
          "{2} {3}",
          ShellQuote(cookie_jar), auth_arg, ShellQuote(tmp_path.string()), ShellQuote(url));

      bool ok = RunShellCommand(cmd) && std::filesystem::exists(tmp_path)
                && std::filesystem::file_size(tmp_path) > 0;
      if (!ok) {
        std::filesystem::remove(tmp_path, ec);
        return false;
      }
      std::filesystem::rename(tmp_path, dest_path, ec);
      if (ec) {
        ec.clear();
        std::filesystem::copy_file(tmp_path, dest_path,
                                   std::filesystem::copy_options::overwrite_existing, ec);
        std::filesystem::remove(tmp_path, ec);
      }
      return !ec;
    }

    void GunzipToFile(const std::filesystem::path& gzip_path,
                      const std::filesystem::path& dest_path) {
      std::filesystem::path tmp_path = dest_path;
      tmp_path += ".part";
      std::error_code ec;
      std::filesystem::remove(tmp_path, ec);

      std::string cmd = fmt::format("gzip -dc {} > {}", ShellQuote(gzip_path.string()),
                                    ShellQuote(tmp_path.string()));
      LUPNT_CHECK(RunShellCommand(cmd) && std::filesystem::exists(tmp_path)
                      && std::filesystem::file_size(tmp_path) > 0,
                  "Failed to decompress BRDC gzip file: " + gzip_path.string(), "RinexNavLoader");

      std::filesystem::rename(tmp_path, dest_path, ec);
      if (ec) {
        ec.clear();
        std::filesystem::copy_file(tmp_path, dest_path,
                                   std::filesystem::copy_options::overwrite_existing, ec);
        std::filesystem::remove(tmp_path, ec);
      }
      LUPNT_CHECK(!ec, "Failed to write decompressed BRDC file: " + dest_path.string(),
                  "RinexNavLoader");
    }

    int DayOfYear(int year, int month, int day) {
      static constexpr std::array<int, 12> kMonthStartsCommon
          = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
      bool leap = (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
      return kMonthStartsCommon.at(static_cast<size_t>(month - 1)) + day
             + (leap && month > 2 ? 1 : 0);
    }

    struct BrdcProductInfo {
      int year = 0;
      int doy = 0;
      std::string filename;         // uncompressed ".rnx"
      std::string zipped_filename;  // ".rnx.gz"
      std::string url_primary;      // .../daily/YYYY/DDD/YYp/<gz>
      std::string url_fallback;     // .../daily/YYYY/brdc/<gz>
    };

    BrdcProductInfo BuildBrdcProductInfo(Real epoch, Time time_scale) {
      // Derive the calendar day the same way Sp3Loader does (via GPS time), so
      // SP3 and BRDC pick the same daily product for a given query epoch.
      const Real t_gps = ConvertTime(epoch, time_scale, Time::GPS);
      auto [year, month, day, hour, minute, second] = MjdToGregorian(TimeToMjd(t_gps));
      (void)hour;
      (void)minute;
      (void)second;

      BrdcProductInfo info;
      info.year = year;
      info.doy = DayOfYear(year, month, day);
      const std::string yy = fmt::format("{:02d}", info.year % 100);

      info.filename = fmt::format("BRDC00IGS_R_{:04d}{:03d}0000_01D_MN.rnx", info.year, info.doy);
      info.zipped_filename = info.filename + ".gz";
      const std::string base = "https://cddis.nasa.gov/archive/gnss/data/daily";
      info.url_primary = fmt::format("{}/{:04d}/{:03d}/{}p/{}", base, info.year, info.doy, yy,
                                     info.zipped_filename);
      info.url_fallback = fmt::format("{}/{:04d}/brdc/{}", base, info.year, info.zipped_filename);
      return info;
    }
  }  // namespace

  // Constructors ***************************************************************

  RinexNavLoader::RinexNavLoader(const std::filesystem::path& filepath) { LoadFile(filepath); }

  RinexNavLoader::RinexNavLoader(const std::vector<std::filesystem::path>& filepaths) {
    for (const auto& filepath : filepaths) LoadFile(filepath);
  }

  // Download helpers ***********************************************************

  std::string RinexNavLoader::FilenameForEpoch(Real epoch, Time time_scale) {
    return BuildBrdcProductInfo(epoch, time_scale).filename;
  }

  std::filesystem::path RinexNavLoader::DownloadFileForEpoch(
      Real epoch, Time time_scale, const std::filesystem::path& cache_dir) {
    const BrdcProductInfo info = BuildBrdcProductInfo(epoch, time_scale);
    const std::filesystem::path brdc_dir
        = cache_dir.empty() ? (GetOutputDir("gnss_files") / "brdc") : cache_dir;
    const std::filesystem::path rnx_path = brdc_dir / info.filename;

    if (std::filesystem::exists(rnx_path)) return rnx_path;

    const std::filesystem::path gzip_path = brdc_dir / info.zipped_filename;
    if (!std::filesystem::exists(gzip_path)) {
      // CDDIS hosts BRDC under two layouts; try the dated one first, then the
      // per-year `brdc/` alias (mirrors BRDCLoader.load_brdc).
      bool ok = DownloadToFileEarthdata(info.url_primary, gzip_path);
      if (ok && LooksLikeHtml(gzip_path)) {
        std::error_code ec;
        std::filesystem::remove(gzip_path, ec);
        ok = false;
      }
      if (!ok) {
        std::error_code ec;
        std::filesystem::remove(gzip_path, ec);
        ok = DownloadToFileEarthdata(info.url_fallback, gzip_path);
      }
      LUPNT_CHECK(ok,
                  "Failed to download BRDC file from CDDIS. Configure Earthdata credentials with "
                  "~/.netrc or EARTHDATA_USERNAME/EARTHDATA_PASSWORD, then retry. URLs: "
                      + info.url_primary + " , " + info.url_fallback,
                  "RinexNavLoader");
    }

    if (LooksLikeHtml(gzip_path)) {
      std::error_code ec;
      std::filesystem::remove(gzip_path, ec);
      LUPNT_CHECK(false,
                  "CDDIS returned an Earthdata Login HTML page instead of the requested BRDC file. "
                  "Configure Earthdata credentials with ~/.netrc or EARTHDATA_USERNAME/"
                  "EARTHDATA_PASSWORD, then retry. URL: "
                      + info.url_primary,
                  "RinexNavLoader");
    }

    GunzipToFile(gzip_path, rnx_path);
    return rnx_path;
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

  void RinexNavLoader::LoadYumaFile(const std::filesystem::path& filepath) {
    LUPNT_CHECK(std::filesystem::exists(filepath),
                "YUMA almanac file not found: " + filepath.string(), "RinexNavLoader");
    std::ifstream file = OpenFile<std::ifstream>(filepath);

    // YUMA blocks are "key: value" lines separated by a header line ("******** Week ...").
    // Accumulate fields per block (keyed by a lowercase substring of the label) and flush on
    // the next header / EOF. Only GPS satellites are produced.
    std::map<std::string, double> f;
    auto flush = [&]() {
      if (f.empty()) return;
      auto has = [&](const std::string& k) { return f.find(k) != f.end(); };
      if (has("id") && has("sqrt")) {
        int prn = static_cast<int>(std::llround(f["id"]));
        NavMessage msg;
        msg.ecc = f.count("eccentricity") ? f["eccentricity"] : 0.0;
        msg.toe = f.count("applicability") ? f["applicability"] : 0.0;
        msg.i0 = f.count("inclination") ? f["inclination"] : 0.0;
        msg.omega_dot = f.count("ascen(r/s)") ? f["ascen(r/s)"] : 0.0;
        msg.sqrt_a = f["sqrt"];
        msg.omega0 = f.count("week(rad)") ? f["week(rad)"] : 0.0;
        msg.omega = f.count("perigee") ? f["perigee"] : 0.0;
        msg.m0 = f.count("anom") ? f["anom"] : 0.0;
        msg.af0 = f.count("af0") ? f["af0"] : 0.0;
        msg.af1 = f.count("af1") ? f["af1"] : 0.0;
        msg.week = f.count("week") ? f["week"] : 0.0;
        // Harmonic corrections, delta_n and idot are not carried by YUMA.
        msg.epoch_tai = GpsEpochTai() + msg.week * kSecWeek + msg.toe;
        if (prn > 0 && msg.sqrt_a > 0.0) nav_[SatId('G', prn)].push_back(msg);
      }
      f.clear();
    };

    std::string line;
    while (std::getline(file, line)) {
      const std::string trimmed = Trim(line);
      if (trimmed.empty()) continue;
      if (trimmed.rfind("****", 0) == 0) {  // block header
        flush();
        continue;
      }
      const size_t colon = trimmed.find(':');
      if (colon == std::string::npos) continue;
      std::string label = trimmed.substr(0, colon);
      const std::string value = Trim(trimmed.substr(colon + 1));
      std::transform(label.begin(), label.end(), label.begin(),
                     [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
      double v = 0.0;
      try {
        v = std::stod(value);
      } catch (...) {
        continue;
      }
      // Store under a stable substring key that survives the exact-label wording, e.g.
      // "Rate of Right Ascen(r/s)" -> match "ascen(r/s)"; "SQRT(A)  (m 1/2)" -> "sqrt".
      auto keyed = [&](const std::string& needle, const std::string& key) {
        if (label.find(needle) != std::string::npos) f[key] = v;
      };
      keyed("id", "id");
      keyed("eccentricity", "eccentricity");
      keyed("applicability", "applicability");
      keyed("inclination", "inclination");
      keyed("ascen(r/s)", "ascen(r/s)");
      keyed("sqrt", "sqrt");
      keyed("week(rad)", "week(rad)");
      keyed("perigee", "perigee");
      keyed("anom", "anom");
      keyed("af0", "af0");
      keyed("af1", "af1");
      // Plain "week:" (no "(rad)") is the GPS week number; guard against the RAAN label.
      if (label.find("week") != std::string::npos && label.find("rad") == std::string::npos)
        f["week"] = v;
    }
    flush();

    sats_.clear();
    for (const auto& [sat, _] : nav_) sats_.push_back(sat);
    std::sort(sats_.begin(), sats_.end());
  }

  // Queries ********************************************************************

  bool RinexNavLoader::HasSatellite(const std::string& sat_id) const {
    return nav_.find(sat_id) != nav_.end();
  }

  double RinexNavLoader::GetLatestEpochTai(const std::string& sat_id) const {
    auto it = nav_.find(sat_id);
    LUPNT_CHECK(it != nav_.end() && !it->second.empty(),
                "Satellite '" + sat_id + "' not found in navigation data", "RinexNavLoader");
    double latest = it->second.front().epoch_tai;
    for (const auto& msg : it->second) latest = std::max(latest, msg.epoch_tai);
    return latest;
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

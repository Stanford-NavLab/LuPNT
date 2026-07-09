/**
 * @file sp3_loader.cc
 * @author Stanford NAV LAB
 * @brief SP3 (precise GNSS ephemeris) file loader
 * @version 0.1
 * @date 2025-06-07
 *
 * @copyright Copyright (c) 2025
 */
#include "lupnt/interfaces/sp3_loader.h"

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <numeric>
#include <sstream>
#include <system_error>

#include "lupnt/conversions/time_conversions.h"
#include "lupnt/core/error.h"
#include "lupnt/core/file.h"
#include "lupnt/numerics/interpolation.h"

namespace lupnt {

  namespace {
    /// @brief Speed of light [m/s] (matches `gnss_file_loader.SP3Loader`).
    constexpr double kC = 299792458.0;
    constexpr double kSecsWeek = 7.0 * SECS_DAY;
    /// @brief SP3 sentinel for unavailable clock data: 999999.999999 µs.
    /// Values above this threshold (in µs, before ×1e-6 conversion) are invalid.
    constexpr double kSp3ClockSentinelUs = 999990.0;

    struct Sp3ProductInfo {
      int gps_week = 0;
      int day_of_week = 0;
      int year = 0;
      int doy = 0;
      std::string filename;
      std::string zipped_filename;
      std::string url;
    };

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

    bool IsSp3DownloadDisabled() { return IsTruthyEnv(std::getenv("LUPNT_SKIP_SP3_DOWNLOAD")); }

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
      if (IsSp3DownloadDisabled()) return false;

      std::filesystem::path tmp_path = dest_path;
      tmp_path += ".part";
      std::error_code ec;
      std::filesystem::remove(tmp_path, ec);
      std::filesystem::create_directories(dest_path.parent_path(), ec);

      const bool has_env_auth = std::getenv("EARTHDATA_USERNAME") != nullptr
                                && std::getenv("EARTHDATA_PASSWORD") != nullptr;
      std::string auth_arg;
      if (has_env_auth) auth_arg = "-u \"$EARTHDATA_USERNAME:$EARTHDATA_PASSWORD\" ";

      std::string cmd = fmt::format(
          "curl -fsSL --netrc-optional {}--connect-timeout 10 --max-time 600 -o {} {}", auth_arg,
          ShellQuote(tmp_path.string()), ShellQuote(url));

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
                  "Failed to decompress SP3 gzip file: " + gzip_path.string(), "Sp3Loader");

      std::filesystem::rename(tmp_path, dest_path, ec);
      if (ec) {
        ec.clear();
        std::filesystem::copy_file(tmp_path, dest_path,
                                   std::filesystem::copy_options::overwrite_existing, ec);
        std::filesystem::remove(tmp_path, ec);
      }
      LUPNT_CHECK(!ec, "Failed to write decompressed SP3 file: " + dest_path.string(), "Sp3Loader");
    }

    int DayOfYear(int year, int month, int day) {
      static constexpr std::array<int, 12> kMonthStartsCommon
          = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
      bool leap = (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
      return kMonthStartsCommon.at(static_cast<size_t>(month - 1)) + day
             + (leap && month > 2 ? 1 : 0);
    }

    Sp3ProductInfo BuildSp3ProductInfo(Real epoch, Time time_scale) {
      const Real t_gps = ConvertTime(epoch, time_scale, Time::GPS);
      const Real gps_epoch = GregorianToTime(1980, 1, 6, 0, 0, 0.0);
      const double delta = t_gps.val() - gps_epoch.val();
      LUPNT_CHECK(delta >= 0.0, "SP3 download helper does not support epochs before GPS epoch",
                  "Sp3Loader");

      Sp3ProductInfo info;
      info.gps_week = static_cast<int>(std::floor(delta / kSecsWeek));
      info.day_of_week
          = static_cast<int>(std::floor((delta - info.gps_week * kSecsWeek) / SECS_DAY));

      LUPNT_CHECK(info.gps_week >= 1962,
                  "C++ SP3 download helper supports modern CDDIS COD/MGEX `.SP3.gz` products "
                  "(GPS week >= 1962) only",
                  "Sp3Loader");

      auto [year, month, day, hour, minute, second] = MjdToGregorian(TimeToMjd(t_gps));
      (void)hour;
      (void)minute;
      (void)second;
      info.year = year;
      info.doy = DayOfYear(year, month, day);

      info.filename
          = fmt::format("COD0MGXFIN_{:04d}{:03d}0000_01D_05M_ORB.SP3", info.year, info.doy);
      info.zipped_filename = info.filename + ".gz";
      info.url = fmt::format("https://cddis.nasa.gov/archive/gnss/products/{:04d}/{}",
                             info.gps_week, info.zipped_filename);
      return info;
    }

  }  // namespace

  // Constructors ***************************************************************

  Sp3Loader::Sp3Loader(const std::filesystem::path& filepath) { LoadFile(filepath); }

  Sp3Loader::Sp3Loader(const std::vector<std::filesystem::path>& filepaths) {
    for (const auto& filepath : filepaths) LoadFile(filepath);
  }

  std::string Sp3Loader::FilenameForEpoch(Real epoch, Time time_scale) {
    return BuildSp3ProductInfo(epoch, time_scale).filename;
  }

  std::string Sp3Loader::UrlForEpoch(Real epoch, Time time_scale) {
    return BuildSp3ProductInfo(epoch, time_scale).url;
  }

  std::filesystem::path Sp3Loader::DownloadFileForEpoch(Real epoch, Time time_scale,
                                                        const std::filesystem::path& cache_dir) {
    const Sp3ProductInfo info = BuildSp3ProductInfo(epoch, time_scale);
    const std::filesystem::path sp3_dir
        = cache_dir.empty() ? (GetOutputDir("gnss_files") / "sp3") : cache_dir;
    const std::filesystem::path sp3_path = sp3_dir / info.filename;

    if (std::filesystem::exists(sp3_path)) return sp3_path;

    const std::filesystem::path gzip_path = sp3_dir / info.zipped_filename;
    if (!std::filesystem::exists(gzip_path)) {
      LUPNT_CHECK(DownloadToFileEarthdata(info.url, gzip_path),
                  "Failed to download SP3 file from CDDIS. Configure Earthdata credentials with "
                  "~/.netrc or EARTHDATA_USERNAME/EARTHDATA_PASSWORD, then retry. URL: "
                      + info.url,
                  "Sp3Loader");
    }

    if (LooksLikeHtml(gzip_path)) {
      std::error_code ec;
      std::filesystem::remove(gzip_path, ec);
      LUPNT_CHECK(false,
                  "CDDIS returned an Earthdata Login HTML page instead of the requested SP3 file. "
                  "Configure Earthdata credentials with ~/.netrc or EARTHDATA_USERNAME/"
                  "EARTHDATA_PASSWORD, then retry. URL: "
                      + info.url,
                  "Sp3Loader");
    }

    GunzipToFile(gzip_path, sp3_path);
    std::error_code ec;
    std::filesystem::remove(gzip_path, ec);
    return sp3_path;
  }

  // Loading ********************************************************************

  void Sp3Loader::LoadFile(const std::filesystem::path& filepath) {
    ParseFile(filepath);

    // (Re-)build sorted, de-duplicated per-satellite epoch/position tables.
    // Mirrors the sort + `np.unique` de-duplication in `SP3Loader.get_posvelclock`.
    sats_.clear();
    for (const auto& [sat, _] : epochs_tai_) sats_.push_back(sat);
    std::sort(sats_.begin(), sats_.end());

    for (const auto& sat : sats_) {
      VecXd& epochs = epochs_tai_.at(sat);
      MatXd& pc = pos_clock_.at(sat);

      std::vector<int> order(epochs.size());
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(),
                [&epochs](int a, int b) { return epochs(a) < epochs(b); });

      VecXd epochs_sorted(epochs.size());
      MatXd pc_sorted(pc.rows(), pc.cols());
      for (int i = 0; i < static_cast<int>(order.size()); i++) {
        epochs_sorted(i) = epochs(order[i]);
        pc_sorted.row(i) = pc.row(order[i]);
      }

      // De-duplicate consecutive equal epochs.
      // Policy: at day-file boundaries the same epoch appears in both files;
      // the earlier file typically has a sentinel clock (999999.999999 µs → NaN)
      // while the later file has the real value.  Prefer the record with a valid
      // (non-NaN) clock; if both have the same clock status, keep the later one
      // (from the newer file, which supersedes the older boundary epoch).
      std::vector<int> keep;
      keep.reserve(epochs_sorted.size());
      for (int i = 0; i < static_cast<int>(epochs_sorted.size()); i++) {
        if (i == 0 || epochs_sorted(i) != epochs_sorted(i - 1)) {
          keep.push_back(i);
        } else {
          bool prev_nan = std::isnan(pc_sorted(keep.back(), 3));
          bool curr_nan = std::isnan(pc_sorted(i, 3));
          if (!prev_nan && curr_nan) {
            // Previous has valid clock; current is sentinel → keep previous.
          } else {
            // Current has valid clock (and previous is sentinel), or both have
            // the same status → prefer the later record.
            keep.back() = i;
          }
        }
      }

      VecXd epochs_unique(keep.size());
      MatXd pc_unique(keep.size(), pc.cols());
      for (size_t i = 0; i < keep.size(); i++) {
        epochs_unique(i) = epochs_sorted(keep[i]);
        pc_unique.row(i) = pc_sorted.row(keep[i]);
      }

      epochs = epochs_unique;
      pc = pc_unique;
    }

    RebuildChebyshevModels();
  }

  void Sp3Loader::ParseFile(const std::filesystem::path& filepath) {
    LUPNT_CHECK(std::filesystem::exists(filepath), "SP3 file not found: " + filepath.string(),
                "Sp3Loader");

    std::ifstream file = OpenFile<std::ifstream>(filepath);

    // Accumulate raw samples per-satellite, then concatenate with any
    // previously-loaded data before sorting (done in `LoadFile`).
    std::map<std::string, std::vector<double>> new_epochs;
    std::map<std::string, std::vector<std::array<double, 4>>> new_pos_clock;

    std::string line;
    Real current_epoch = 0.0;
    bool have_epoch = false;

    while (std::getline(file, line)) {
      if (line.empty()) continue;

      if (line[0] == '*') {
        // Epoch line: "* yyyy mm dd hh mm ss.sssss"
        std::istringstream ss(line.substr(1));
        int year, month, day, hour, minute;
        double second;
        ss >> year >> month >> day >> hour >> minute >> second;

        current_epoch = ConvertTime(GregorianToTime(year, month, day, hour, minute, second),
                                    Time::GPS, Time::TAI);
        have_epoch = true;

      } else if (line[0] == 'P' && have_epoch) {
        // Position line: "P{sv}{x:14.6f}{y:14.6f}{z:14.6f}{clock:14.6f}" (fixed columns)
        if (line.size() < 60) continue;
        std::string sv = line.substr(1, 3);
        sv.erase(sv.find_last_not_of(" \t") + 1);
        sv.erase(0, sv.find_first_not_of(" \t"));
        if (sv.empty()) continue;

        double x = std::stod(line.substr(4, 14)) * 1000.0;   // km -> m
        double y = std::stod(line.substr(18, 14)) * 1000.0;  // km -> m
        double z = std::stod(line.substr(32, 14)) * 1000.0;  // km -> m
        double clock_us = std::stod(line.substr(46, 14));
        // Store sentinel (999999.999999 µs) as NaN; valid biases are << 1 s.
        double clock = (clock_us > kSp3ClockSentinelUs) ? std::numeric_limits<double>::quiet_NaN()
                                                        : clock_us * 1e-6;

        new_epochs[sv].push_back(current_epoch.val());
        new_pos_clock[sv].push_back({x, y, z, clock});
      }
    }

    for (auto& [sat, ep_vec] : new_epochs) {
      const auto& pc_vec = new_pos_clock.at(sat);
      int n_new = static_cast<int>(ep_vec.size());

      auto it_ep = epochs_tai_.find(sat);
      if (it_ep == epochs_tai_.end()) {
        VecXd epochs(n_new);
        MatXd pc(n_new, 4);
        for (int i = 0; i < n_new; i++) {
          epochs(i) = ep_vec[i];
          for (int k = 0; k < 4; k++) pc(i, k) = pc_vec[i][k];
        }
        epochs_tai_[sat] = epochs;
        pos_clock_[sat] = pc;
      } else {
        VecXd& epochs = it_ep->second;
        MatXd& pc = pos_clock_.at(sat);
        int n_old = static_cast<int>(epochs.size());

        VecXd epochs_cat(n_old + n_new);
        MatXd pc_cat(n_old + n_new, 4);
        epochs_cat.head(n_old) = epochs;
        pc_cat.topRows(n_old) = pc;
        for (int i = 0; i < n_new; i++) {
          epochs_cat(n_old + i) = ep_vec[i];
          for (int k = 0; k < 4; k++) pc_cat(n_old + i, k) = pc_vec[i][k];
        }
        epochs = epochs_cat;
        pc = pc_cat;
      }
    }
  }

  void Sp3Loader::RebuildChebyshevModels() {
    pos_cheby_.clear();
    for (const auto& sat : sats_) {
      const VecXd& epochs = epochs_tai_.at(sat);
      const MatXd& pc = pos_clock_.at(sat);
      if (epochs.size() < 2) continue;

      double min_step = epochs(epochs.size() - 1) - epochs(0);
      for (int i = 1; i < epochs.size(); ++i) {
        if (epochs(i) > epochs(i - 1)) min_step = std::min(min_step, epochs(i) - epochs(i - 1));
      }
      const double segment_length = std::max(1.0, 8.0 * min_step);
      const int num_coeffs = std::min(12, std::max(2, static_cast<int>(epochs.size())));

      // Fit Chebyshev only for position (x, y, z); clock uses linear interpolation
      // at query time to avoid ringing / Gibbs-like jumps at day boundaries. Sample the
      // tabulated positions with a local high-order Lagrange interpolant (the standard SP3
      // interpolation) rather than a linear one: a GPS/Galileo arc is strongly curved over the
      // 5-min product interval, so linear sampling incurs kilometre-level chord error that the
      // Chebyshev fit would faithfully reproduce.
      const int lagrange_order = std::min<int>(10, std::max<int>(2, static_cast<int>(epochs.size())));
      pos_cheby_[sat] = FitChebyshevModel(
          [&epochs, &pc, lagrange_order](double t) {
            LagrangeInterpolator interp(epochs, t, lagrange_order);
            VecXd sample(3);
            for (int k = 0; k < 3; ++k) {
              VecXd col = pc.col(k);
              sample(k) = interp.Interpolate(col);
            }
            return sample;
          },
          epochs(0), epochs(epochs.size() - 1), 3, segment_length, num_coeffs);
    }
  }

  // Queries ********************************************************************

  bool Sp3Loader::HasSatellite(const std::string& sat_id) const {
    return epochs_tai_.find(sat_id) != epochs_tai_.end();
  }

  std::pair<double, double> Sp3Loader::GetTimeSpan(const std::string& sat_id) const {
    auto it = epochs_tai_.find(sat_id);
    LUPNT_CHECK(it != epochs_tai_.end(), "Satellite '" + sat_id + "' not found in SP3 ephemeris",
                "Sp3Loader");
    const VecXd& epochs = it->second;
    return {epochs(0), epochs(epochs.size() - 1)};
  }

  void Sp3Loader::GetPosVelClock(const std::string& sat_id, Real t_tai, Vec6& rv_ecef,
                                 Real& clock_bias_s) const {
    auto it_ep = epochs_tai_.find(sat_id);
    LUPNT_CHECK(it_ep != epochs_tai_.end(), "Satellite '" + sat_id + "' not found in SP3 ephemeris",
                "Sp3Loader");
    const VecXd& epochs = it_ep->second;

    LUPNT_CHECK(epochs.size() >= 2,
                "SP3 ephemeris for satellite '" + sat_id + "' has fewer than 2 samples",
                "Sp3Loader");

    double t = t_tai.val();
    LUPNT_CHECK(t >= epochs(0) && t <= epochs(epochs.size() - 1),
                "Requested epoch is outside the loaded SP3 ephemeris span for satellite '" + sat_id
                    + "' (no orbital-propagation fallback)",
                "Sp3Loader");

    auto it_model = pos_cheby_.find(sat_id);
    LUPNT_CHECK(
        it_model != pos_cheby_.end(),
        "Chebyshev SP3 interpolation model for satellite '" + sat_id + "' has not been built",
        "Sp3Loader");

    VecX pos;
    VecX pos_dot;
    LUPNT_CHECK(it_model->second.Eval(t_tai, &pos, &pos_dot),
                "Requested epoch is outside the Chebyshev SP3 interpolation span for satellite '"
                    + sat_id + "'",
                "Sp3Loader");

    rv_ecef.head(3) = pos.head(3);
    rv_ecef.tail(3) = pos_dot.head(3);

    // Clock: linear interpolation over raw SP3 samples, excluding sentinel
    // (NaN) epochs so that day-boundary 999999.999999 µs values do not corrupt
    // the interpolation.
    VecXd clk_col = pos_clock_.at(sat_id).col(3);
    std::vector<double> valid_ep, valid_clk;
    valid_ep.reserve(static_cast<size_t>(epochs.size()));
    valid_clk.reserve(static_cast<size_t>(epochs.size()));
    for (int i = 0; i < epochs.size(); ++i) {
      if (!std::isnan(clk_col(i))) {
        valid_ep.push_back(epochs(i));
        valid_clk.push_back(clk_col(i));
      }
    }
    if (valid_ep.size() >= 2) {
      VecXd vep = Eigen::Map<VecXd>(valid_ep.data(), static_cast<int>(valid_ep.size()));
      VecXd vclk = Eigen::Map<VecXd>(valid_clk.data(), static_cast<int>(valid_clk.size()));
      clock_bias_s = LinearInterp1d(vep, vclk, t);
    } else {
      clock_bias_s = std::numeric_limits<double>::quiet_NaN();
    }

    // Relativistic clock correction: t_corr_rel = -2/c^2 * dot(r, v)
    Real dot_rv = rv_ecef.head(3).dot(rv_ecef.tail(3));
    Real t_corr_rel = -2.0 / (kC * kC) * dot_rv;
    clock_bias_s = clock_bias_s + t_corr_rel;
  }

  Vec6 Sp3Loader::GetPosVel(const std::string& sat_id, Real t_tai) const {
    Vec6 rv_ecef;
    Real clock_bias_s;
    GetPosVelClock(sat_id, t_tai, rv_ecef, clock_bias_s);
    return rv_ecef;
  }

}  // namespace lupnt

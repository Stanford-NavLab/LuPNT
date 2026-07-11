#include "lupnt/simulations/lunar_gnss_odts/lunar_gnss_odts_simulation.h"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>

#include "lupnt/lupnt.h"

namespace lupnt {
  using namespace lupnt;

  namespace {
    template <typename T>
    T ReadYaml(const YAML::Node& node, const std::string& key, const T& fallback) {
      return node && node[key] ? node[key].as<T>() : fallback;
    }

    std::filesystem::path ResolvePath(const std::filesystem::path& base,
                                      const std::filesystem::path& path) {
      if (path.empty() || path.is_absolute()) return path;
      return base / path;
    }

    bool FileHasContent(const std::filesystem::path& path) {
      std::error_code ec;
      return std::filesystem::exists(path, ec) && std::filesystem::is_regular_file(path, ec)
             && std::filesystem::file_size(path, ec) > 0;
    }

    std::filesystem::path LinkMetadataPath(const std::filesystem::path& links_file) {
      return std::filesystem::path(links_file.string() + ".meta");
    }

    std::string PathListString(const std::vector<std::filesystem::path>& paths) {
      std::vector<std::string> normalized;
      normalized.reserve(paths.size());
      for (const auto& path : paths) normalized.push_back(path.lexically_normal().string());
      std::sort(normalized.begin(), normalized.end());
      std::ostringstream oss;
      for (const auto& path : normalized) oss << path << ";";
      return oss.str();
    }

    template <typename T>
    void AppendFingerprintField(std::ostringstream& oss, const std::string& key, const T& value) {
      oss << key << "=" << value << "\n";
    }

    template <typename T> void AppendFingerprintVector(std::ostringstream& oss,
                                                       const std::string& key,
                                                       const std::vector<T>& values) {
      oss << key << "=";
      for (const auto& value : values) oss << value << ",";
      oss << "\n";
    }

    std::string LinkCacheFingerprint(const LunarGnssODTSConfig& cfg) {
      std::ostringstream oss;
      oss << std::setprecision(17);
      AppendFingerprintField(oss, "version", 17);  // 17: high-order (Lagrange) SP3 interpolation
      AppendFingerprintField(oss, "seed", cfg.seed);
      AppendFingerprintField(oss, "duration_s", cfg.duration_s);
      AppendFingerprintField(oss, "dt_s", cfg.dt_s);
      AppendFingerprintField(oss, "ephemeris_dt_s", cfg.ephemeris_dt_s);
      AppendFingerprintField(oss, "start_epoch_utc", cfg.start_epoch_utc);
      AppendFingerprintField(oss, "receiver_rate_hz", cfg.receiver_app.rate_hz);
      AppendFingerprintField(oss, "receiver_a_m", cfg.receiver_a_m);
      AppendFingerprintField(oss, "receiver_ecc", cfg.receiver_ecc);
      AppendFingerprintField(oss, "receiver_inc_rad", cfg.receiver_inc_rad);
      AppendFingerprintField(oss, "receiver_raan_rad", cfg.receiver_raan_rad);
      AppendFingerprintField(oss, "receiver_argp_rad", cfg.receiver_argp_rad);
      AppendFingerprintField(oss, "receiver_mean_anomaly_rad", cfg.receiver_mean_anomaly_rad);
      AppendFingerprintField(oss, "clock_bias_s", cfg.clock_bias_s);
      AppendFingerprintField(oss, "clock_drift_sps", cfg.clock_drift_sps);
      AppendFingerprintField(oss, "sp3_files", PathListString(cfg.constellation.sp3_files));
      AppendFingerprintField(oss, "antex_file",
                             cfg.constellation.antex_file.lexically_normal().string());
      AppendFingerprintField(oss, "use_all_gps", cfg.constellation.use_all_gps);
      AppendFingerprintField(oss, "include_galileo", cfg.constellation.include_galileo);
      AppendFingerprintVector(oss, "gps_prns", cfg.constellation.gps_prns);
      AppendFingerprintVector(oss, "galileo_prns", cfg.constellation.galileo_prns);
      AppendFingerprintField(oss, "design_name", cfg.design.name);
      AppendFingerprintField(oss, "setup_transmitters", cfg.design.setup_transmitters);
      AppendFingerprintField(oss, "receiver_antenna_name", cfg.design.receiver_antenna_name);
      AppendFingerprintField(oss, "tx_yaw_dedicated", cfg.design.tx_yaw_dedicated);
      AppendFingerprintField(oss, "apply_cn0_threshold", cfg.design.apply_cn0_threshold);
      AppendFingerprintField(oss, "cn0_threshold_dbhz", cfg.design.cn0_threshold_dbhz);
      AppendFingerprintField(oss, "cn0_acquisition_threshold_dbhz",
                             cfg.design.cn0_acquisition_threshold_dbhz);
      AppendFingerprintField(oss, "cn0_tracking_threshold_dbhz",
                             cfg.design.cn0_tracking_threshold_dbhz);
      // Measurement-selection and filter-only settings are intentionally absent. Precompute
      // writes a superset link cache (primary + secondary frequencies, pseudorange + Doppler
      // metadata); the EKF run chooses pseudorange/Doppler/TDCP and IF/single-frequency usage
      // from the current config without changing the generated measurement geometry.
      AppendFingerprintField(oss, "moon_gravity_degree_truth", cfg.moon_gravity_degree_truth);
      AppendFingerprintField(oss, "moon_gravity_order_truth", cfg.moon_gravity_order_truth);
      AppendFingerprintField(oss, "include_earth", cfg.include_earth);
      AppendFingerprintField(oss, "include_sun", cfg.include_sun);
      AppendFingerprintField(oss, "use_relativity", cfg.use_relativity);
      AppendFingerprintField(oss, "use_srp_truth", cfg.use_srp_truth);
      AppendFingerprintField(oss, "srp_coeff_truth_m2_kg", cfg.srp_coeff_truth_m2_kg);
      return std::to_string(std::hash<std::string>{}(oss.str()));
    }

    bool LinkCacheMatchesConfig(const LunarGnssODTSConfig& cfg, const std::string& fingerprint) {
      if (!FileHasContent(cfg.links_file)) return false;
      std::ifstream meta(LinkMetadataPath(cfg.links_file));
      std::string cached;
      return static_cast<bool>(std::getline(meta, cached)) && cached == fingerprint;
    }

    void WriteLinkCacheMetadata(const LunarGnssODTSConfig& cfg, const std::string& fingerprint) {
      std::ofstream meta(LinkMetadataPath(cfg.links_file));
      meta << fingerprint << "\n";
    }

    void RemoveStaleDelayCache(const LunarGnssODTSConfig& cfg) {
      std::error_code ec;
      std::filesystem::remove(cfg.delays_file, ec);
      std::filesystem::remove(std::filesystem::path(cfg.delays_file.string() + ".tmp"), ec);
    }

    void RequireDelayTable(const LunarGnssODTSConfig& cfg) {
      if (!cfg.plasma.simulate_truth) return;
      LUPNT_CHECK(FileHasContent(cfg.delays_file),
                  "Plasma/ionosphere truth is enabled, but the delay table was not generated: "
                      + cfg.delays_file.string()
                      + ". Run `pixi run precompute-gnss-delays` before the C++ Monte Carlo stage.",
                  "LunarGnssODTS");
    }

    std::vector<std::filesystem::path> ReadPathVector(
        const YAML::Node& node, const std::string& key,
        const std::vector<std::filesystem::path>& fallback, const std::filesystem::path& base) {
      if (!node || !node[key]) return fallback;
      std::vector<std::filesystem::path> paths;
      for (const auto& item : node[key]) paths.push_back(ResolvePath(base, item.as<std::string>()));
      return paths;
    }

    enum class AppEventType { RECEIVER_GNSS_MEASUREMENT };

    struct ScheduledAppCall {
      AppEventType type = AppEventType::RECEIVER_GNSS_MEASUREMENT;
      double receiver_clock_s = 0.0;
      double elapsed_coordinate_s = 0.0;
      int sequence = 0;
    };

    struct ScheduledAppCallLater {
      bool operator()(const ScheduledAppCall& a, const ScheduledAppCall& b) const {
        if (a.elapsed_coordinate_s == b.elapsed_coordinate_s) return a.sequence > b.sequence;
        return a.elapsed_coordinate_s > b.elapsed_coordinate_s;
      }
    };

    class ReceiverApp {
    public:
      explicit ReceiverApp(const ReceiverAppConfig& config) : config_(config) {}

      std::vector<ScheduledAppCall> BuildCalls(const LunarGnssODTSConfig& cfg) const {
        const double dt_local = 1.0 / config_.rate_hz;
        const double local_start = cfg.clock_bias_s;
        const double local_end = cfg.clock_bias_s + (1.0 + cfg.clock_drift_sps) * cfg.duration_s;
        const int n = static_cast<int>(std::floor((local_end - local_start) / dt_local)) + 1;
        std::vector<ScheduledAppCall> calls;
        calls.reserve(n);
        for (int i = 0; i < n; ++i) {
          const double clock_reading_s = local_start + i * dt_local;
          calls.push_back({AppEventType::RECEIVER_GNSS_MEASUREMENT, clock_reading_s,
                           (clock_reading_s - cfg.clock_bias_s) / (1.0 + cfg.clock_drift_sps), i});
        }
        return calls;
      }

    private:
      ReceiverAppConfig config_;
    };

    class AppScheduler {
    public:
      void AddCalls(const std::vector<ScheduledAppCall>& calls) {
        for (const auto& call : calls) queue_.push(call);
      }

      std::vector<ScheduledAppCall> Drain() {
        std::vector<ScheduledAppCall> calls;
        while (!queue_.empty()) {
          calls.push_back(queue_.top());
          queue_.pop();
        }
        return calls;
      }

    private:
      std::priority_queue<ScheduledAppCall, std::vector<ScheduledAppCall>, ScheduledAppCallLater>
          queue_;
    };

    std::vector<ScheduledAppCall> BuildReceiverAppSchedule(const LunarGnssODTSConfig& cfg) {
      AppScheduler scheduler;
      scheduler.AddCalls(ReceiverApp(cfg.receiver_app).BuildCalls(cfg));
      return scheduler.Drain();
    }

    VecXd ElapsedTimesFromSchedule(const std::vector<ScheduledAppCall>& calls) {
      VecXd t(calls.size());
      for (int i = 0; i < static_cast<int>(calls.size()); ++i) {
        t(i) = calls[i].elapsed_coordinate_s;
      }
      return t;
    }

    VecXd ReceiverClockTimesFromSchedule(const std::vector<ScheduledAppCall>& calls) {
      VecXd t(calls.size());
      for (int i = 0; i < static_cast<int>(calls.size()); ++i) t(i) = calls[i].receiver_clock_s;
      return t;
    }

    std::string GnssConstName(GnssConst gnss_const) {
      switch (gnss_const) {
        case GnssConst::GPS: return "GPS";
        case GnssConst::GLONASS: return "GLONASS";
        case GnssConst::GALILEO: return "GALILEO";
        case GnssConst::BEIDOU: return "BEIDOU";
        case GnssConst::QZSS: return "QZSS";
      }
      return "UNKNOWN";
    }

    std::string GnssFreqName(GnssFreq freq) {
      switch (freq) {
        case GnssFreq::L1: return "L1";
        case GnssFreq::L2: return "L2";
        case GnssFreq::L5: return "L5";
        case GnssFreq::E1: return "E1";
        case GnssFreq::E6: return "E6";
        case GnssFreq::E5: return "E5";
        case GnssFreq::E5a: return "E5a";
        case GnssFreq::E5b: return "E5b";
      }
      return "UNKNOWN";
    }

    std::vector<GnssObservable> CurrentEpochObservables(const LunarGnssODTSConfig& cfg) {
      std::vector<GnssObservable> observables;
      if (cfg.use_pseudorange) observables.push_back(GnssObservable::PSEUDORANGE);
      if (cfg.use_doppler) observables.push_back(GnssObservable::DOPPLER);
      LUPNT_CHECK(!observables.empty() || cfg.use_tdcp,
                  "At least one GNSS measurement type must be enabled", "LunarGnssODTS");
      return observables;
    }

    std::vector<GnssObservable> PrecomputeObservables() {
      return {GnssObservable::PSEUDORANGE, GnssObservable::DOPPLER};
    }

    std::string ObservableName(GnssObservable observable) {
      switch (observable) {
        case GnssObservable::PSEUDORANGE: return "pseudorange";
        case GnssObservable::DOPPLER: return "doppler";
        case GnssObservable::CARRIER_PHASE: return "carrier_phase";
        default: return "unknown";
      }
    }

    std::string MatrixString(const MatXd& matrix) {
      std::ostringstream oss;
      Eigen::IOFormat fmt(6, Eigen::DontAlignCols, ", ", ";\n", "", "", "[", "]");
      oss << std::scientific << matrix.format(fmt);
      return oss.str();
    }

    struct LinkKey {
      int epoch_index = 0;
      int gnss_const = 0;
      int prn = 0;
      int frequency = 0;

      bool operator==(const LinkKey& other) const {
        return epoch_index == other.epoch_index && gnss_const == other.gnss_const
               && prn == other.prn && frequency == other.frequency;
      }
    };

    struct LinkKeyHash {
      std::size_t operator()(const LinkKey& key) const {
        std::size_t h = std::hash<int>{}(key.epoch_index);
        h ^= std::hash<int>{}(key.gnss_const + 0x9e3779b9 + (h << 6) + (h >> 2));
        h ^= std::hash<int>{}(key.prn + 0x9e3779b9 + (h << 6) + (h >> 2));
        h ^= std::hash<int>{}(key.frequency + 0x9e3779b9 + (h << 6) + (h >> 2));
        return h;
      }
    };

    struct DelayRecord {
      double ionosphere_plasma_delay_m = 0.0;
    };

    std::vector<std::string> SplitCsvLine(const std::string& line) {
      std::vector<std::string> cols;
      std::stringstream ss(line);
      std::string item;
      while (std::getline(ss, item, ',')) cols.push_back(item);
      return cols;
    }

    std::unordered_map<std::string, int> CsvColumnIndex(const std::string& header) {
      std::vector<std::string> columns = SplitCsvLine(header);
      std::unordered_map<std::string, int> idx;
      for (int i = 0; i < static_cast<int>(columns.size()); ++i) idx[columns[i]] = i;
      return idx;
    }

    std::string CsvColumn(const std::vector<std::string>& row,
                          const std::unordered_map<std::string, int>& idx,
                          const std::string& name) {
      auto it = idx.find(name);
      if (it == idx.end() || it->second >= static_cast<int>(row.size())) return "";
      return row[it->second];
    }

    void RequireCsvColumns(const std::unordered_map<std::string, int>& idx,
                           const std::vector<std::string>& names,
                           const std::filesystem::path& path) {
      for (const auto& name : names) {
        LUPNT_CHECK(
            idx.count(name) > 0,
            "GNSS link cache is missing required runtime column `" + name + "`: " + path.string()
                + ". Regenerate it with `pixi run python python/examples/ex6_precompute.py`.",
            "LunarGnssODTS");
      }
    }

    VecXd MakeEphemerisTimeVector(const LunarGnssODTSConfig& cfg) {
      const double margin_s = std::max(600.0, 2.0 * cfg.ephemeris_dt_s);
      const int n
          = static_cast<int>(std::floor((cfg.duration_s + 2.0 * margin_s) / cfg.ephemeris_dt_s))
            + 1;
      VecXd t(n);
      for (int i = 0; i < n; ++i) t(i) = -margin_s + i * cfg.ephemeris_dt_s;
      return t;
    }

    VecXd ShiftTimes(Real t0, const VecXd& elapsed_s) {
      VecXd times(elapsed_s.size());
      for (int i = 0; i < elapsed_s.size(); ++i) times(i) = (t0 + elapsed_s(i)).val();
      return times;
    }

    std::vector<Real> ToRealVector(const VecXd& values) {
      std::vector<Real> out;
      out.reserve(values.size());
      for (int i = 0; i < values.size(); ++i) out.push_back(values(i));
      return out;
    }

    VecXd ConvertTimeVector(const VecXd& times, Time from, Time to) {
      VecXd converted(times.size());
      for (int i = 0; i < times.size(); ++i)
        converted(i) = ConvertTime(Real(times(i)), from, to).val();
      return converted;
    }

    std::string Lowercase(std::string s) {
      std::transform(s.begin(), s.end(), s.begin(),
                     [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
      return s;
    }

    std::pair<double, double> EphemerisWindowTai(const LunarGnssODTSConfig& cfg) {
      Real t0_tdb = ConvertTime(GregorianToTime(cfg.start_epoch_utc), Time::UTC, Time::TDB);
      VecXd ephem_elapsed_s = MakeEphemerisTimeVector(cfg);
      const double t_start_tdb = (t0_tdb + ephem_elapsed_s.minCoeff()).val();
      const double t_end_tdb = (t0_tdb + ephem_elapsed_s.maxCoeff()).val();
      return {ConvertTime(t_start_tdb, Time::TDB, Time::TAI).val(),
              ConvertTime(t_end_tdb, Time::TDB, Time::TAI).val()};
    }

    std::vector<std::filesystem::path> SelectSp3FilesForEpochWindow(
        const std::filesystem::path& directory, double t_start_tai, double t_end_tai) {
      constexpr double kSp3EdgeToleranceS = 2.0 * 3600.0;
      LUPNT_CHECK(std::filesystem::exists(directory),
                  "GNSS SP3 directory not found: " + directory.string(), "LunarGnssODTS");

      std::vector<std::filesystem::path> candidates;
      for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (!entry.is_regular_file()) continue;
        if (Lowercase(entry.path().extension().string()) == ".sp3") {
          candidates.push_back(entry.path());
        }
      }
      std::sort(candidates.begin(), candidates.end());

      std::vector<std::filesystem::path> selected;
      for (const auto& path : candidates) {
        Sp3Loader loader(path);
        double file_start = std::numeric_limits<double>::infinity();
        double file_end = -std::numeric_limits<double>::infinity();
        for (const std::string& sat : loader.GetSatellites()) {
          const auto span = loader.GetTimeSpan(sat);
          file_start = std::min(file_start, span.first);
          file_end = std::max(file_end, span.second);
        }
        if (file_start <= t_end_tai + kSp3EdgeToleranceS
            && file_end >= t_start_tai - kSp3EdgeToleranceS) {
          selected.push_back(path);
        }
      }

      LUPNT_CHECK(!selected.empty(),
                  "No SP3 files overlap the requested TAI ephemeris window ["
                      + std::to_string(t_start_tai) + ", " + std::to_string(t_end_tai) + "]",
                  "LunarGnssODTS");

      Sp3Loader combined(selected);
      bool covers_window = false;
      for (const std::string& sat : combined.GetSatellites()) {
        const auto span = combined.GetTimeSpan(sat);
        if (span.first <= t_start_tai && span.second >= t_end_tai) {
          covers_window = true;
          break;
        }
      }
      LUPNT_CHECK(covers_window,
                  "Selected SP3 files do not cover the full requested TAI ephemeris window ["
                      + std::to_string(t_start_tai) + ", " + std::to_string(t_end_tai) + "]",
                  "LunarGnssODTS");
      return selected;
    }

    // When auto_select_sp3 is set and no explicit SP3 files were provided, pick the SP3
    // products in sp3_directory that cover the run's ephemeris window. Runs for both the YAML
    // and the struct-based (Python) config paths so the notebook can rely on it too.
    void ResolveAutoSelectSp3Files(LunarGnssODTSConfig& cfg) {
      if (!cfg.constellation.auto_select_sp3 || !cfg.constellation.sp3_files.empty()) return;
      const auto [t_start_tai, t_end_tai] = EphemerisWindowTai(cfg);
      cfg.constellation.sp3_files
          = SelectSp3FilesForEpochWindow(cfg.constellation.sp3_directory, t_start_tai, t_end_tai);
    }

    LunarGnssODTSConfig ResolveRuntimeConfig(LunarGnssODTSConfig cfg) {
      if (cfg.constellation.use_all_gps) cfg.constellation.gps_prns.clear();
      ResolveAutoSelectSp3Files(cfg);
      return cfg;
    }

    // Resolve the RINEX-nav (BRDC) files: an explicit list, else every `*.rnx`/`*.nav` file in
    // brdc_directory. The RinexNavLoader picks, per epoch, the navigation message with the
    // closest time-of-ephemeris, so passing the whole directory is safe.
    std::vector<std::filesystem::path> ResolveBrdcFiles(const ConstellationSourceConfig& c) {
      if (!c.brdc_files.empty()) return c.brdc_files;
      std::vector<std::filesystem::path> files;
      if (c.brdc_directory.empty() || !std::filesystem::is_directory(c.brdc_directory))
        return files;
      for (const auto& entry : std::filesystem::directory_iterator(c.brdc_directory)) {
        if (!entry.is_regular_file()) continue;
        const std::string ext = entry.path().extension().string();
        const std::string name = entry.path().filename().string();
        if (ext == ".rnx" || ext == ".nav" || name.find("_MN.") != std::string::npos)
          files.push_back(entry.path());
      }
      std::sort(files.begin(), files.end());
      return files;
    }

    // Broadcast (RINEX-nav) transmitter-ephemeris error model for the ODTS filter. The truth
    // measurements keep the precise SP3 transmitter states, while the filter (receiver) model
    // is fed the broadcast position and clock -- evaluated live from the loaded navigation
    // message parameters at each signal transmit epoch, exactly as a real receiver decodes and
    // propagates the broadcast message. The injected broadcast-minus-precise error is debiased:
    // a per-constellation systematic clock offset (median over all satellites and epochs of
    // broadcast-minus-precise clock) and, for QZSS, a per-satellite median radial orbit offset
    // are removed (see Montenbruck & Steigenberger, J. Navigation, 2018).
    class BroadcastEphemerisError {
    public:
      BroadcastEphemerisError(const std::vector<std::filesystem::path>& sp3_files,
                              const std::vector<std::filesystem::path>& brdc_files,
                              const std::filesystem::path& antex_file,
                              const std::vector<std::pair<GnssConst, int>>& sats,
                              double t_start_tai, double t_end_tai, double sample_dt_s,
                              bool debias_clock, bool debias_qzss_radial)
          : sp3_(sp3_files),
            brdc_(brdc_files),
            antex_(antex_file),
            debias_clock_(debias_clock),
            debias_qzss_radial_(debias_qzss_radial) {
        // Debiasing pass: sample precise and broadcast over the ephemeris window and reduce to
        // a per-constellation clock median and a per-QZSS-satellite radial-orbit median. A
        // coarse sampling is sufficient: these systematic offsets are near-constant in time.
        std::map<GnssConst, std::vector<double>> clock_diffs;                   // brdc - sp3 [s]
        std::map<std::pair<GnssConst, int>, std::vector<double>> radial_diffs;  // (brdc-sp3).r_hat
        const double span = std::max(t_end_tai - t_start_tai, 0.0);
        const int n = std::max(2, static_cast<int>(span / std::max(sample_dt_s, 1.0)) + 1);
        for (const auto& [gc, prn] : sats) {
          for (int k = 0; k < n; k++) {
            const double t = t_start_tai + span * k / (n - 1);
            Vec3 dr, r_brdc;
            Real dc;
            if (!RawDelta(gc, prn, GnssFreq::L1, Real(t), dr, dc, r_brdc)) continue;
            clock_diffs[gc].push_back(dc.val());
            if (gc == GnssConst::QZSS)
              radial_diffs[{gc, prn}].push_back(dr.dot(r_brdc.normalized()).val());
          }
        }
        for (auto& [gc, v] : clock_diffs) clock_median_[gc] = Median(v);
        for (auto& [key, v] : radial_diffs) radial_median_[key] = Median(v);
      }

      // Debiased broadcast-minus-precise transmitter error at transmit epoch `t_tai`. `dr_eci`
      // is the position delta [m] to add to the SP3 tx_state -- identical in any celestial-
      // inertial frame that shares J2000 axes (e.g. MOON_CI), since it is a difference of two
      // positions. `dc_s` is the clock delta [s]. Returns false if the satellite has no
      // broadcast message (then the filter keeps the precise SP3 state, i.e. no injected error).
      bool GetDebiasedDelta(GnssConst gc, int prn, GnssFreq freq, Real t_tai, Vec3& dr_eci,
                            Real& dc_s) const {
        Vec3 r_brdc;
        if (!RawDelta(gc, prn, freq, t_tai, dr_eci, dc_s, r_brdc)) return false;
        if (debias_clock_) {
          auto it = clock_median_.find(gc);
          if (it != clock_median_.end()) dc_s -= it->second;
        }
        if (debias_qzss_radial_ && gc == GnssConst::QZSS) {
          auto it = radial_median_.find({gc, prn});
          if (it != radial_median_.end()) dr_eci -= Real(it->second) * r_brdc.normalized();
        }
        return true;
      }

    private:
      // Raw broadcast-minus-precise delta (no debiasing), both antenna-phase-center, in ECI.
      // Also returns the broadcast position `r_brdc_eci` (for the QZSS radial direction).
      bool RawDelta(GnssConst gc, int prn, GnssFreq freq, Real t_tai, Vec3& dr_eci, Real& dc_s,
                    Vec3& r_brdc_eci) const {
        const std::string sat_id = AntexLoader::SatId(gc, prn);
        if (!brdc_.HasSatellite(sat_id)) return false;
        Vec6 rv_sp3_ecef, rv_brdc_ecef, tmp;
        Real clk_sp3, clk_brdc;
        try {
          rv_sp3_ecef = sp3_.GetPosVel(sat_id, t_tai);
          sp3_.GetPosVelClock(sat_id, t_tai, tmp, clk_sp3);
          brdc_.GetPosVelClock(sat_id, t_tai, rv_brdc_ecef, clk_brdc);
        } catch (const std::exception&) {
          return false;  // outside the ephemeris span, or no covering nav message
        }
        // Precise SP3 center-of-mass -> antenna phase center, matching the constellation build.
        const Vec3d pos_sp3_ecef(rv_sp3_ecef(0).val(), rv_sp3_ecef(1).val(), rv_sp3_ecef(2).val());
        const Vec3d pco = antex_.HasPco(gc, prn, freq, t_tai) ? antex_.GetPco(gc, prn, freq, t_tai)
                                                              : Vec3d::Zero();
        const Vec3d pos_sp3_apc = AntexLoader::ApplyPcoCorrectionEcef(t_tai, pos_sp3_ecef, pco);
        Vec6 sp3_apc_ecef;
        sp3_apc_ecef << Real(pos_sp3_apc(0)), Real(pos_sp3_apc(1)), Real(pos_sp3_apc(2)),
            rv_sp3_ecef(3), rv_sp3_ecef(4), rv_sp3_ecef(5);

        const Real t_tdb = ConvertTime(t_tai, Time::TAI, Time::TDB);
        const Vec6 rv_sp3_eci = ConvertFrame(t_tdb, sp3_apc_ecef, Frame::ECEF, Frame::ECI, false);
        const Vec6 rv_brdc_eci = ConvertFrame(t_tdb, rv_brdc_ecef, Frame::ECEF, Frame::ECI, false);
        r_brdc_eci = rv_brdc_eci.head(3);
        dr_eci = rv_brdc_eci.head(3) - rv_sp3_eci.head(3);
        dc_s = clk_brdc - clk_sp3;
        return true;
      }

      static double Median(std::vector<double> v) {
        if (v.empty()) return 0.0;
        std::sort(v.begin(), v.end());
        const size_t m = v.size() / 2;
        return v.size() % 2 == 1 ? v[m] : 0.5 * (v[m - 1] + v[m]);
      }

      Sp3Loader sp3_;
      RinexNavLoader brdc_;
      AntexLoader antex_;
      bool debias_clock_;
      bool debias_qzss_radial_;
      std::map<GnssConst, double> clock_median_;
      std::map<std::pair<GnssConst, int>, double> radial_median_;
    };

    bool EstimateSrp(const LunarGnssODTSConfig& cfg) { return cfg.estimate_srp_coefficient; }

    bool UseFilterSrp(const LunarGnssODTSConfig& cfg) {
      return cfg.use_srp_filter || cfg.estimate_srp_coefficient;
    }

    double TruthSrpCoeff(const LunarGnssODTSConfig& cfg) {
      return cfg.use_srp_truth ? cfg.srp_coeff_truth_m2_kg : 0.0;
    }

    Ptr<NBodyDynamics> CreateLunarDynamics(int moon_degree, int moon_order, bool use_autodiff,
                                           const LunarGnssODTSConfig& cfg, bool use_srp = false,
                                           double srp_coeff_m2_kg = 0.0) {
      auto dynamics = MakePtr<NBodyDynamics>();
      dynamics->SetIntegrator(IntegratorType::RKF45);
      dynamics->SetIntegratorParams(IntegratorParams(20, 1.0e-12, 1.0e-12));
      dynamics->AddBody(Body::Moon(moon_degree, moon_order));
      if (cfg.include_earth) dynamics->AddBody(Body::Earth());
      if (cfg.include_sun) dynamics->AddBody(Body::Sun());
      dynamics->SetFrame(Frame::MOON_CI);
      dynamics->SetTimeStep(cfg.integration_step_s);
      dynamics->SetAutodiff(use_autodiff);
      dynamics->SetUseRelativity(cfg.use_relativity);
      if (use_srp) dynamics->SetSrpCoeff(srp_coeff_m2_kg);
      return dynamics;
    }

    Vec6 ReceiverElfoInitialState(Real t0_tdb, const LunarGnssODTSConfig& cfg) {
      Vec6 coe_moon_op(cfg.receiver_a_m, cfg.receiver_ecc, cfg.receiver_inc_rad,
                       cfg.receiver_raan_rad, cfg.receiver_argp_rad, cfg.receiver_mean_anomaly_rad);
      Vec6 rv0_op = ClassicalToCart(coe_moon_op, GM_MOON);
      return ConvertFrame(t0_tdb, rv0_op, Frame::MOON_OP, Frame::MOON_CI);
    }

    // Receiver clock model from the config string (e.g. "OCXO", "RAFS", "CSAC", ...). The
    // noise coefficients per model live in ClockDynamics::GetClockValues.
    ClockModel ClockModelFromConfig(const LunarGnssODTSConfig& cfg) {
      auto model = enum_cast<ClockModel>(cfg.clock_model);
      LUPNT_CHECK(model.has_value() && model.value() != ClockModel::UNDEFINED,
                  "Unknown clock_model '" + cfg.clock_model
                      + "' (expected one of OCXO, USO, CSAC, MINI_RAFS, RAFS, DSAC)",
                  "LunarGnssODTS");
      return model.value();
    }

    Ptr<ClockDynamics> CreateClockDynamics(ClockModel model, bool add_noise, int seed = 0) {
      auto clock = MakePtr<ClockDynamics>();
      clock->SetModel(model);
      clock->SetClockBiasUnit(ClockBiasUnit::METERS);
      clock->SetAddNoise(add_noise);
      clock->SetSeed(seed);
      return clock;
    }

    Ptr<JointOrbitClockDynamics> CreateJointDynamics(Ptr<NumericalDynamics> orbit_dynamics,
                                                     ClockModel clock_model, bool add_clock_noise,
                                                     int seed = 0) {
      auto joint = MakePtr<JointOrbitClockDynamics>();
      joint->SetOrbitDynamics(orbit_dynamics);
      joint->SetClockDynamics(CreateClockDynamics(clock_model, add_clock_noise, seed));
      joint->SetUseClockRelativity(true);
      joint->SetRelativityCenterBody(BodyId::MOON);
      joint->SetAddClockNoise(add_clock_noise);
      joint->SetFrame(Frame::MOON_CI);
      joint->SetTimeStep(orbit_dynamics->GetTimeStep());
      joint->SetIntegrator(IntegratorType::RKF45);
      joint->SetIntegratorParams(IntegratorParams(20, 1.0e-12, 1.0e-12));
      return joint;
    }

    JointOrbitClockState MakeJointState(const Vec6& rv, double clock_bias_s, double clock_drift_sps,
                                        bool three_state = false,
                                        double clock_drift_rate_sps2 = 0.0) {
      Cart6 orbit(rv, Frame::MOON_CI);
      if (three_state) {
        ClockState3 clock;  // [bias, drift, drift-rate] in range-equivalent units
        clock.b() = ClockDynamics::SecondsToBiasUnits(clock_bias_s, ClockBiasUnit::METERS);
        clock.d() = ClockDynamics::SecondsToBiasUnits(clock_drift_sps, ClockBiasUnit::METERS);
        clock.dr()
            = ClockDynamics::SecondsToBiasUnits(clock_drift_rate_sps2, ClockBiasUnit::METERS);
        clock.SetUnits(ClockDynamics::GetStateUnits(3, ClockBiasUnit::METERS));
        return JointOrbitClockState(orbit, clock);
      }
      ClockState2 clock;
      clock.b() = ClockDynamics::SecondsToBiasUnits(clock_bias_s, ClockBiasUnit::METERS);
      clock.d() = ClockDynamics::SecondsToBiasUnits(clock_drift_sps, ClockBiasUnit::METERS);
      clock.SetUnits(ClockDynamics::GetStateUnits(2, ClockBiasUnit::METERS));
      return JointOrbitClockState(orbit, clock);
    }

    VecXd SampleInitialError(const LunarGnssODTSConfig& cfg);

    JointOrbitClockState ExtractJointOrbitClockState(const State& x) {
      if (x.size() == 8) return JointOrbitClockState(x);
      LUPNT_CHECK(x.size() == 9, "GNSS filter state must be orbit/clock or orbit/clock/SRP",
                  "LunarGnssODTS");
      State x_joint(8);
      x_joint = x.head(8);
      x_joint.SetFrame(x.GetFrame());
      x_joint.SetName("JointOrbitClock");
      x_joint.SetNames({"r_x", "r_y", "r_z", "v_x", "v_y", "v_z", "b", "d"});
      x_joint.SetUnits({"m", "m", "m", "m/s", "m/s", "m/s", "m", "m/s"});
      return JointOrbitClockState(x_joint);
    }

    State MakeAugmentedSrpState(const State& joint_state, double srp_coeff_m2_kg) {
      State x(9);
      x.head(8) = joint_state.head(8);
      x(8) = srp_coeff_m2_kg;
      x.SetFrame(joint_state.GetFrame());
      x.SetName("JointOrbitClockSrp");
      x.SetNames({"r_x", "r_y", "r_z", "v_x", "v_y", "v_z", "b", "d", "bcoeff_srp"});
      x.SetUnits({"m", "m", "m", "m/s", "m/s", "m/s", "m", "m/s", "m^2/kg"});
      return x;
    }

    State MakeInitialEstimateState(const State& truth0, const LunarGnssODTSConfig& cfg) {
      State x0 = EstimateSrp(cfg) ? MakeAugmentedSrpState(truth0, cfg.initial_srp_coeff_m2_kg)
                                  : State(truth0);
      x0 += SampleInitialError(cfg).cast<Real>();
      return x0;
    }

    ParamState SrpParamState(double srp_coeff_m2_kg) {
      Vec2 params;
      params << srp_coeff_m2_kg, 0.0;
      return ParamState(params, {"bcoeff_srp", "bcoeff_drag"});
    }

    State PropagateOrbitClock(const State& x0, Real t0, Real tf, JointOrbitClockDynamics& dynamics,
                              const LunarGnssODTSConfig& cfg, MatXd* F) {
      const bool estimate_srp = EstimateSrp(cfg);
      const bool use_srp = estimate_srp || cfg.use_srp_filter;
      const double srp_coeff = estimate_srp ? x0(8).val() : cfg.srp_coeff_filter_m2_kg;
      JointOrbitClockState x_joint = ExtractJointOrbitClockState(x0);
      State xf_state;

      if (F != nullptr) {
        F->setZero(x0.size(), x0.size());
        MatXd stm_state;
        if (use_srp) {
          MatXd stm_param;
          xf_state = dynamics.PropagateWithParams(x_joint, t0, tf, SrpParamState(srp_coeff),
                                                  nullptr, &stm_state, &stm_param);
          F->block(0, 0, 8, 8) = stm_state;
          if (estimate_srp) F->block(0, 8, 8, 1) = stm_param.col(0);
        } else {
          xf_state = dynamics.Propagate(x_joint, t0, tf, nullptr, &stm_state);
          F->block(0, 0, 8, 8) = stm_state;
        }
        if (estimate_srp) {
          (*F)(8, 8) = 1.0;
        }
      } else if (use_srp) {
        xf_state = dynamics.PropagateWithParams(x_joint, t0, tf, SrpParamState(srp_coeff), nullptr);
      } else {
        xf_state = dynamics.Propagate(x_joint, t0, tf, nullptr);
      }

      State xf = estimate_srp ? MakeAugmentedSrpState(xf_state, srp_coeff) : State(xf_state);
      return State(x0, xf);
    }

    FilterDynamicsFunction CreateDynamicsFunction(Ptr<JointOrbitClockDynamics> dynamics,
                                                  LunarGnssODTSConfig cfg) {
      return [dynamics, cfg](const State& x, Real t0, Real tf, const State*, MatXd* F) {
        LUPNT_CHECK(dynamics != nullptr, "Filter dynamics not configured", "LunarGnssODTS");
        State xf = PropagateOrbitClock(x, t0, tf, *dynamics, cfg, F);
        if (F != nullptr) {
          LUPNT_CHECK(F->rows() == x.size() && F->cols() == x.size(),
                      "Filter STM has incorrect size", "LunarGnssODTS");
        }
        return xf;
      };
    }

    // Process noise for the UDU filter. The position/velocity state-noise-compensation block
    // and the 2-state clock block are both *correlated* (dense), but the UDU time update needs
    // a diagonal Q with correlated noise routed through a mapping matrix G (Q_block = U D U^T,
    // so G <- U and diagonal Q <- D). We UDU-decompose each correlated block and install G on
    // the filter via SetProcessNoiseMappingMatrix each step (G depends on dt). This mirrors the
    // reference implementation in projects/Plasmasphere_Delay_Filtering. `filter` is a raw
    // pointer to avoid a shared_ptr cycle (the filter owns this closure).
    ProcessNoiseFunction CreateProcessNoiseFunction(const LunarGnssODTSConfig& cfg,
                                                    UDUEKF* filter) {
      const ClockModel clock_model = ClockModelFromConfig(cfg);
      return [cfg, clock_model, filter](const State& x, Real t0, Real tf) {
        const double dt = std::abs((tf - t0).val());
        const int n = static_cast<int>(x.size());

        // Correlated blocks: position/velocity (SNC) and 2-state clock.
        Mat3d Q_acc = std::pow(cfg.process_accel_sigma_mps2, 2) * Mat3d::Identity();
        MatXd Q_rv = ProcessNoisePosVel(Q_acc, dt);  // 6x6, dense
        MatXd Q_clk
            = ClockDynamics::TwoStateNoise(clock_model, dt, ClockBiasUnit::METERS).cast<double>();

        // Factor each correlated block as U diag(D) U^T: D -> diagonal Q, U -> mapping G.
        MatXd G = MatXd::Identity(n, n);
        VecXd Q_diag = VecXd::Zero(n);

        VecMatPair du_rv = UDUDecomposition(Q_rv);
        Q_diag.head(6) = du_rv.first;
        G.block(0, 0, 6, 6) = du_rv.second;

        VecMatPair du_clk = UDUDecomposition(Q_clk);
        Q_diag.segment(6, 2) = du_clk.first;
        G.block(6, 6, 2, 2) = du_clk.second;

        // SRP coefficient: scalar noise, already diagonal (identity mapping).
        if (EstimateSrp(cfg)) {
          Q_diag(8) = std::pow(cfg.process_srp_coeff_sigma_m2_kg_sqrt_s, 2) * dt;
        }

        filter->SetProcessNoiseMappingMatrix(G);
        return MatXd(Q_diag.asDiagonal());
      };
    }

    MatXd InitialCovariance(const LunarGnssODTSConfig& cfg) {
      const int n = EstimateSrp(cfg) ? 9 : 8;
      MatXd P = MatXd::Zero(n, n);
      P.block(0, 0, 3, 3) = std::pow(cfg.initial_position_sigma_m, 2) * Mat3d::Identity();
      P.block(3, 3, 3, 3) = std::pow(cfg.initial_velocity_sigma_mps, 2) * Mat3d::Identity();
      P(6, 6) = std::pow(C * cfg.initial_clock_bias_sigma_s, 2);
      P(7, 7) = std::pow(C * cfg.initial_clock_drift_sigma_sps, 2);
      if (EstimateSrp(cfg)) P(8, 8) = std::pow(cfg.initial_srp_coeff_sigma_m2_kg, 2);
      return P;
    }

    VecXd SampleInitialError(const LunarGnssODTSConfig& cfg) {
      VecXd err(EstimateSrp(cfg) ? 9 : 8);
      err.setZero();
      err(0) = SampleNormal(0.0, cfg.initial_position_sigma_m).val();
      err(1) = SampleNormal(0.0, cfg.initial_position_sigma_m).val();
      err(2) = SampleNormal(0.0, cfg.initial_position_sigma_m).val();
      err(3) = SampleNormal(0.0, cfg.initial_velocity_sigma_mps).val();
      err(4) = SampleNormal(0.0, cfg.initial_velocity_sigma_mps).val();
      err(5) = SampleNormal(0.0, cfg.initial_velocity_sigma_mps).val();
      err(6) = SampleNormal(0.0, C * cfg.initial_clock_bias_sigma_s).val();
      err(7) = SampleNormal(0.0, C * cfg.initial_clock_drift_sigma_sps).val();
      if (EstimateSrp(cfg)) err(8) = SampleNormal(0.0, cfg.initial_srp_coeff_sigma_m2_kg).val();
      return err;
    }

    std::string FormatDuration(double seconds);  // defined below (progress helpers)

    struct RuntimeConstellation {
      Ptr<GnssConstellation> constellation;
      GnssFreq frequency = GnssFreq::L1;
    };

    RuntimeConstellation BuildConstellationFromFiles(GnssConst gnss_const, GnssFreq frequency,
                                                     const std::vector<int>& prns,
                                                     const LunarGnssODTSConfig& cfg,
                                                     const VecXd& ephem_times_tai) {
      auto constellation = MakePtr<GnssConstellation>(gnss_const);
      constellation->SetupSatelliteStatesFromFiles(cfg.constellation.sp3_files,
                                                   cfg.constellation.antex_file, ephem_times_tai,
                                                   frequency, prns);
      if (cfg.design.setup_transmitters) constellation->SetupTransmitters();
      return {constellation, frequency};
    }

    std::vector<RuntimeConstellation> BuildConstellations(const LunarGnssODTSConfig& cfg,
                                                          const VecXd& ephem_times_tai) {
      std::vector<RuntimeConstellation> out;
      Logger::Info(
          "Setting up GNSS constellations from SP3 (Chebyshev fit; grows with run "
          "length)...",
          "LunarGnssODTS");
      const auto constellation_setup_start = std::chrono::steady_clock::now();
      LUPNT_CHECK(!cfg.constellation.sp3_files.empty(), "GNSS filtering requires SP3 files",
                  "LunarGnssODTS");
      for (const auto& sp3 : cfg.constellation.sp3_files) {
        LUPNT_CHECK(std::filesystem::exists(sp3), "GNSS SP3 file not found: " + sp3.string(),
                    "LunarGnssODTS");
      }
      LUPNT_CHECK(std::filesystem::exists(cfg.constellation.antex_file),
                  "GNSS ANTEX file not found: " + cfg.constellation.antex_file.string(),
                  "LunarGnssODTS");

      // Collect the (constellation, frequency) sets to build, then set them up in parallel:
      // GPS L1/L5 and Galileo E1/E5a are always prepared so a link cache can be reused when the
      // EKF switches between ionosphere-free and single-frequency measurement combinations.
      // Each SP3 Chebyshev fit / transmitter setup is independent; SPICE calls inside are
      // guarded by `omp critical`. Order is preserved by writing into a pre-sized vector by
      // index (matches the serial build).
      struct ConstellationSpec {
        GnssConst gnss_const;
        GnssFreq frequency;
        std::vector<int> prns;
      };
      std::vector<ConstellationSpec> specs;
      const std::vector<int> gps_prns = cfg.constellation.gps_prns;
      specs.push_back({GnssConst::GPS, GnssFreq::L1, gps_prns});
      specs.push_back({GnssConst::GPS, GnssFreq::L5, gps_prns});
      if (cfg.constellation.include_galileo) {
        const std::vector<int> galileo_prns = cfg.constellation.galileo_prns;
        specs.push_back({GnssConst::GALILEO, GnssFreq::E1, galileo_prns});
        specs.push_back({GnssConst::GALILEO, GnssFreq::E5a, galileo_prns});
      }

      const int n_specs = static_cast<int>(specs.size());
      out.resize(n_specs);
      const int n_threads = cfg.precompute_num_threads > 0
                                ? std::min(cfg.precompute_num_threads, n_specs)
                                : n_specs;
#pragma omp parallel for schedule(dynamic) num_threads(n_threads) if (n_specs > 1)
      for (int i = 0; i < n_specs; ++i) {
        out[i] = BuildConstellationFromFiles(specs[i].gnss_const, specs[i].frequency, specs[i].prns,
                                             cfg, ephem_times_tai);
      }

      const double setup_s = std::chrono::duration<double>(std::chrono::steady_clock::now()
                                                           - constellation_setup_start)
                                 .count();
      Logger::Info(std::to_string(out.size()) + " GNSS constellation-frequency set(s) ready in "
                       + FormatDuration(setup_s),
                   "LunarGnssODTS");
      return out;
    }

    std::string FormatDuration(double seconds) {
      std::ostringstream oss;
      oss << std::fixed << std::setprecision(1);
      if (seconds >= 90.0)
        oss << (seconds / 60.0) << " min";
      else
        oss << seconds << " s";
      return oss.str();
    }

    // Time-throttled progress logger with an ETA and a wall-clock finish estimate. Prints a
    // "starting" line, then at most one line every `interval_s` seconds, then a final line.
    class ProgressReporter {
    public:
      ProgressReporter(std::string label, long total, double interval_s = 3.0)
          : label_(std::move(label)),
            total_(total),
            interval_s_(interval_s),
            start_(std::chrono::steady_clock::now()),
            last_(start_) {
        Logger::Info(label_ + ": starting (" + std::to_string(total_) + " steps)", "LunarGnssODTS");
      }

      void Update(long done) {
        const auto now = std::chrono::steady_clock::now();
        const bool finished = done >= total_;
        if (!finished && std::chrono::duration<double>(now - last_).count() < interval_s_) return;
        last_ = now;
        const double elapsed = std::chrono::duration<double>(now - start_).count();
        const double frac
            = total_ > 0 ? static_cast<double>(done) / static_cast<double>(total_) : 1.0;
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(1) << label_ << ": " << done << "/" << total_ << " ("
            << (100.0 * frac) << "%), elapsed " << FormatDuration(elapsed);
        if (finished) {
          oss << " -- done";
        } else if (frac > 1e-6) {
          const double eta = elapsed * (1.0 - frac) / frac;
          const std::time_t finish = std::time(nullptr) + static_cast<std::time_t>(eta);
          char buf[16];
          std::strftime(buf, sizeof(buf), "%H:%M:%S", std::localtime(&finish));
          oss << ", ETA " << FormatDuration(eta) << " (finish ~" << buf << ")";
        }
        Logger::Info(oss.str(), "LunarGnssODTS");
      }

    private:
      std::string label_;
      long total_;
      double interval_s_;
      std::chrono::steady_clock::time_point start_;
      std::chrono::steady_clock::time_point last_;
    };

    std::vector<State> BuildTruthStates(const LunarGnssODTSConfig& cfg, Real t0_tdb,
                                        const VecXd& elapsed_s, const VecXd& times_tdb,
                                        bool add_clock_noise, int seed,
                                        const std::string& progress_label = "") {
      (void)elapsed_s;
      Ptr<NBodyDynamics> orbit_dynamics
          = CreateLunarDynamics(cfg.moon_gravity_degree_truth, cfg.moon_gravity_order_truth, false,
                                cfg, cfg.use_srp_truth, cfg.srp_coeff_truth_m2_kg);
      Ptr<JointOrbitClockDynamics> dynamics
          = CreateJointDynamics(orbit_dynamics, ClockModelFromConfig(cfg), add_clock_noise, seed);
      Vec6 rv0 = ReceiverElfoInitialState(t0_tdb, cfg);
      State x = MakeJointState(rv0, cfg.clock_bias_s, cfg.clock_drift_sps,
                               cfg.use_three_state_clock_truth, cfg.clock_drift_rate_sps2);
      std::vector<State> states;
      states.reserve(times_tdb.size());
      std::unique_ptr<ProgressReporter> prog;
      if (!progress_label.empty())
        prog = std::make_unique<ProgressReporter>(progress_label, times_tdb.size(),
                                                  cfg.precompute_progress_interval_s);
      for (int i = 0; i < times_tdb.size(); ++i) {
        if (i == 0) {
          states.push_back(x);
        } else {
          x = dynamics->Propagate(JointOrbitClockState(x), times_tdb(i - 1), times_tdb(i), nullptr);
          states.push_back(x);
        }
        if (prog) prog->Update(i + 1);
      }
      return states;
    }

    GnssMeasurementOptions MeasurementOptions(const LunarGnssODTSConfig& cfg, bool truth_model) {
      GnssMeasurementOptions options;
      options.observables = CurrentEpochObservables(cfg);
      options.clock_bias_unit = ClockBiasUnit::METERS;
      options.receive_time_scale = Time::TDB;
      options.ephemeris_time_scale = Time::TAI;
      options.frame = Frame::MOON_CI;
      options.solve_light_time = true;
      options.recompute_transmit_time_from_ephemeris = !truth_model;
      options.apply_transmitter_relativity = false;
      options.apply_shapiro_delay = truth_model;
      options.apply_visibility = true;
      options.apply_cn0_threshold = cfg.design.apply_cn0_threshold;
      // Separate acquisition/tracking C/N0 thresholds implement the standard lock hysteresis
      // (a satellite must exceed the acquisition threshold to be acquired, then stays tracked
      // while above the lower tracking threshold). Each defaults to the deprecated single
      // `cn0_threshold_dbhz` when left < 0, preserving the old single-threshold behavior.
      const double acq = cfg.design.cn0_acquisition_threshold_dbhz >= 0.0
                             ? cfg.design.cn0_acquisition_threshold_dbhz
                             : cfg.design.cn0_threshold_dbhz;
      const double trk = cfg.design.cn0_tracking_threshold_dbhz >= 0.0
                             ? cfg.design.cn0_tracking_threshold_dbhz
                             : cfg.design.cn0_threshold_dbhz;
      options.cn0_threshold_dbhz = cfg.design.cn0_threshold_dbhz;
      options.cn0_acquisition_threshold_dbhz = acq;
      options.cn0_tracking_threshold_dbhz = trk;
      options.apply_ionosphere_plasma_delay = truth_model && cfg.plasma.simulate_truth;
      options.tx_yaw_model = cfg.design.tx_yaw_dedicated
                                 ? GnssMeasurementOptions::TxYawModel::DEDICATED
                                 : GnssMeasurementOptions::TxYawModel::NOMINAL;
      return options;
    }

    std::vector<GnssOccludingBody> LunarGnssOccludingBodies(Real t_tdb) {
      return {
          {Real(R_MOON), Vec3::Zero()},
          {Real(R_EARTH), GetBodyPos(t_tdb, BodyId::MOON, BodyId::EARTH, Frame::MOON_CI)},
      };
    }

    GnssMeasurementOptions CarrierMeasurementOptions(const LunarGnssODTSConfig& cfg,
                                                     bool truth_model) {
      GnssMeasurementOptions options = MeasurementOptions(cfg, truth_model);
      options.observables = {GnssObservable::CARRIER_PHASE};
      return options;
    }

    GnssIonospherePlasmaRayTraceOptions RayTraceOptions(const LunarGnssODTSConfig& cfg) {
      GnssIonospherePlasmaRayTraceOptions options;
      options.config.step_size = cfg.plasma.raytrace_step_size_km;
      options.config.correction = cfg.plasma.raytrace_correction;
      options.config.fine_correction = cfg.plasma.raytrace_fine_correction;
      options.config.cutoff_r = cfg.plasma.raytrace_cutoff_radius_re * pecsim::RE;
      options.config.gradn_dx = cfg.plasma.raytrace_gradient_step_km;
      options.config.integ_method = cfg.plasma.raytrace_integrator;
      options.config.correction_method = cfg.plasma.raytrace_correction_method;
      options.config.kp = cfg.plasma.raytrace_kp;
      options.config.rz12 = cfg.plasma.raytrace_rz12;
      options.config.use_fortran_gcpm = cfg.plasma.raytrace_use_fortran_gcpm;
      options.config.corr_tol = cfg.plasma.raytrace_correction_tolerance_m;
      options.config.compute_higher_order = cfg.plasma.raytrace_compute_higher_order;
      options.config.use_adaptive_step = cfg.plasma.raytrace_use_adaptive_step;
      options.config.straight_ray = cfg.plasma.raytrace_straight_ray;
      options.raytrace_frame = Frame::ECEF;
      options.raytrace_epoch_scale = Time::UTC;
      options.use_channel_frequency = true;
      options.iri_model = pecsim::IRIModel::IRI_2007;
      options.delay_mode = cfg.plasma.raytrace_compute_higher_order
                               ? GnssIonospherePlasmaRayTraceDelayMode::TEC_PLUS_HIGHER_ORDER
                               : GnssIonospherePlasmaRayTraceDelayMode::TEC_ONLY;
      return options;
    }

    std::vector<GnssChannel> MakeFilterChannels(std::vector<GnssChannel> channels,
                                                const LunarGnssODTSConfig& cfg) {
      if (!cfg.plasma.model_in_filter) {
        for (auto& channel : channels) channel.ionosphere_plasma_delay_m = 0.0;
      }
      // The Shapiro (gravitational) delay is modeled by the filter: it is deterministic and
      // essentially insensitive to the ~100 m filter position error, so the value precomputed
      // from the link geometry is kept rather than zeroed.
      return channels;
    }

    void ApplyMeasurementSigmas(std::vector<GnssChannel>& channels, const LunarGnssODTSConfig& cfg,
                                bool filter_covariance) {
      for (auto& channel : channels) {
        if (cfg.design.use_cn0_measurement_sigmas) {
          LUPNT_CHECK(std::isfinite(channel.sigma_pseudorange_m.val())
                          && std::isfinite(channel.sigma_doppler_hz.val())
                          && std::isfinite(channel.sigma_carrier_phase_cycles.val()),
                      "C/N0-derived GNSS measurement sigmas are missing. Ensure the selected "
                      "design sets link_budget.setup_transmitters=true and that transmitter "
                      "power/antenna metadata exist for every selected GNSS signal.",
                      "LunarGnssODTS");
        } else {
          if (!std::isfinite(channel.sigma_pseudorange_m.val()))
            channel.sigma_pseudorange_m = cfg.pseudorange_sigma_m;
          if (!std::isfinite(channel.sigma_doppler_hz.val()))
            channel.sigma_doppler_hz = cfg.doppler_sigma_hz;
          if (!std::isfinite(channel.sigma_carrier_phase_cycles.val()))
            channel.sigma_carrier_phase_cycles = cfg.carrier_phase_sigma_m / channel.Wavelength();
        }
        if (filter_covariance) {
          // In ionosphere-free mode the pseudorange inflation is applied to the *combined*
          // IF sigma in Step (via InflatePseudorangeSigma), not to the raw L1/L5 sigmas here
          // (otherwise it would be amplified by the IF combination). For single-frequency
          // pseudoranges the plasma-mismatch and the general filter inflation apply here.
          if (!cfg.use_ionosphere_free) {
            channel.sigma_pseudorange_m = std::hypot(
                channel.sigma_pseudorange_m.val(), cfg.plasma.filter_pseudorange_noise_inflation_m);
            channel.sigma_pseudorange_m = std::hypot(channel.sigma_pseudorange_m.val(),
                                                     cfg.filter_pseudorange_noise_inflation_m);
          }
          channel.sigma_doppler_hz = std::hypot(channel.sigma_doppler_hz.val(),
                                                cfg.plasma.filter_doppler_noise_inflation_hz);
        }
      }
    }

    // The lunar-GNSS combined + TDCP measurement model lives in
    // `LunarGnssCombinedMeasurement`; the truth-data path below reuses its static helpers.
    using TdcpPair = LunarGnssCombinedMeasurement::TdcpPair;

    // One row per filter measurement: its type (pseudorange/doppler/carrier_phase/TDCP) and the
    // satellite it comes from, in the same order they are stacked into the measurement vector.
    struct MeasRowDesc {
      std::string type;
      std::string sat;      // e.g. "GPS PRN 18 L1" or "GPS PRN 18 L1+L5" for dual-frequency
      double iono_m = 0.0;  // ionosphere/plasma delay on the signal: per-channel for the current
                            // observables (ionosphere-free -> ~0), epoch-to-epoch change for TDCP.
      // Range-direction transmitter ephemeris error: the (true SP3 - broadcast) transmitter
      // position error projected onto the receiver->transmitter line of sight [m] -- i.e. the
      // broadcast orbit error the pseudorange actually sees. NaN for TDCP rows.
      double eph_range_err_m = std::numeric_limits<double>::quiet_NaN();
      // Broadcast transmitter clock error in range units [m], true - broadcast (the other half
      // of the broadcast SISE). NaN for TDCP rows.
      double eph_clock_err_m = std::numeric_limits<double>::quiet_NaN();
    };

    // Per-constellation dual-frequency pair for the ionosphere-free combination and TDCP:
    // GPS uses (L1, L5); Galileo uses (E1, E5a). E1/E5a share the L1/L5 carrier frequencies,
    // so the IF math is identical -- only the frequency-enum selection differs.
    GnssFreq PrimaryFrequency(GnssConst gnss_const) {
      return gnss_const == GnssConst::GALILEO ? GnssFreq::E1 : GnssFreq::L1;
    }
    GnssFreq SecondaryFrequency(GnssConst gnss_const) {
      return gnss_const == GnssConst::GALILEO ? GnssFreq::E5a : GnssFreq::L5;
    }

    // Range-direction transmitter ephemeris error: the (true SP3 - broadcast) transmitter position
    // error projected onto the receiver->transmitter line of sight [m] -- the broadcast orbit error
    // the pseudorange actually sees. All positions are in MOON_CI (the measurement options frame),
    // so this is a same-frame projection. `filter_ch` is the truth channel with the broadcast delta
    // injected, so (truth_tx - filter_tx) is exactly the negated broadcast orbit error.
    double LosEphError(const GnssChannel& truth_ch, const GnssChannel& filter_ch, const Vec3d& rx) {
      Vec3d tx_true = truth_ch.tx_state.head(3).cast<double>();
      Vec3d tx_eph = filter_ch.tx_state.head(3).cast<double>();
      Vec3d los = (tx_true - rx).normalized();
      return (tx_true - tx_eph).dot(los);
    }

    // Broadcast transmitter clock error in range units [m], same (true - broadcast) sign
    // convention as LosEphError: C * (precise - broadcast tx clock bias). This is the *other*
    // half of the broadcast signal-in-space error the pseudorange sees (it enters rho as the
    // c*dt_tx term); the debiasing removes only the per-constellation median, not this residual.
    double ClockEphError(const GnssChannel& truth_ch, const GnssChannel& filter_ch) {
      return C * (truth_ch.tx_clock_bias_s - filter_ch.tx_clock_bias_s).val();
    }

    // Descriptors are built from the *truth* channels/pairs so `iono_m` reflects the actual
    // signal delay (the filter channels have it zeroed unless plasma.model_in_filter). Truth and
    // filter channel sets share the same satellites in the same order, so the row order matches
    // the stacked measurement vector. `dual_frequency` labels the current-epoch (pseudorange/
    // Doppler) rows as the ionosphere-free L1+L5 / E1+E5a combination; TDCP is always L1/E1.
    std::vector<MeasRowDesc> MeasurementRowDescriptors(
        const std::vector<GnssChannel>& channels, const std::vector<GnssChannel>& filter_channels,
        const Vec3d& rx_pos_moonci, const GnssMeasurementOptions& current_options,
        const std::vector<TdcpPair>& tdcp_pairs, bool dual_frequency) {
      auto sat_name = [](const GnssChannel& c, bool dual) {
        const std::string freq = dual ? GnssFreqName(PrimaryFrequency(c.gnss_const)) + "+"
                                            + GnssFreqName(SecondaryFrequency(c.gnss_const))
                                      : GnssFreqName(c.frequency);
        return GnssConstName(c.gnss_const) + " PRN " + std::to_string(c.prn) + " " + freq;
      };
      std::vector<MeasRowDesc> rows;
      rows.reserve(channels.size() * current_options.observables.size() + tdcp_pairs.size());
      const double kNaN = std::numeric_limits<double>::quiet_NaN();
      for (std::size_t c = 0; c < channels.size(); ++c) {
        const bool have_bc = c < filter_channels.size();
        const double eph
            = have_bc ? LosEphError(channels[c], filter_channels[c], rx_pos_moonci) : kNaN;
        const double eph_clk = have_bc ? ClockEphError(channels[c], filter_channels[c]) : kNaN;
        for (const auto observable : current_options.observables)
          rows.push_back({ObservableName(observable), sat_name(channels[c], dual_frequency),
                          channels[c].ionosphere_plasma_delay_m.val(), eph, eph_clk});
      }
      for (const auto& pair : tdcp_pairs)
        rows.push_back(
            {"TDCP", sat_name(pair.current, false),
             (pair.current.ionosphere_plasma_delay_m - pair.previous.ionosphere_plasma_delay_m)
                 .val()});
      return rows;
    }

    void PrintMatrixDiagnostics(int epoch_index, const LunarGnssODTSConfig& cfg,
                                const Ptr<UDUEKF>& filter,
                                const std::vector<GnssChannel>& truth_channels,
                                const std::vector<GnssChannel>& filter_channels,
                                const Vec3d& rx_pos_moonci,
                                const GnssMeasurementOptions& current_options,
                                const std::vector<TdcpPair>& truth_tdcp_pairs) {
      if (cfg.debug_print_matrix_epochs <= 0 || epoch_index >= cfg.debug_print_matrix_epochs)
        return;

      const MatXd F = filter->GetStateJacobian();
      const MatXd H = filter->GetMeasurementJacobian();
      const MatXd R = filter->GetMeasurementNoiseCov();
      const VecXd z_obs = filter->GetTrueMeasurement();
      const VecXd z_pred = filter->GetPredictedMeasurement();
      const VecXd dz = filter->GetMeasurementResidual();
      if (z_obs.size() == 0) {
        Logger::Info("\n[GNSS filter matrix diagnostics] epoch " + std::to_string(epoch_index)
                         + "\nNo measurements available at this epoch.",
                     "LunarGnssODTS");
        return;
      }
      const int n_rows = H.rows();
      const int n_print = std::min(n_rows, std::max(0, cfg.debug_print_matrix_max_rows));

      std::ostringstream oss;
      oss << "\n[GNSS filter matrix diagnostics] epoch " << epoch_index << "\n";
      if (F.size() > 0) {
        oss << "STM F (" << F.rows() << "x" << F.cols() << "):\n" << MatrixString(F) << "\n";
      } else {
        oss << "STM F: not available at epoch 0 before the first predict step\n";
      }
      oss << "Measurement Jacobian H: " << H.rows() << "x" << H.cols() << " (printing first "
          << n_print << " rows)\n";
      if (n_print > 0) oss << "H rows:\n" << MatrixString(H.topRows(n_print)) << "\n";

      // Per-measurement table. S = H P H^T + R is the innovation covariance (post-update S_).
      // `iono_m` is the ionosphere/plasma delay on the signal (from truth): ~0 for the
      // ionosphere-free pseudorange, and the small epoch-to-epoch change for TDCP.
      const auto rows
          = MeasurementRowDescriptors(truth_channels, filter_channels, rx_pos_moonci,
                                      current_options, truth_tdcp_pairs, cfg.use_ionosphere_free);
      const MatXd S = filter->GetInnovationCov();
      const double kNaN = std::numeric_limits<double>::quiet_NaN();
      auto cell = [](double v) {
        std::ostringstream s;
        s << std::scientific << std::setprecision(4) << v;
        return s.str();
      };
      auto diag_sqrt = [&](const MatXd& M, int i) {
        return (M.rows() > i && M.cols() > i) ? std::sqrt(std::max(0.0, M(i, i))) : kNaN;
      };
      oss << "measurements (" << n_rows << " rows):\n";
      oss << std::left << std::setw(4) << "#" << std::setw(13) << "type" << std::setw(16) << "sat"
          << std::right << std::setw(13) << "iono_m" << std::setw(14) << "eph_los_m"
          << std::setw(14) << "eph_clk_m" << std::setw(16) << "z_obs" << std::setw(16) << "z_pred"
          << std::setw(14) << "dz" << std::setw(13) << "R_sqrt" << std::setw(13) << "S_sqrt"
          << "\n";
      for (int i = 0; i < n_rows; ++i) {
        const bool have = i < static_cast<int>(rows.size());
        const std::string type = have ? rows[i].type : "?";
        const std::string sat = have ? rows[i].sat : "?";
        oss << std::left << std::setw(4) << i << std::setw(13) << type << std::setw(16) << sat
            << std::right << std::setw(13) << cell(have ? rows[i].iono_m : kNaN) << std::setw(14)
            << cell(have ? rows[i].eph_range_err_m : kNaN) << std::setw(14)
            << cell(have ? rows[i].eph_clock_err_m : kNaN) << std::setw(16)
            << cell(i < z_obs.size() ? z_obs(i) : kNaN) << std::setw(16)
            << cell(i < z_pred.size() ? z_pred(i) : kNaN) << std::setw(14)
            << cell(i < dz.size() ? dz(i) : kNaN) << std::setw(13) << cell(diag_sqrt(R, i))
            << std::setw(13) << cell(diag_sqrt(S, i)) << "\n";
      }
      Logger::Info(oss.str(), "LunarGnssODTS");
    }

    State MakeClonedState(const State& x) {
      State x_clone(2 * x.size());
      x_clone.head(x.size()) = x;
      x_clone.tail(x.size()) = x;
      x_clone.SetFrame(x.GetFrame());
      return x_clone;
    }

    MatXd MakeClonedCovariance(const MatXd& P) {
      const int n = static_cast<int>(P.rows());
      MatXd P_clone = MatXd::Zero(2 * n, 2 * n);
      P_clone.topLeftCorner(n, n) = P;
      P_clone.topRightCorner(n, n) = P;
      P_clone.bottomLeftCorner(n, n) = P;
      P_clone.bottomRightCorner(n, n) = P;
      return P_clone;
    }

    State CurrentFilterState(const State& x, const LunarGnssODTSConfig& cfg) {
      return cfg.use_tdcp ? State(x.head(x.size() / 2)) : x;
    }

    MatXd CurrentFilterCovariance(const MatXd& P, const LunarGnssODTSConfig& cfg) {
      return cfg.use_tdcp ? P.topLeftCorner(P.rows() / 2, P.cols() / 2) : P;
    }

    VecXd AddMeasurementNoise(const VecXd& y, const std::vector<GnssChannel>& channels,
                              const GnssMeasurementOptions& options) {
      VecXd noisy = y;
      const int n_obs = static_cast<int>(options.observables.size());
      if (n_obs == 0) return noisy;
      LUPNT_CHECK(y.size() == n_obs * static_cast<int>(channels.size()),
                  "GNSS measurement vector size does not match channel count", "LunarGnssODTS");
      for (int i = 0; i < y.size(); ++i) {
        const GnssChannel& channel = channels[i / n_obs];
        const GnssObservable obs = options.observables[i % n_obs];
        const double sigma
            = obs == GnssObservable::PSEUDORANGE
                  ? channel.sigma_pseudorange_m.val()
                  : (obs == GnssObservable::DOPPLER ? channel.sigma_doppler_hz.val()
                                                    : channel.sigma_carrier_phase_cycles.val());
        LUPNT_CHECK(std::isfinite(sigma) && sigma >= 0.0,
                    "GNSS measurement sigma must be finite and non-negative", "LunarGnssODTS");
        noisy(i) += SampleNormal(0.0, sigma).val();
      }
      return noisy;
    }

    VecXd AddTdcpNoise(const VecXd& y, const std::vector<TdcpPair>& pairs,
                       const LunarGnssODTSConfig& cfg) {
      VecXd noisy = y;
      MatXd R = LunarGnssCombinedMeasurement::TdcpCovariance(pairs, cfg.tdcp_sigma_m,
                                                             /*filter_tdcp_inflation_m=*/0.0);
      for (int i = 0; i < noisy.size(); ++i)
        noisy(i) += SampleNormal(0.0, std::sqrt(R(i, i))).val();
      return noisy;
    }

    Mat3d RotCartToRtn(const Vec3d& r, const Vec3d& v) {
      Vec3d r_hat = r.normalized();
      Vec3d h_hat = r_hat.cross(v.normalized()).normalized();
      Vec3d theta_hat = h_hat.cross(r_hat).normalized();
      Mat3d R;
      R.col(0) = r_hat;
      R.col(1) = theta_hat;
      R.col(2) = h_hat;
      return R;
    }

    int CountTrackedSatellites(const std::vector<GnssChannel>& channels) {
      std::set<std::pair<int, int>> sats;
      for (const auto& channel : channels) {
        sats.insert({static_cast<int>(channel.gnss_const), channel.prn});
      }
      return static_cast<int>(sats.size());
    }

    std::pair<double, double> EarthTangentRadiusAltitude(const Vec3& rx_ecef, const Vec3& tx_ecef) {
      Vec3d rx(rx_ecef(0).val(), rx_ecef(1).val(), rx_ecef(2).val());
      Vec3d tx(tx_ecef(0).val(), tx_ecef(1).val(), tx_ecef(2).val());
      Vec3d line = tx - rx;
      const double line_norm2 = line.squaredNorm();
      double u = 0.0;
      if (line_norm2 > 0.0) {
        u = std::clamp(-rx.dot(line) / line_norm2, 0.0, 1.0);
      }
      const double radius_m = (rx + u * line).norm();
      return {radius_m, radius_m - pecsim::RE * 1000.0};
    }

    // Tangent altitude [m] of the receiver->transmitter line of sight above Earth's surface.
    // The tangent *radius* is Earth-centered and therefore frame-independent, so we evaluate it
    // in ECI (cheaper than converting to ECEF). Reuses EarthTangentRadiusAltitude.
    double LinkTangentAltitudeM(Real t_tdb, const Vec3& rx_mci, const GnssChannel& channel) {
      Vec3 rx_eci = ConvertFrame(t_tdb, rx_mci, Frame::MOON_CI, Frame::ECI);
      Vec3 tx_eci = ConvertFrame(t_tdb, Vec3(channel.tx_state.head(3)), channel.frame, Frame::ECI);
      return EarthTangentRadiusAltitude(rx_eci, tx_eci).second;
    }

    // Drop channels whose LOS tangent point passes below `min_alt_m` of Earth altitude.
    std::vector<GnssChannel> FilterByTangentAltitude(const std::vector<GnssChannel>& channels,
                                                     Real t_tdb, const Vec3& rx_mci,
                                                     double min_alt_m) {
      if (min_alt_m <= 0.0) return channels;
      std::vector<GnssChannel> out;
      out.reserve(channels.size());
      for (const auto& channel : channels) {
        if (LinkTangentAltitudeM(t_tdb, rx_mci, channel) >= min_alt_m) out.push_back(channel);
      }
      return out;
    }

    // Channels on each constellation's primary frequency (L1/E1) -- used for TDCP.
    std::vector<GnssChannel> SelectPrimaryFrequency(const std::vector<GnssChannel>& channels) {
      std::vector<GnssChannel> out;
      for (const auto& channel : channels)
        if (channel.frequency == PrimaryFrequency(channel.gnss_const)) out.push_back(channel);
      return out;
    }

    // Ionosphere-free pseudorange channels: one per satellite (GPS or Galileo) that has both a
    // primary (L1/E1) and secondary (L5/E5a) channel. Geometry/ephemeris are copied from the
    // primary (identical for both frequencies). Only the primary is ray-traced; the secondary
    // delay is *derived* from the primary by first-order rescaling d_2 = d_1 * (f1/f2)^2, and
    // the IF-combined delay
    //   d_IF = f1^2/(f1^2-f2^2) * d_1 - f2^2/(f1^2-f2^2) * d_2
    // then cancels to 0 for a pure first-order model (kept explicit rather than hard-coded).
    // The measurement noise is the IF-combined C/N0 sigma
    //   sigma_IF = hypot( f1^2/(f1^2-f2^2) * sigma_1 , f2^2/(f1^2-f2^2) * sigma_2 ).
    // Satellites with only one frequency are dropped (no IF combination possible). Keys include
    // the constellation, so GPS PRN n and Galileo PRN n are treated as distinct satellites.
    std::vector<GnssChannel> MakeIonosphereFreeChannels(const std::vector<GnssChannel>& channels) {
      using SatKey = std::pair<int, int>;  // (gnss_const, prn)
      std::map<SatKey, const GnssChannel*> primary_by_sat;
      std::map<SatKey, const GnssChannel*> secondary_by_sat;
      for (const auto& channel : channels) {
        const SatKey key{static_cast<int>(channel.gnss_const), channel.prn};
        if (channel.frequency == PrimaryFrequency(channel.gnss_const))
          primary_by_sat[key] = &channel;
        else if (channel.frequency == SecondaryFrequency(channel.gnss_const))
          secondary_by_sat[key] = &channel;
      }
      std::vector<GnssChannel> out;
      for (const auto& [key, prim] : primary_by_sat) {
        auto it = secondary_by_sat.find(key);
        if (it == secondary_by_sat.end()) continue;
        const GnssChannel& sec = *it->second;
        const double f1 = (C / prim->Wavelength()).val();
        const double f2 = (C / sec.Wavelength()).val();
        const double denom = f1 * f1 - f2 * f2;
        const double a1 = f1 * f1 / denom;
        const double a2 = f2 * f2 / denom;
        // Secondary delay rescaled from the ray-traced primary delay (first-order: ~1/f^2).
        const double d_1 = prim->ionosphere_plasma_delay_m.val();
        const double d_2 = d_1 * (f1 * f1) / (f2 * f2);
        GnssChannel ifc = *prim;
        ifc.ionosphere_plasma_delay_m = a1 * d_1 - a2 * d_2;
        ifc.sigma_pseudorange_m
            = std::hypot(a1 * prim->sigma_pseudorange_m.val(), a2 * sec.sigma_pseudorange_m.val());
        out.push_back(ifc);
      }
      std::sort(out.begin(), out.end(), [](const GnssChannel& a, const GnssChannel& b) {
        if (a.gnss_const != b.gnss_const) return a.gnss_const < b.gnss_const;
        return a.prn < b.prn;
      });
      return out;
    }

    // Add the filter-only pseudorange noise inflation (quadrature) to a channel set.
    void InflatePseudorangeSigma(std::vector<GnssChannel>& channels, double inflation_m) {
      if (inflation_m <= 0.0) return;
      for (auto& channel : channels)
        channel.sigma_pseudorange_m = std::hypot(channel.sigma_pseudorange_m.val(), inflation_m);
    }

    double TransmitterBoresightAngleDeg(Real t_tdb, const Vec3& rx_mci,
                                        const GnssChannel& channel) {
      Vec3 rx_eci = ConvertFrame(t_tdb, rx_mci, Frame::MOON_CI, Frame::ECI);
      Vec3 tx_eci = ConvertFrame(t_tdb, Vec3(channel.tx_state.head(3)), channel.frame, Frame::ECI);
      Vec3 u_tx2rx = (rx_eci - tx_eci).normalized();
      Vec3 u_tx2earth = (-tx_eci).normalized();
      return (safe_acos(u_tx2rx.dot(u_tx2earth)) * 180.0 / PI).val();
    }

    bool LunarGnssLinkVisible(Real t_tdb, const Vec3& rx_mci, const GnssChannel& channel) {
      Vec3 rx_eci = ConvertFrame(t_tdb, rx_mci, Frame::MOON_CI, Frame::ECI);
      Vec3 tx_eci = ConvertFrame(t_tdb, Vec3(channel.tx_state.head(3)), channel.frame, Frame::ECI);
      if (!ComputeVisibility(rx_eci, tx_eci, R_EARTH, Vec3::Zero())) return false;

      Vec3 tx_mci = ConvertFrame(t_tdb, tx_eci, Frame::ECI, Frame::MOON_CI);
      return ComputeVisibility(rx_mci, tx_mci, R_MOON, Vec3::Zero());
    }

    GnssChannel ConvertChannelToMoonCi(const GnssChannel& channel) {
      GnssChannel out = channel;
      Real t_tx_tdb = ConvertTime(channel.transmit_time, channel.transmit_time_scale, Time::TDB);
      out.tx_state = ConvertFrame(t_tx_tdb, channel.tx_state, channel.frame, Frame::MOON_CI);
      out.frame = Frame::MOON_CI;
      out.ephemeris_chebyshev = ChebyshevFitModel{};
      if (channel.ephemeris_times.size() > 0 && channel.ephemeris_tx_states.rows() > 0) {
        out.ephemeris_tx_states.resize(channel.ephemeris_tx_states.rows(),
                                       channel.ephemeris_tx_states.cols());
        for (int i = 0; i < channel.ephemeris_tx_states.rows(); ++i) {
          Real t_tdb
              = ConvertTime(channel.ephemeris_times(i), channel.ephemeris_time_scale, Time::TDB);
          Vec6 rv_in = channel.ephemeris_tx_states.row(i).transpose();
          Vec6 rv_out = ConvertFrame(t_tdb, rv_in, channel.frame, Frame::MOON_CI);
          out.ephemeris_tx_states.row(i) = rv_out.transpose();
        }
      }
      return out;
    }

    std::vector<GNSSMeasurementsEpoch> PrecomputeConstellationChannels(
        const LunarGnssODTSConfig& cfg, const std::vector<RuntimeConstellation>& constellations,
        const std::vector<Real>& receive_times, const std::vector<State>& truth_states) {
      std::vector<GNSSMeasurementsEpoch> combined(receive_times.size());
      for (size_t i = 0; i < receive_times.size(); ++i) {
        combined[i].receive_time = receive_times[i];
        combined[i].time = receive_times[i];
        combined[i].receive_time_scale = Time::TDB;
      }

      // Parallelize across constellations (GPS L1/L5, Galileo E1/E5a): each is fully
      // independent (its own GNSSMeasurements + channel set), so the per-constellation C/N0
      // acquisition/tracking hysteresis stays sequential within each thread and the results
      // are merged deterministically in constellation order -- byte-identical to the serial
      // build. Frame/ephemeris (SPICE) calls in the hot path are guarded by `omp critical`
      // in the conversion layer, so this is thread-safe. Set precompute_num_threads>0 to cap
      // the thread count (0 = OpenMP default, i.e. all cores; capped to the constellation
      // count either way).
      const int n_const = static_cast<int>(constellations.size());
      std::vector<std::vector<GNSSMeasurementsEpoch>> per_const(n_const);
      const int n_threads = cfg.precompute_num_threads > 0
                                ? std::min(cfg.precompute_num_threads, n_const)
                                : n_const;
      ProgressReporter prog("Precompute GNSS links",
                            static_cast<long>(n_const) * static_cast<long>(receive_times.size()),
                            cfg.precompute_progress_interval_s);
      std::atomic<long> done{0};

#pragma omp parallel for schedule(dynamic) num_threads(n_threads) if (n_const > 1)
      for (int c = 0; c < n_const; ++c) {
        const auto& runtime = constellations[c];
        GNSSMeasurements gnss(runtime.constellation);
        gnss.SetFrequency(runtime.frequency);
        GnssMeasurementOptions candidate_options = MeasurementOptions(cfg, true);
        candidate_options.observables = PrecomputeObservables();
        candidate_options.frame = Frame::ECI;
        candidate_options.apply_visibility = false;
        gnss.SetOptions(candidate_options);
        gnss.SetReceiverParams(cfg.design.receiver_params);
        if (!cfg.design.receiver_antenna_name.empty())
          gnss.SetReceiverAntenna(Antenna(cfg.design.receiver_antenna_name));
        gnss.SetSunPositionProvider(
            [](Real t_tdb) { return GetBodyPos(t_tdb, BodyId::EARTH, BodyId::SUN, Frame::ECI); });
        if (cfg.plasma.simulate_truth)
          gnss.SetIonospherePlasmaRayTraceOptions(RayTraceOptions(cfg));
        std::vector<GNSSMeasurementsEpoch> epochs;
        epochs.reserve(receive_times.size());
        for (size_t i = 0; i < receive_times.size(); ++i) {
          gnss.SetOccludingBodies(LunarGnssOccludingBodies(receive_times[i]));
          State receiver_eci = truth_states[i];
          receiver_eci.head(6) = ConvertFrame(receive_times[i], Vec6(truth_states[i].head(6)),
                                              Frame::MOON_CI, Frame::ECI);
          GNSSMeasurementsEpoch epoch = gnss.Compute(receive_times[i], receiver_eci, nullptr);
          Vec3 rx_mci = truth_states[i].head(3);
          epoch.channels.erase(std::remove_if(epoch.channels.begin(), epoch.channels.end(),
                                              [&](const GnssChannel& channel) {
                                                return !LunarGnssLinkVisible(receive_times[i],
                                                                             rx_mci, channel);
                                              }),
                               epoch.channels.end());
          for (auto& channel : epoch.channels) channel = ConvertChannelToMoonCi(channel);
          epochs.push_back(epoch);
          const long d = ++done;
          if ((i & 0xFF) == 0 || i + 1 == receive_times.size()) {
#pragma omp critical(precompute_progress)
            prog.Update(d);
          }
        }
        per_const[c] = std::move(epochs);
      }
      prog.Update(static_cast<long>(n_const) * static_cast<long>(receive_times.size()));

      // Deterministic merge: constellation order, then epoch order (matches the serial build).
      for (int c = 0; c < n_const; ++c) {
        for (size_t i = 0; i < combined.size(); ++i) {
          combined[i].channels.insert(combined[i].channels.end(), per_const[c][i].channels.begin(),
                                      per_const[c][i].channels.end());
        }
      }
      return combined;
    }

    std::unordered_map<LinkKey, DelayRecord, LinkKeyHash> LoadDelayFile(
        const std::filesystem::path& path) {
      std::unordered_map<LinkKey, DelayRecord, LinkKeyHash> delays;
      if (path.empty() || !std::filesystem::exists(path)) return delays;

      std::ifstream in(path);
      std::string header;
      if (!std::getline(in, header)) return delays;
      std::unordered_map<std::string, int> idx = CsvColumnIndex(header);

      std::string line;
      while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row = SplitCsvLine(line);
        LinkKey key;
        key.epoch_index = std::stoi(CsvColumn(row, idx, "epoch_index"));
        key.gnss_const = std::stoi(CsvColumn(row, idx, "gnss_const_id"));
        key.prn = std::stoi(CsvColumn(row, idx, "prn"));
        key.frequency = std::stoi(CsvColumn(row, idx, "frequency_id"));
        DelayRecord record;
        std::string delay = CsvColumn(row, idx, "ionosphere_plasma_delay_m");
        if (delay.empty()) delay = CsvColumn(row, idx, "total_delay_m");
        record.ionosphere_plasma_delay_m = delay.empty() ? 0.0 : std::stod(delay);
        delays[key] = record;
      }
      return delays;
    }

    std::vector<GNSSMeasurementsEpoch> LoadLinksFile(const LunarGnssODTSConfig& cfg,
                                                     const VecXd& times_tdb) {
      std::ifstream in(cfg.links_file);
      LUPNT_CHECK(static_cast<bool>(in),
                  "Unable to open GNSS link cache: " + cfg.links_file.string(), "LunarGnssODTS");

      std::string header;
      LUPNT_CHECK(static_cast<bool>(std::getline(in, header)),
                  "GNSS link cache is empty: " + cfg.links_file.string(), "LunarGnssODTS");
      std::unordered_map<std::string, int> idx = CsvColumnIndex(header);
      RequireCsvColumns(idx,
                        {"epoch_index",
                         "t_tdb",
                         "gnss_const_id",
                         "prn",
                         "frequency_id",
                         "tx_x_mci_m",
                         "tx_y_mci_m",
                         "tx_z_mci_m",
                         "tx_vx_mci_mps",
                         "tx_vy_mci_mps",
                         "tx_vz_mci_mps",
                         "transmit_time_tdb",
                         "tx_clock_bias_s",
                         "tx_clock_drift_sps",
                         "relativistic_correction_s",
                         "group_delay_s",
                         "shapiro_delay_m",
                         "cn0_dbhz",
                         "sigma_pseudorange_m",
                         "sigma_doppler_hz",
                         "sigma_carrier_phase_cycles"},
                        cfg.links_file);

      std::vector<GNSSMeasurementsEpoch> epochs(times_tdb.size());
      for (int i = 0; i < times_tdb.size(); ++i) {
        epochs[i].receive_time = times_tdb(i);
        epochs[i].time = times_tdb(i);
        epochs[i].receive_time_scale = Time::TDB;
      }

      int n_links = 0;
      std::string line;
      while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row = SplitCsvLine(line);
        int epoch_index = std::stoi(CsvColumn(row, idx, "epoch_index"));
        LUPNT_CHECK(epoch_index >= 0 && epoch_index < static_cast<int>(epochs.size()),
                    "GNSS link cache epoch_index is outside the configured run duration: "
                        + std::to_string(epoch_index),
                    "LunarGnssODTS");

        GnssChannel channel;
        channel.gnss_const
            = static_cast<GnssConst>(std::stoi(CsvColumn(row, idx, "gnss_const_id")));
        channel.prn = std::stoi(CsvColumn(row, idx, "prn"));
        channel.frequency = static_cast<GnssFreq>(std::stoi(CsvColumn(row, idx, "frequency_id")));
        channel.receive_time = std::stod(CsvColumn(row, idx, "t_tdb"));
        channel.receive_time_scale = Time::TDB;
        channel.transmit_time = std::stod(CsvColumn(row, idx, "transmit_time_tdb"));
        channel.transmit_time_scale = Time::TDB;
        channel.ephemeris_time_scale = Time::TDB;
        channel.frame = Frame::MOON_CI;
        channel.tx_state << std::stod(CsvColumn(row, idx, "tx_x_mci_m")),
            std::stod(CsvColumn(row, idx, "tx_y_mci_m")),
            std::stod(CsvColumn(row, idx, "tx_z_mci_m")),
            std::stod(CsvColumn(row, idx, "tx_vx_mci_mps")),
            std::stod(CsvColumn(row, idx, "tx_vy_mci_mps")),
            std::stod(CsvColumn(row, idx, "tx_vz_mci_mps"));
        channel.tx_clock_bias_s = std::stod(CsvColumn(row, idx, "tx_clock_bias_s"));
        channel.tx_clock_drift = std::stod(CsvColumn(row, idx, "tx_clock_drift_sps"));
        channel.relativistic_correction_s
            = std::stod(CsvColumn(row, idx, "relativistic_correction_s"));
        channel.group_delay_s = std::stod(CsvColumn(row, idx, "group_delay_s"));
        channel.shapiro_delay_m = std::stod(CsvColumn(row, idx, "shapiro_delay_m"));
        channel.cn0_dbhz = std::stod(CsvColumn(row, idx, "cn0_dbhz"));
        channel.sigma_pseudorange_m = std::stod(CsvColumn(row, idx, "sigma_pseudorange_m"));
        channel.sigma_doppler_hz = std::stod(CsvColumn(row, idx, "sigma_doppler_hz"));
        channel.sigma_carrier_phase_cycles
            = std::stod(CsvColumn(row, idx, "sigma_carrier_phase_cycles"));
        epochs[epoch_index].channels.push_back(channel);
        ++n_links;
      }

      LUPNT_CHECK(n_links > 0,
                  "GNSS link cache contains no usable rows: " + cfg.links_file.string(),
                  "LunarGnssODTS");
      Logger::Info("Loaded cached GNSS links from " + cfg.links_file.string() + " ("
                       + std::to_string(n_links) + " links)",
                   "LunarGnssODTS");
      return epochs;
    }

    void ApplyDelayFile(std::vector<GNSSMeasurementsEpoch>& epochs,
                        const std::filesystem::path& delay_file) {
      std::unordered_map<LinkKey, DelayRecord, LinkKeyHash> delays = LoadDelayFile(delay_file);
      LUPNT_CHECK(!delays.empty(),
                  "Delay table contains no usable GNSS delay rows: " + delay_file.string(),
                  "LunarGnssODTS");

      for (int epoch_index = 0; epoch_index < static_cast<int>(epochs.size()); ++epoch_index) {
        for (auto& channel : epochs[epoch_index].channels) {
          LinkKey key{epoch_index, static_cast<int>(channel.gnss_const), channel.prn,
                      static_cast<int>(channel.frequency)};
          auto it = delays.find(key);
          if (it != delays.end()) {
            channel.ionosphere_plasma_delay_m = it->second.ionosphere_plasma_delay_m;
          }
        }
      }
    }

    void WriteLinksCsv(const LunarGnssODTSConfig& cfg,
                       const std::vector<GNSSMeasurementsEpoch>& epochs,
                       const std::vector<State>& truth_states, int epoch_index_offset = 0) {
      std::filesystem::create_directories(cfg.links_file.parent_path());
      std::ofstream out(cfg.links_file);
      out << std::scientific << std::setprecision(12);
      out << "link_id,epoch_index,t_tdb,t_utc,gnss_const_id,gnss_const,prn,frequency_id,frequency,"
             "rx_x_ecef_m,rx_y_ecef_m,rx_z_ecef_m,tx_x_ecef_m,tx_y_ecef_m,tx_z_ecef_m,"
             "tx_x_mci_m,tx_y_mci_m,tx_z_mci_m,tx_vx_mci_mps,tx_vy_mci_mps,tx_vz_mci_mps,"
             "transmit_time_tdb,tx_clock_bias_s,tx_clock_drift_sps,relativistic_correction_s,"
             "group_delay_s,shapiro_delay_m,sigma_pseudorange_m,sigma_doppler_hz,"
             "sigma_carrier_phase_cycles,"
             "tangent_radius_earth_center_m,tangent_altitude_m,tx_boresight_angle_deg,cn0_dbhz,"
             "ionosphere_plasma_delay_m,tecu,tec_delay_m,second_delay_m,third_delay_m,"
             "dist_bend_m,tec_delay_bend_m,max_sep_line_m,final_pos_err_m\n";

      int link_id = 0;
      for (int k = 0; k < static_cast<int>(epochs.size()); ++k) {
        const Real t_tdb = epochs[k].receive_time;
        const Real t_utc = ConvertTime(t_tdb, Time::TDB, Time::UTC);
        Vec3 rx_mci = truth_states[k].head(3);
        Vec3 rx_eci = ConvertFrame(t_tdb, rx_mci, Frame::MOON_CI, Frame::ECI);
        Vec3 rx_ecef = ConvertFrame(t_tdb, rx_eci, Frame::ECI, Frame::ECEF);
        for (const auto& channel : epochs[k].channels) {
          Vec3 tx_eci
              = ConvertFrame(t_tdb, Vec3(channel.tx_state.head(3)), channel.frame, Frame::ECI);
          Vec3 tx_ecef = ConvertFrame(t_tdb, tx_eci, Frame::ECI, Frame::ECEF);
          Real transmit_time_tdb
              = ConvertTime(channel.transmit_time, channel.transmit_time_scale, Time::TDB);
          auto [tangent_radius_m, tangent_altitude_m]
              = EarthTangentRadiusAltitude(rx_ecef, tx_ecef);
          double tx_boresight_angle_deg = TransmitterBoresightAngleDeg(t_tdb, rx_mci, channel);
          out << link_id++ << "," << (epoch_index_offset + k) << "," << t_tdb.val() << ","
              << t_utc.val() << "," << static_cast<int>(channel.gnss_const) << ","
              << GnssConstName(channel.gnss_const) << "," << channel.prn << ","
              << static_cast<int>(channel.frequency) << "," << GnssFreqName(channel.frequency)
              << "," << rx_ecef(0).val() << "," << rx_ecef(1).val() << "," << rx_ecef(2).val()
              << "," << tx_ecef(0).val() << "," << tx_ecef(1).val() << "," << tx_ecef(2).val()
              << "," << channel.tx_state(0).val() << "," << channel.tx_state(1).val() << ","
              << channel.tx_state(2).val() << "," << channel.tx_state(3).val() << ","
              << channel.tx_state(4).val() << "," << channel.tx_state(5).val() << ","
              << transmit_time_tdb.val() << "," << channel.tx_clock_bias_s.val() << ","
              << channel.tx_clock_drift.val() << "," << channel.relativistic_correction_s.val()
              << "," << channel.group_delay_s.val() << "," << channel.shapiro_delay_m.val() << ","
              << channel.sigma_pseudorange_m.val() << "," << channel.sigma_doppler_hz.val() << ","
              << channel.sigma_carrier_phase_cycles.val() << "," << tangent_radius_m << ","
              << tangent_altitude_m << "," << tx_boresight_angle_deg << ","
              << channel.cn0_dbhz.val() << ",0,0,0,0,0,0,0,0,0\n";
        }
      }
    }

    struct LunarGnssODTSRuntimeContext {
      const LunarGnssODTSConfig& cfg;
      int mc_index = 0;
      Real t0_tdb = 0.0;
      const VecXd& elapsed_s;
      const VecXd& times_tdb;
      const VecXd& receiver_clock_s;
      const std::vector<State>& truth_states;
      const std::vector<GNSSMeasurementsEpoch>& precomputed;
    };

    class LunarODTSApp : public LunaNetSubApp {
    public:
      explicit LunarODTSApp(const LunarGnssODTSRuntimeContext& context)
          : LunaNetSubApp("lunar_odts"), context_(context) {}

      void Setup(LunaNetSatApp& app) override {
        LunaNetSubApp::Setup(app);
        RandomEngine::SetSeed(static_cast<unsigned int>(context_.cfg.seed + context_.mc_index));

        truth_options_ = MeasurementOptions(context_.cfg, true);
        filter_options_ = MeasurementOptions(context_.cfg, false);
        truth_carrier_options_ = CarrierMeasurementOptions(context_.cfg, true);
        filter_carrier_options_ = CarrierMeasurementOptions(context_.cfg, false);

        filter_ = context_.cfg.use_tdcp
                      ? std::static_pointer_cast<UDUEKF>(MakePtr<UDUStochasticCloningEKF>())
                      : MakePtr<UDUEKF>();

        State x0 = MakeInitialEstimateState(context_.truth_states.front(), context_.cfg);
        MatXd P0 = InitialCovariance(context_.cfg);
        filter_->SetTime(context_.times_tdb(0));
        if (context_.cfg.use_tdcp) {
          auto cloned_filter = std::dynamic_pointer_cast<UDUStochasticCloningEKF>(filter_);
          cloned_filter->SetBaseStateSize(x0.size());
          filter_->SetState(MakeClonedState(x0));
          filter_->SetCovariance(MakeClonedCovariance(P0));
        } else {
          filter_->SetState(x0);
          filter_->SetCovariance(P0);
        }

        Ptr<NBodyDynamics> orbit_dynamics_filter = CreateLunarDynamics(
            context_.cfg.moon_gravity_degree_filter, context_.cfg.moon_gravity_order_filter, true,
            context_.cfg, UseFilterSrp(context_.cfg), context_.cfg.srp_coeff_filter_m2_kg);
        Ptr<JointOrbitClockDynamics> dynamics_filter
            = CreateJointDynamics(orbit_dynamics_filter, ClockModelFromConfig(context_.cfg), false,
                                  context_.cfg.seed + 2000 + context_.mc_index);
        filter_->SetDynamicsFunction(CreateDynamicsFunction(dynamics_filter, context_.cfg));
        filter_->SetProcessNoiseFunction(CreateProcessNoiseFunction(context_.cfg, filter_.get()));
        filter_->SetOutlierThreshold(1.0e12);

        std::filesystem::create_directories(context_.cfg.output_dir);
        trajectory_.open(context_.cfg.output_dir
                         / ("trajectory_mc" + std::to_string(context_.mc_index) + ".csv"));
        trajectory_ << "mc,t,receiver_clock_s,pos_error_m,vel_error_mps,"
                       "x_error_m,y_error_m,z_error_m,vx_error_mps,vy_error_mps,vz_error_mps,"
                       "clock_bias_error_m,clock_drift_error_mps,"
                       "x_3sigma_m,y_3sigma_m,z_3sigma_m,vx_3sigma_mps,vy_3sigma_mps,vz_3sigma_mps,"
                       "r_error_m,t_error_m,n_error_m,rdot_error_mps,tdot_error_mps,ndot_error_mps,"
                       "r_3sigma_m,t_3sigma_m,n_3sigma_m,rdot_3sigma_mps,tdot_3sigma_mps,"
                       "ndot_3sigma_mps,"
                       "clock_bias_3sigma_m,clock_drift_3sigma_mps,"
                       "srp_coeff_est_m2_kg,srp_coeff_error_m2_kg,srp_coeff_3sigma_m2_kg,"
                       "num_channels,num_tracked_satellites,num_measurements\n";

        // Per-satellite, per-epoch range-direction transmitter ephemeris (broadcast-vs-precise)
        // residuals for the whole run, for distribution analysis (histogram). Deterministic, so
        // identical across Monte-Carlo runs.
        eph_residuals_.open(
            context_.cfg.output_dir
            / ("ephemeris_residuals_mc" + std::to_string(context_.mc_index) + ".csv"));
        eph_residuals_ << "mc,epoch,t_s,gnss,prn,eph_los_m,eph_clk_m\n";

        summary_.monte_carlo_index = context_.mc_index;
        summary_.num_epochs = static_cast<int>(context_.times_tdb.size());
        pos_err2_sum_ = 0.0;
        vel_err2_sum_ = 0.0;
        next_epoch_index_ = 0;

        // Progress reporting: configured epoch stride, or ~50 update lines over the run when
        // run_progress_interval_epochs <= 0. Always includes the first/last epoch.
        const int n_epochs = static_cast<int>(context_.times_tdb.size());
        progress_stride_ = context_.cfg.run_progress_interval_epochs > 0
                               ? context_.cfg.run_progress_interval_epochs
                               : std::max(1, n_epochs / 50);
        diverged_warned_ = false;
        start_time_ = std::chrono::steady_clock::now();
        Logger::Info("Running EKF over " + std::to_string(n_epochs) + " epochs (mc "
                         + std::to_string(context_.mc_index) + ")",
                     "LunarGnssODTS");

        // Broadcast-ephemeris injection: build the live broadcast (RINEX-nav) transmitter-error
        // model once. Truth keeps the precise SP3 transmitter states; the filter measurement
        // model is fed the debiased broadcast position/clock (see BroadcastEphemerisError).
        broadcast_error_.reset();
        if (context_.cfg.constellation.use_broadcast_ephemeris) {
          std::set<std::pair<GnssConst, int>> sat_set;
          for (const auto& epoch : context_.precomputed)
            for (const auto& ch : epoch.channels) sat_set.emplace(ch.gnss_const, ch.prn);
          const std::vector<std::pair<GnssConst, int>> sats(sat_set.begin(), sat_set.end());
          const std::vector<std::filesystem::path> brdc_files
              = ResolveBrdcFiles(context_.cfg.constellation);
          LUPNT_CHECK(!brdc_files.empty(),
                      "use_broadcast_ephemeris is set but no BRDC (RINEX-nav) files were found; "
                      "set constellation.brdc_directory (or brdc_files)",
                      "LunarGnssODTS");
          const auto [t_start_tai, t_end_tai] = EphemerisWindowTai(context_.cfg);
          broadcast_error_ = std::make_shared<BroadcastEphemerisError>(
              context_.cfg.constellation.sp3_files, brdc_files,
              context_.cfg.constellation.antex_file, sats, t_start_tai, t_end_tai,
              /*sample_dt_s=*/300.0, context_.cfg.constellation.debias_broadcast_clock,
              context_.cfg.constellation.debias_qzss_radial);
          Logger::Info("Broadcast-ephemeris injection enabled (" + std::to_string(sats.size())
                           + " satellites, " + std::to_string(brdc_files.size()) + " BRDC files)",
                       "LunarGnssODTS");
        }
      }

      void Step(Real) override {
        const int k = next_epoch_index_++;
        LUPNT_CHECK(k < context_.times_tdb.size(), "Lunar ODTS app called after final epoch",
                    "LunarGnssODTS");
        std::vector<GnssChannel> raw_truth_channels = context_.precomputed[k].channels;
        ApplyMeasurementSigmas(raw_truth_channels, context_.cfg, false);

        const Real t_tdb = context_.times_tdb(k);
        const Vec3 rx_mci = context_.truth_states[k].head(3);

        // Pseudorange channels: ionosphere-free (L1+L5) combination when enabled, else
        // primary-frequency only; then apply the pseudorange tangent-altitude cutoff.
        std::vector<GnssChannel> truth_channels
            = context_.cfg.use_ionosphere_free ? MakeIonosphereFreeChannels(raw_truth_channels)
                                               : SelectPrimaryFrequency(raw_truth_channels);
        truth_channels = FilterByTangentAltitude(truth_channels, t_tdb, rx_mci,
                                                 context_.cfg.pseudorange_min_tangent_altitude_m);
        std::vector<GnssChannel> filter_channels = MakeFilterChannels(truth_channels, context_.cfg);
        if (context_.cfg.use_ionosphere_free) {
          // IF sigma is already the C/N0 combination; add the filter-only inflation here.
          InflatePseudorangeSigma(filter_channels,
                                  context_.cfg.filter_pseudorange_noise_inflation_m);
        } else {
          ApplyMeasurementSigmas(filter_channels, context_.cfg, true);
        }

        // TDCP channels use the primary frequency (GPS L1, Galileo E1), with their own
        // tangent-altitude cutoff. These retain the ray-traced primary-frequency plasmaspheric
        // delay, so truth TDCP carrier phase is simulated with the real delay change even when
        // ionosphere-free pseudorange cancels the first-order delay. Carrier-phase sigmas stay
        // at the C/N0 floor; filter-only TDCP inflation is added in TdcpCovariance.
        std::vector<GnssChannel> truth_primary = SelectPrimaryFrequency(raw_truth_channels);
        truth_primary = FilterByTangentAltitude(truth_primary, t_tdb, rx_mci,
                                                context_.cfg.tdcp_min_tangent_altitude_m);
        std::vector<GnssChannel> filter_primary = MakeFilterChannels(truth_primary, context_.cfg);

        // Feed the filter (receiver) measurement model the broadcast transmitter ephemeris,
        // evaluated live from the navigation-message parameters at each signal transmit epoch;
        // the truth channels keep the precise SP3 states. The debiased broadcast-minus-precise
        // error (per satellite) is thus injected as an unmodeled measurement error. Cache per
        // satellite within the epoch so the pseudorange and TDCP channel sets share one eval.
        if (broadcast_error_) {
          std::map<std::pair<GnssConst, int>, std::pair<Vec3, Real>> cache;
          auto inject = [&](std::vector<GnssChannel>& chans) {
            for (auto& ch : chans) {
              const std::pair<GnssConst, int> key{ch.gnss_const, ch.prn};
              auto it = cache.find(key);
              if (it == cache.end()) {
                Vec3 dr;
                Real dc;
                const Real t_tx_tai
                    = ConvertTime(ch.transmit_time, ch.transmit_time_scale, Time::TAI);
                if (!broadcast_error_->GetDebiasedDelta(ch.gnss_const, ch.prn, ch.frequency,
                                                        t_tx_tai, dr, dc))
                  continue;  // no broadcast message -> leave this satellite on the precise state
                it = cache.emplace(key, std::make_pair(dr, dc)).first;
              }
              ch.tx_state.head(3) += it->second.first;
              ch.tx_clock_bias_s += it->second.second;
            }
          };
          inject(filter_channels);
          inject(filter_primary);
        }

        if (k > 0) filter_->Predict(context_.times_tdb(k));

        VecXd y_true = LunarGnssCombinedMeasurement::ComputeMeasurementVector(
            context_.truth_states[k], truth_channels, truth_options_);
        VecXd y_obs = AddMeasurementNoise(y_true, truth_channels, truth_options_);
        std::vector<TdcpPair> truth_tdcp_pairs;
        std::vector<TdcpPair> filter_tdcp_pairs;
        if (context_.cfg.use_tdcp && k > 0) {
          truth_tdcp_pairs
              = LunarGnssCombinedMeasurement::MakeTdcpPairs(truth_primary, previous_truth_primary_);
          filter_tdcp_pairs = LunarGnssCombinedMeasurement::MakeTdcpPairs(filter_primary,
                                                                          previous_filter_primary_);
          if (!truth_tdcp_pairs.empty() && !filter_tdcp_pairs.empty()) {
            VecXd y_tdcp_true = LunarGnssCombinedMeasurement::ComputeTdcpVector(
                context_.truth_states[k], context_.truth_states[k - 1], truth_tdcp_pairs,
                truth_carrier_options_);
            VecXd y_tdcp_obs = AddTdcpNoise(y_tdcp_true, truth_tdcp_pairs, context_.cfg);
            VecXd y_combined(y_obs.size() + y_tdcp_obs.size());
            y_combined << y_obs, y_tdcp_obs;
            y_obs = y_combined;
          }
        }

        LunarGnssCombinedMeasurement::Config meas_cfg;
        meas_cfg.channels = filter_channels;
        meas_cfg.current_options = filter_options_;
        meas_cfg.carrier_options = filter_carrier_options_;
        meas_cfg.tdcp_pairs = filter_tdcp_pairs;
        meas_cfg.use_tdcp = context_.cfg.use_tdcp;
        meas_cfg.tdcp_sigma_m = context_.cfg.tdcp_sigma_m;
        meas_cfg.filter_tdcp_noise_inflation_m = context_.cfg.filter_tdcp_noise_inflation_m;
        filter_->SetMeasurementFunction(LunarGnssCombinedMeasurement(meas_cfg).CreateFunction());
        filter_->Update(y_obs);
        PrintMatrixDiagnostics(k, context_.cfg, filter_, truth_channels, filter_channels,
                               context_.truth_states[k].head(3).cast<double>(), truth_options_,
                               truth_tdcp_pairs);

        // Persist the per-satellite range-direction ephemeris residuals every epoch (for the
        // distribution / histogram in the notebook).
        if (eph_residuals_.is_open()) {
          const Vec3d rx = context_.truth_states[k].head(3).cast<double>();
          for (std::size_t c = 0; c < truth_channels.size() && c < filter_channels.size(); ++c)
            eph_residuals_ << context_.mc_index << "," << k << "," << context_.elapsed_s(k) << ","
                           << GnssConstName(truth_channels[c].gnss_const) << ","
                           << truth_channels[c].prn << ","
                           << LosEphError(truth_channels[c], filter_channels[c], rx) << ","
                           << ClockEphError(truth_channels[c], filter_channels[c]) << "\n";
        }

        State x_est = CurrentFilterState(filter_->GetState(), context_.cfg);
        Vec3d dr = (context_.truth_states[k].head(3) - x_est.head(3)).cast<double>();
        Vec3d dv = (context_.truth_states[k].segment(3, 3) - x_est.segment(3, 3)).cast<double>();
        const double pos_err = dr.norm();
        const double vel_err = dv.norm();
        const double clk_b_err_m = (context_.truth_states[k](6) - x_est(6)).val();
        const double clk_d_err_mps = (context_.truth_states[k](7) - x_est(7)).val();
        const double srp_est
            = EstimateSrp(context_.cfg) ? x_est(8).val() : std::numeric_limits<double>::quiet_NaN();
        const double srp_err = EstimateSrp(context_.cfg)
                                   ? (x_est(8).val() - TruthSrpCoeff(context_.cfg))
                                   : std::numeric_limits<double>::quiet_NaN();

        MatXd P = CurrentFilterCovariance(filter_->GetCovariance(), context_.cfg);
        const double sx = 3.0 * std::sqrt(std::max(0.0, P(0, 0)));
        const double sy = 3.0 * std::sqrt(std::max(0.0, P(1, 1)));
        const double sz = 3.0 * std::sqrt(std::max(0.0, P(2, 2)));
        const double svx = 3.0 * std::sqrt(std::max(0.0, P(3, 3)));
        const double svy = 3.0 * std::sqrt(std::max(0.0, P(4, 4)));
        const double svz = 3.0 * std::sqrt(std::max(0.0, P(5, 5)));
        const double sb = 3.0 * std::sqrt(std::max(0.0, P(6, 6)));
        const double sd = 3.0 * std::sqrt(std::max(0.0, P(7, 7)));
        const double ssrp = EstimateSrp(context_.cfg) ? 3.0 * std::sqrt(std::max(0.0, P(8, 8)))
                                                      : std::numeric_limits<double>::quiet_NaN();
        const Mat3d R_rtn = RotCartToRtn(context_.truth_states[k].head(3).cast<double>(),
                                         context_.truth_states[k].segment(3, 3).cast<double>());
        const Vec3d dr_rtn = R_rtn.transpose() * dr;
        const Vec3d dv_rtn = R_rtn.transpose() * dv;
        const Mat3d Pr_rtn = R_rtn.transpose() * P.block(0, 0, 3, 3) * R_rtn;
        const Mat3d Pv_rtn = R_rtn.transpose() * P.block(3, 3, 3, 3) * R_rtn;
        const Vec3d srtn = 3.0 * Pr_rtn.diagonal().cwiseMax(0.0).cwiseSqrt();
        const Vec3d svrtn = 3.0 * Pv_rtn.diagonal().cwiseMax(0.0).cwiseSqrt();
        const int num_tracked_satellites = CountTrackedSatellites(filter_channels);

        pos_err2_sum_ += pos_err * pos_err;
        vel_err2_sum_ += vel_err * vel_err;

        trajectory_ << context_.mc_index << "," << context_.elapsed_s(k) << ","
                    << context_.receiver_clock_s(k) << "," << pos_err << "," << vel_err << ","
                    << dr(0) << "," << dr(1) << "," << dr(2) << "," << dv(0) << "," << dv(1) << ","
                    << dv(2) << "," << clk_b_err_m << "," << clk_d_err_mps << "," << sx << "," << sy
                    << "," << sz << "," << svx << "," << svy << "," << svz << "," << dr_rtn(0)
                    << "," << dr_rtn(1) << "," << dr_rtn(2) << "," << dv_rtn(0) << "," << dv_rtn(1)
                    << "," << dv_rtn(2) << "," << srtn(0) << "," << srtn(1) << "," << srtn(2) << ","
                    << svrtn(0) << "," << svrtn(1) << "," << svrtn(2) << "," << sb << "," << sd
                    << "," << srp_est << "," << srp_err << "," << ssrp << ","
                    << filter_channels.size() << "," << num_tracked_satellites << ","
                    << y_obs.size() << "\n";

        previous_truth_primary_ = truth_primary;
        previous_filter_primary_ = filter_primary;

        // --- Live progress + error reporting -----------------------------------------
        const int n_epochs = static_cast<int>(context_.times_tdb.size());
        const bool nonfinite = !std::isfinite(pos_err) || !std::isfinite(clk_b_err_m);
        const bool diverged = pos_err > 1.0e6;  // >1000 km: clearly diverged for an ELFO
        if ((nonfinite || diverged) && !diverged_warned_) {
          Logger::Warn("EKF diverging at epoch " + std::to_string(k + 1) + "/"
                           + std::to_string(n_epochs) + ": position error "
                           + (nonfinite ? std::string("non-finite (NaN/Inf)")
                                        : std::to_string(static_cast<long long>(pos_err)) + " m")
                           + " -- check process/measurement noise tuning",
                       "LunarGnssODTS");
          diverged_warned_ = true;  // warn once per run to avoid flooding the log
        }
        if ((k % progress_stride_ == 0) || (k == n_epochs - 1)) {
          const double wall_s
              = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time_)
                    .count();
          std::ostringstream msg;
          msg << std::fixed << std::setprecision(1) << "epoch " << (k + 1) << "/" << n_epochs
              << " (" << (100.0 * (k + 1) / n_epochs) << "%): pos " << pos_err << " m, clk "
              << clk_b_err_m << " m, " << num_tracked_satellites << " sat, " << y_obs.size()
              << " meas";
          if (EstimateSrp(context_.cfg))
            msg << ", SRP " << std::setprecision(3) << srp_est * 1e3 << "e-3";
          msg << " [" << std::setprecision(1) << wall_s << "s]";
          Logger::Info(msg.str(), "LunarGnssODTS");
        }

        if (k == n_epochs - 1) {
          summary_.final_position_error_m = pos_err;
          summary_.final_velocity_error_mps = vel_err;
          summary_.final_clock_bias_error_m = clk_b_err_m;
          summary_.final_clock_drift_error_mps = clk_d_err_mps;
          summary_.final_srp_coeff_error_m2_kg = srp_err;
        }
      }

      void Finish() override {
        summary_.rms_position_error_m
            = std::sqrt(pos_err2_sum_ / static_cast<double>(context_.times_tdb.size()));
        summary_.rms_velocity_error_mps
            = std::sqrt(vel_err2_sum_ / static_cast<double>(context_.times_tdb.size()));
        trajectory_.close();
        eph_residuals_.close();

        const double wall_s
            = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time_).count();
        std::ostringstream msg;
        msg << std::fixed << std::setprecision(2) << "EKF done (mc " << context_.mc_index << ") in "
            << std::setprecision(1) << wall_s << "s: final pos error " << std::setprecision(2)
            << summary_.final_position_error_m << " m, RMS " << summary_.rms_position_error_m
            << " m";
        Logger::Info(msg.str(), "LunarGnssODTS");
      }

      const LunarGnssODTSSummary& GetSummary() const { return summary_; }

    private:
      const LunarGnssODTSRuntimeContext& context_;
      GnssMeasurementOptions truth_options_;
      GnssMeasurementOptions filter_options_;
      GnssMeasurementOptions truth_carrier_options_;
      GnssMeasurementOptions filter_carrier_options_;
      Ptr<UDUEKF> filter_;
      std::ofstream trajectory_;
      std::ofstream eph_residuals_;
      std::vector<GnssChannel> previous_truth_primary_;
      std::vector<GnssChannel> previous_filter_primary_;
      std::shared_ptr<BroadcastEphemerisError> broadcast_error_;
      LunarGnssODTSSummary summary_;
      double pos_err2_sum_ = 0.0;
      double vel_err2_sum_ = 0.0;
      int next_epoch_index_ = 0;
      int progress_stride_ = 1;
      bool diverged_warned_ = false;
      std::chrono::steady_clock::time_point start_time_;
    };

    LunarGnssODTSSummary RunOneMonteCarlo(const LunarGnssODTSConfig& cfg, int mc_index, Real t0_tdb,
                                          const VecXd& elapsed_s, const VecXd& times_tdb,
                                          const VecXd& receiver_clock_s,
                                          std::vector<GNSSMeasurementsEpoch> precomputed,
                                          const std::vector<State>* receiver_truth = nullptr) {
      // Truth trajectory: use the physical receiver `Spacecraft`'s self-propagated truth when
      // the hosting app provides it (agent-driven path); otherwise build it internally from the
      // config (legacy struct API / standalone). Both land on the same absolute TDB epochs.
      std::vector<State> truth_states = receiver_truth
                                            ? *receiver_truth
                                            : BuildTruthStates(cfg, t0_tdb, elapsed_s, times_tdb,
                                                               true, cfg.seed + 1000 + mc_index);

      LunarGnssODTSRuntimeContext context{cfg,       mc_index,         t0_tdb,       elapsed_s,
                                          times_tdb, receiver_clock_s, truth_states, precomputed};
      Ptr<LunarODTSApp> odts_app = MakePtr<LunarODTSApp>(context);
      LunarODTSApp* odts_app_ptr = odts_app.get();

      LunaNetSatApp sat_app;
      sat_app.SetName("lunanet_sat");
      sat_app.AddSubApp(odts_app);
      sat_app.Setup();
      for (int k = 0; k < times_tdb.size(); ++k) sat_app.Step(elapsed_s(k));
      sat_app.Finish();
      return odts_app_ptr->GetSummary();
    }

    void WriteSummary(const LunarGnssODTSConfig& cfg,
                      const std::vector<LunarGnssODTSSummary>& summaries) {
      std::filesystem::create_directories(cfg.output_dir);
      std::ofstream out(cfg.output_dir / "summary.csv");
      out << "mc,epochs,final_position_error_m,final_velocity_error_mps,"
             "final_clock_bias_error_m,final_clock_drift_error_mps,"
             "final_srp_coeff_error_m2_kg,rms_position_error_m,rms_velocity_error_mps\n";
      out << std::scientific << std::setprecision(12);
      for (const auto& s : summaries) {
        out << s.monte_carlo_index << "," << s.num_epochs << "," << s.final_position_error_m << ","
            << s.final_velocity_error_mps << "," << s.final_clock_bias_error_m << ","
            << s.final_clock_drift_error_mps << "," << s.final_srp_coeff_error_m2_kg << ","
            << s.rms_position_error_m << "," << s.rms_velocity_error_mps << "\n";
      }
    }

    void ClearPreviousOutputs(const LunarGnssODTSConfig& cfg) {
      if (!std::filesystem::exists(cfg.output_dir)) return;
      for (const auto& entry : std::filesystem::directory_iterator(cfg.output_dir)) {
        if (!entry.is_regular_file()) continue;
        const std::string name = entry.path().filename().string();
        if ((name.rfind("trajectory_mc", 0) == 0 && entry.path().extension() == ".csv")
            || name == "summary.csv") {
          std::filesystem::remove(entry.path());
        }
      }
    }

    void ReadReceiverParams(const YAML::Node& node, GnssReceiverParams& params) {
      params.Bp = ReadYaml(node, "carrier_loop_bandwidth_hz", params.Bp);
      params.T = ReadYaml(node, "integration_time_s", params.T);
      params.b = ReadYaml(node, "front_end_bandwidth_factor", params.b);
      params.Bn = ReadYaml(node, "code_loop_bandwidth_hz", params.Bn);
      params.Bf = ReadYaml(node, "frequency_loop_bandwidth_hz", params.Bf);
      params.D = ReadYaml(node, "early_late_spacing_chips", params.D);
      params.L_ad = ReadYaml(node, "ad_converter_loss_db", params.L_ad);
      params.L_pol = ReadYaml(node, "polarization_loss_db", params.L_pol);
      params.L_atm = ReadYaml(node, "atmospheric_loss_db", params.L_atm);
      params.T_eff = ReadYaml(node, "effective_noise_temperature_k", params.T_eff);
    }

    void ApplyDesignNode(const YAML::Node& node, LunarGnssODTSConfig& cfg) {
      if (!node) return;
      ReadReceiverParams(node["receiver_tracking"], cfg.design.receiver_params);

      const YAML::Node link = node["link_budget"];
      cfg.design.cn0_threshold_dbhz
          = ReadYaml(link, "cn0_threshold_dbhz", cfg.design.cn0_threshold_dbhz);
      cfg.design.cn0_acquisition_threshold_dbhz = ReadYaml(
          link, "cn0_acquisition_threshold_dbhz", cfg.design.cn0_acquisition_threshold_dbhz);
      cfg.design.cn0_tracking_threshold_dbhz
          = ReadYaml(link, "cn0_tracking_threshold_dbhz", cfg.design.cn0_tracking_threshold_dbhz);
      cfg.design.apply_cn0_threshold
          = ReadYaml(link, "apply_cn0_threshold", cfg.design.apply_cn0_threshold);
      cfg.design.setup_transmitters
          = ReadYaml(link, "setup_transmitters", cfg.design.setup_transmitters);
      cfg.design.receiver_antenna_name
          = ReadYaml(link, "receiver_antenna_name", cfg.design.receiver_antenna_name);
      cfg.design.receiver_antenna_name
          = ReadYaml(link, "receiver_antenna", cfg.design.receiver_antenna_name);
      cfg.design.use_cn0_measurement_sigmas
          = ReadYaml(link, "use_cn0_measurement_sigmas", cfg.design.use_cn0_measurement_sigmas);
      cfg.design.tx_yaw_dedicated = ReadYaml(link, "tx_yaw_dedicated", cfg.design.tx_yaw_dedicated);

      const YAML::Node meas = node["measurement_noise"];
      cfg.pseudorange_sigma_m = ReadYaml(meas, "pseudorange_sigma_m", cfg.pseudorange_sigma_m);
      cfg.doppler_sigma_hz = ReadYaml(meas, "doppler_sigma_hz", cfg.doppler_sigma_hz);
      cfg.carrier_phase_sigma_m
          = ReadYaml(meas, "carrier_phase_sigma_m", cfg.carrier_phase_sigma_m);
      cfg.tdcp_sigma_m = ReadYaml(meas, "tdcp_sigma_m", cfg.tdcp_sigma_m);

      const YAML::Node plasma = node["plasma_mismatch"];
      cfg.plasma.filter_pseudorange_noise_inflation_m
          = ReadYaml(plasma, "filter_pseudorange_noise_inflation_m",
                     cfg.plasma.filter_pseudorange_noise_inflation_m);
      cfg.plasma.filter_doppler_noise_inflation_hz
          = ReadYaml(plasma, "filter_doppler_noise_inflation_hz",
                     cfg.plasma.filter_doppler_noise_inflation_hz);
    }

    void LoadDesignDatabase(LunarGnssODTSConfig& cfg) {
      if (cfg.design.database_path.empty() || !std::filesystem::exists(cfg.design.database_path)) {
        return;
      }

      YAML::Node db = YAML::LoadFile(cfg.design.database_path.string());
      YAML::Node design = db["designs"] ? db["designs"][cfg.design.name] : YAML::Node();
      LUPNT_CHECK(design, "GNSS design not found in database: " + cfg.design.name, "LunarGnssODTS");
      ApplyDesignNode(design, cfg);
    }

    void ReadDesignSection(const YAML::Node& root, LunarGnssODTSConfig& cfg,
                           const std::filesystem::path& config_dir) {
      const YAML::Node design = root["design"];
      cfg.design.database_path = ResolvePath(
          config_dir, ReadYaml<std::string>(design, "database", cfg.design.database_path.string()));
      cfg.design.name = ReadYaml(design, "name", cfg.design.name);
      LoadDesignDatabase(cfg);
      ApplyDesignNode(design["custom"], cfg);
    }
  }  // namespace

  void ParsePlasmaDelayConfig(const Config& plasma, PlasmaDelayConfig& cfg) {
    // No `plasma:` block present (e.g. it now lives under `world:`): leave defaults untouched.
    if (!plasma || !plasma.IsMap()) return;
    cfg.simulate_truth = ReadYaml(plasma, "simulate_truth", cfg.simulate_truth);
    cfg.model_in_filter = ReadYaml(plasma, "model_in_filter", cfg.model_in_filter);
    const YAML::Node raytrace = plasma["raytrace"];
    cfg.raytrace_step_size_km = ReadYaml(raytrace, "step_size_km", cfg.raytrace_step_size_km);
    cfg.raytrace_correction = ReadYaml(raytrace, "correction", cfg.raytrace_correction);
    cfg.raytrace_fine_correction
        = ReadYaml(raytrace, "fine_correction", cfg.raytrace_fine_correction);
    cfg.raytrace_straight_ray = ReadYaml(raytrace, "straight_ray", cfg.raytrace_straight_ray);
    cfg.raytrace_compute_higher_order
        = ReadYaml(raytrace, "compute_higher_order", cfg.raytrace_compute_higher_order);
    cfg.raytrace_use_adaptive_step
        = ReadYaml(raytrace, "use_adaptive_step", cfg.raytrace_use_adaptive_step);
    cfg.raytrace_use_fortran_gcpm
        = ReadYaml(raytrace, "use_fortran_gcpm", cfg.raytrace_use_fortran_gcpm);
    cfg.raytrace_cutoff_radius_re
        = ReadYaml(raytrace, "cutoff_radius_re", cfg.raytrace_cutoff_radius_re);
    cfg.raytrace_gradient_step_km
        = ReadYaml(raytrace, "gradient_step_km", cfg.raytrace_gradient_step_km);
    cfg.raytrace_correction_tolerance_m
        = ReadYaml(raytrace, "correction_tolerance_m", cfg.raytrace_correction_tolerance_m);
    cfg.raytrace_kp = ReadYaml(raytrace, "kp", cfg.raytrace_kp);
    cfg.raytrace_rz12 = ReadYaml(raytrace, "rz12", cfg.raytrace_rz12);
    cfg.raytrace_integrator = ReadYaml(raytrace, "integrator", cfg.raytrace_integrator);
    cfg.raytrace_correction_method
        = ReadYaml(raytrace, "correction_method", cfg.raytrace_correction_method);
    cfg.filter_pseudorange_noise_inflation_m = ReadYaml(
        plasma, "filter_pseudorange_noise_inflation_m", cfg.filter_pseudorange_noise_inflation_m);
    cfg.filter_doppler_noise_inflation_hz = ReadYaml(plasma, "filter_doppler_noise_inflation_hz",
                                                     cfg.filter_doppler_noise_inflation_hz);
  }

  LunarGnssODTSConfig ParseLunarGnssODTSConfig(const Config& root,
                                               const std::filesystem::path& config_dir) {
    LunarGnssODTSConfig cfg;

    const YAML::Node sim = root["simulation"];
    cfg.seed = ReadYaml(sim, "seed", cfg.seed);
    cfg.monte_carlo_runs = ReadYaml(sim, "monte_carlo_runs", cfg.monte_carlo_runs);
    cfg.duration_s = ReadYaml(sim, "duration_s", cfg.duration_s);
    cfg.dt_s = ReadYaml(sim, "dt_s", cfg.dt_s);
    cfg.ephemeris_dt_s = ReadYaml(sim, "ephemeris_dt_s", cfg.ephemeris_dt_s);
    cfg.output_dir = ReadYaml<std::string>(sim, "output_dir", cfg.output_dir.string());
    if (cfg.output_dir.is_relative()) cfg.output_dir = std::filesystem::absolute(cfg.output_dir);
    cfg.start_epoch_utc = ReadYaml(sim, "start_epoch_utc", cfg.start_epoch_utc);

    cfg.links_file = cfg.output_dir / "precomputed_links.csv";
    cfg.delays_file = cfg.output_dir / "precomputed_delays.csv";
    const YAML::Node pipeline = root["pipeline"];
    cfg.links_file = ReadYaml<std::string>(pipeline, "links_file", cfg.links_file.string());
    cfg.delays_file = ReadYaml<std::string>(pipeline, "delays_file", cfg.delays_file.string());
    if (cfg.links_file.is_relative()) cfg.links_file = std::filesystem::absolute(cfg.links_file);
    if (cfg.delays_file.is_relative()) cfg.delays_file = std::filesystem::absolute(cfg.delays_file);

    ReadDesignSection(root, cfg, config_dir);

    const YAML::Node app = root["receiver_app"];
    cfg.receiver_app.rate_hz = ReadYaml(app, "rate_hz", cfg.receiver_app.rate_hz);

    const YAML::Node truth = root["truth"];
    cfg.receiver_a_m = ReadYaml(truth, "receiver_a_m", cfg.receiver_a_m);
    cfg.receiver_ecc = ReadYaml(truth, "receiver_ecc", cfg.receiver_ecc);
    cfg.receiver_inc_rad = ReadYaml(truth, "receiver_inc_rad", cfg.receiver_inc_rad);
    cfg.receiver_raan_rad = ReadYaml(truth, "receiver_raan_rad", cfg.receiver_raan_rad);
    cfg.receiver_argp_rad = ReadYaml(truth, "receiver_argp_rad", cfg.receiver_argp_rad);
    cfg.receiver_mean_anomaly_rad
        = ReadYaml(truth, "receiver_mean_anomaly_rad", cfg.receiver_mean_anomaly_rad);
    cfg.clock_bias_s = ReadYaml(truth, "clock_bias_s", cfg.clock_bias_s);
    cfg.clock_drift_sps = ReadYaml(truth, "clock_drift_sps", cfg.clock_drift_sps);
    cfg.clock_model = ReadYaml(truth, "clock_model", cfg.clock_model);
    cfg.use_three_state_clock_truth
        = ReadYaml(truth, "use_three_state_clock_truth", cfg.use_three_state_clock_truth);
    cfg.clock_drift_rate_sps2 = ReadYaml(truth, "clock_drift_rate_sps2", cfg.clock_drift_rate_sps2);

    const YAML::Node constellation = root["constellation"];
    cfg.constellation.sp3_directory
        = ResolvePath(config_dir, ReadYaml<std::string>(constellation, "sp3_directory",
                                                        cfg.constellation.sp3_directory.string()));
    cfg.constellation.auto_select_sp3
        = ReadYaml(constellation, "auto_select_sp3", cfg.constellation.auto_select_sp3);
    cfg.constellation.sp3_files
        = ReadPathVector(constellation, "sp3_files", cfg.constellation.sp3_files, config_dir);
    cfg.constellation.antex_file = ResolvePath(
        config_dir,
        ReadYaml<std::string>(constellation, "antex_file", cfg.constellation.antex_file.string()));
    cfg.constellation.use_all_gps
        = ReadYaml(constellation, "use_all_gps", cfg.constellation.use_all_gps);
    cfg.constellation.include_galileo
        = ReadYaml(constellation, "include_galileo", cfg.constellation.include_galileo);
    cfg.constellation.gps_prns = ReadYaml(constellation, "gps_prns", cfg.constellation.gps_prns);
    cfg.constellation.galileo_prns
        = ReadYaml(constellation, "galileo_prns", cfg.constellation.galileo_prns);
    cfg.constellation.brdc_directory
        = ResolvePath(config_dir, ReadYaml<std::string>(constellation, "brdc_directory",
                                                        cfg.constellation.brdc_directory.string()));
    cfg.constellation.brdc_files
        = ReadPathVector(constellation, "brdc_files", cfg.constellation.brdc_files, config_dir);
    cfg.constellation.use_broadcast_ephemeris = ReadYaml(constellation, "use_broadcast_ephemeris",
                                                         cfg.constellation.use_broadcast_ephemeris);
    cfg.constellation.debias_broadcast_clock = ReadYaml(constellation, "debias_broadcast_clock",
                                                        cfg.constellation.debias_broadcast_clock);
    cfg.constellation.debias_qzss_radial
        = ReadYaml(constellation, "debias_qzss_radial", cfg.constellation.debias_qzss_radial);

    // Plasma is a shared truth-environment property declared under `world: plasma:`; the app
    // reads it from the World in Setup(). This app-block parse remains for back-compat (a
    // `plasma:` block directly under `application:`), and is a no-op when absent.
    ParsePlasmaDelayConfig(root["plasma"], cfg.plasma);

    const YAML::Node dynamics = root["dynamics"];
    cfg.moon_gravity_degree_truth
        = ReadYaml(dynamics, "moon_gravity_degree_truth", cfg.moon_gravity_degree_truth);
    cfg.moon_gravity_order_truth
        = ReadYaml(dynamics, "moon_gravity_order_truth", cfg.moon_gravity_order_truth);
    cfg.moon_gravity_degree_filter
        = ReadYaml(dynamics, "moon_gravity_degree_filter", cfg.moon_gravity_degree_filter);
    cfg.moon_gravity_order_filter
        = ReadYaml(dynamics, "moon_gravity_order_filter", cfg.moon_gravity_order_filter);
    cfg.moon_gravity_degree_constellation = ReadYaml(dynamics, "moon_gravity_degree_constellation",
                                                     cfg.moon_gravity_degree_constellation);
    cfg.moon_gravity_order_constellation = ReadYaml(dynamics, "moon_gravity_order_constellation",
                                                    cfg.moon_gravity_order_constellation);
    cfg.include_earth = ReadYaml(dynamics, "include_earth", cfg.include_earth);
    cfg.include_sun = ReadYaml(dynamics, "include_sun", cfg.include_sun);
    cfg.use_relativity = ReadYaml(dynamics, "use_relativity", cfg.use_relativity);
    cfg.use_srp_truth = ReadYaml(dynamics, "use_srp_truth", cfg.use_srp_truth);
    cfg.use_srp_filter = ReadYaml(dynamics, "use_srp_filter", cfg.use_srp_filter);
    cfg.srp_coeff_truth_m2_kg
        = ReadYaml(dynamics, "srp_coeff_truth_m2_kg", cfg.srp_coeff_truth_m2_kg);
    cfg.srp_coeff_filter_m2_kg
        = ReadYaml(dynamics, "srp_coeff_filter_m2_kg", cfg.srp_coeff_filter_m2_kg);

    const YAML::Node meas = root["measurements"];
    cfg.use_pseudorange = ReadYaml(meas, "use_pseudorange", cfg.use_pseudorange);
    cfg.use_doppler = ReadYaml(meas, "use_doppler", cfg.use_doppler);
    cfg.use_tdcp = ReadYaml(meas, "use_tdcp", cfg.use_tdcp);
    cfg.pseudorange_sigma_m = ReadYaml(meas, "pseudorange_sigma_m", cfg.pseudorange_sigma_m);
    cfg.doppler_sigma_hz = ReadYaml(meas, "doppler_sigma_hz", cfg.doppler_sigma_hz);
    cfg.carrier_phase_sigma_m = ReadYaml(meas, "carrier_phase_sigma_m", cfg.carrier_phase_sigma_m);
    cfg.tdcp_sigma_m = ReadYaml(meas, "tdcp_sigma_m", cfg.tdcp_sigma_m);
    cfg.use_ionosphere_free = ReadYaml(meas, "use_ionosphere_free", cfg.use_ionosphere_free);
    cfg.filter_pseudorange_noise_inflation_m = ReadYaml(
        meas, "filter_pseudorange_noise_inflation_m", cfg.filter_pseudorange_noise_inflation_m);
    cfg.filter_tdcp_noise_inflation_m
        = ReadYaml(meas, "filter_tdcp_noise_inflation_m", cfg.filter_tdcp_noise_inflation_m);
    cfg.pseudorange_min_tangent_altitude_m = ReadYaml(meas, "pseudorange_min_tangent_altitude_m",
                                                      cfg.pseudorange_min_tangent_altitude_m);
    cfg.tdcp_min_tangent_altitude_m
        = ReadYaml(meas, "tdcp_min_tangent_altitude_m", cfg.tdcp_min_tangent_altitude_m);

    const YAML::Node filter = root["filter"];
    cfg.estimate_srp_coefficient
        = ReadYaml(filter, "estimate_srp_coefficient", cfg.estimate_srp_coefficient);
    cfg.initial_srp_coeff_m2_kg
        = ReadYaml(filter, "initial_srp_coeff_m2_kg", cfg.initial_srp_coeff_m2_kg);
    cfg.initial_position_sigma_m
        = ReadYaml(filter, "initial_position_sigma_m", cfg.initial_position_sigma_m);
    cfg.initial_velocity_sigma_mps
        = ReadYaml(filter, "initial_velocity_sigma_mps", cfg.initial_velocity_sigma_mps);
    cfg.initial_clock_bias_sigma_s
        = ReadYaml(filter, "initial_clock_bias_sigma_s", cfg.initial_clock_bias_sigma_s);
    cfg.initial_clock_drift_sigma_sps
        = ReadYaml(filter, "initial_clock_drift_sigma_sps", cfg.initial_clock_drift_sigma_sps);
    cfg.initial_srp_coeff_sigma_m2_kg
        = ReadYaml(filter, "initial_srp_coeff_sigma_m2_kg", cfg.initial_srp_coeff_sigma_m2_kg);
    cfg.process_accel_sigma_mps2
        = ReadYaml(filter, "process_accel_sigma_mps2", cfg.process_accel_sigma_mps2);
    cfg.process_srp_coeff_sigma_m2_kg_sqrt_s = ReadYaml(
        filter, "process_srp_coeff_sigma_m2_kg_sqrt_s", cfg.process_srp_coeff_sigma_m2_kg_sqrt_s);
    cfg.integration_step_s = ReadYaml(filter, "integration_step_s", cfg.integration_step_s);
    cfg.precompute_progress_interval_s
        = ReadYaml(sim, "precompute_progress_interval_s", cfg.precompute_progress_interval_s);
    cfg.run_progress_interval_epochs
        = ReadYaml(sim, "run_progress_interval_epochs", cfg.run_progress_interval_epochs);
    cfg.debug_print_matrix_epochs
        = ReadYaml(sim, "debug_print_matrix_epochs", cfg.debug_print_matrix_epochs);
    cfg.debug_print_matrix_max_rows
        = ReadYaml(sim, "debug_print_matrix_max_rows", cfg.debug_print_matrix_max_rows);
    cfg.precompute_num_threads
        = ReadYaml(sim, "precompute_num_threads", cfg.precompute_num_threads);

    LUPNT_CHECK(cfg.dt_s > 0.0, "GNSS filtering dt_s must be positive", "LunarGnssODTS");
    LUPNT_CHECK(cfg.receiver_app.rate_hz > 0.0, "GNSS receiver app rate_hz must be positive",
                "LunarGnssODTS");
    LUPNT_CHECK(cfg.ephemeris_dt_s > 0.0, "GNSS filtering ephemeris_dt_s must be positive",
                "LunarGnssODTS");
    LUPNT_CHECK(cfg.duration_s >= 0.0, "GNSS filtering duration_s must be non-negative",
                "LunarGnssODTS");
    LUPNT_CHECK(cfg.receiver_a_m > R_MOON, "ELFO semi-major axis must exceed Moon radius",
                "LunarGnssODTS");
    LUPNT_CHECK(cfg.receiver_ecc >= 0.0 && cfg.receiver_ecc < 1.0,
                "ELFO eccentricity must be in [0, 1)", "LunarGnssODTS");
    if (!cfg.constellation.auto_select_sp3) {
      LUPNT_CHECK(!cfg.constellation.sp3_files.empty(), "GNSS filtering requires SP3 files",
                  "LunarGnssODTS");
    }
    if (cfg.constellation.use_all_gps) {
      cfg.constellation.gps_prns.clear();
    } else {
      LUPNT_CHECK(!cfg.constellation.gps_prns.empty(),
                  "Set constellation.use_all_gps=true or provide constellation.gps_prns",
                  "LunarGnssODTS");
    }

    ResolveAutoSelectSp3Files(cfg);
    return cfg;
  }

  void ResolveLunarGnssODTSConfigForRun(LunarGnssODTSConfig& cfg) {
    ResolveAutoSelectSp3Files(cfg);
  }

  LunarGnssODTSConfig LoadLunarGnssODTSConfig(const std::filesystem::path& path) {
    YAML::Node root = YAML::LoadFile(path.string());
    std::filesystem::path config_dir
        = path.has_parent_path() ? path.parent_path() : std::filesystem::current_path();
    return ParseLunarGnssODTSConfig(root, config_dir);
  }

  std::vector<LunarGnssODTSSummary> RunLunarGnssODTSMonteCarlo(
      const LunarGnssODTSConfig& cfg, const std::vector<State>* receiver_truth) {
    ClearPreviousOutputs(cfg);
    RequireDelayTable(cfg);

    Real t0_tdb = ConvertTime(GregorianToTime(cfg.start_epoch_utc), Time::UTC, Time::TDB);
    std::vector<ScheduledAppCall> receiver_app_calls = BuildReceiverAppSchedule(cfg);
    VecXd elapsed_s = ElapsedTimesFromSchedule(receiver_app_calls);
    VecXd receiver_clock_s = ReceiverClockTimesFromSchedule(receiver_app_calls);
    VecXd ephem_elapsed_s = MakeEphemerisTimeVector(cfg);
    VecXd times_tdb = ShiftTimes(t0_tdb, elapsed_s);
    VecXd ephem_times_tdb = ShiftTimes(t0_tdb, ephem_elapsed_s);
    VecXd ephem_times_tai = ConvertTimeVector(ephem_times_tdb, Time::TDB, Time::TAI);

    // Receiver truth: prefer the physical Spacecraft's self-propagated grid (agent-driven);
    // else build the nominal (noiseless) trajectory internally. The link geometry only needs
    // the orbit, so a noisy-clock agent grid gives identical links.
    std::vector<State> nominal_truth_states
        = receiver_truth ? *receiver_truth
                         : BuildTruthStates(cfg, t0_tdb, elapsed_s, times_tdb, false, cfg.seed,
                                            "Truth trajectory");
    std::vector<GNSSMeasurementsEpoch> precomputed;
    const std::string link_fingerprint = LinkCacheFingerprint(cfg);
    if (LinkCacheMatchesConfig(cfg, link_fingerprint)) {
      precomputed = LoadLinksFile(cfg, times_tdb);
    } else {
      LUPNT_CHECK(
          !FileHasContent(cfg.links_file),
          "GNSS link cache is stale for this config: " + cfg.links_file.string()
              + ". Regenerate it with `pixi run python python/examples/ex6_precompute.py` before "
                "running the EKF.",
          "LunarGnssODTS");
      Logger::Info("No GNSS link cache found; computing links in-memory for this run.",
                   "LunarGnssODTS");
      std::vector<RuntimeConstellation> constellations = BuildConstellations(cfg, ephem_times_tai);
      std::vector<Real> receive_times = ToRealVector(times_tdb);
      LunarGnssODTSConfig link_cfg = cfg;
      link_cfg.plasma.simulate_truth = false;
      precomputed = PrecomputeConstellationChannels(link_cfg, constellations, receive_times,
                                                    nominal_truth_states);
    }
    if (cfg.plasma.simulate_truth) ApplyDelayFile(precomputed, cfg.delays_file);

    std::vector<LunarGnssODTSSummary> summaries;
    summaries.reserve(cfg.monte_carlo_runs);
    for (int mc = 0; mc < cfg.monte_carlo_runs; ++mc) {
      summaries.push_back(RunOneMonteCarlo(cfg, mc, t0_tdb, elapsed_s, times_tdb, receiver_clock_s,
                                           precomputed, receiver_truth));
    }
    WriteSummary(cfg, summaries);
    return summaries;
  }

  int LunarGnssODTSPrecomputeEpochCount(const LunarGnssODTSConfig& cfg) {
    return static_cast<int>(BuildReceiverAppSchedule(cfg).size());
  }

  VecXd LunarGnssODTSReceiverElapsedTimes(const LunarGnssODTSConfig& cfg) {
    return ElapsedTimesFromSchedule(BuildReceiverAppSchedule(cfg));
  }

  bool LunarGnssODTSLinkCacheValid(const LunarGnssODTSConfig& cfg) {
    LunarGnssODTSConfig resolved = ResolveRuntimeConfig(cfg);
    return LinkCacheMatchesConfig(resolved, LinkCacheFingerprint(resolved));
  }

  void FinalizeLunarGnssODTSLinkCache(const LunarGnssODTSConfig& cfg) {
    LunarGnssODTSConfig resolved = ResolveRuntimeConfig(cfg);
    WriteLinkCacheMetadata(resolved, LinkCacheFingerprint(resolved));
  }

  void PrecomputeLunarGnssODTSLinks(const LunarGnssODTSConfig& input_cfg) {
    LunarGnssODTSConfig cfg = ResolveRuntimeConfig(input_cfg);
    const std::string fingerprint = LinkCacheFingerprint(cfg);
    if (LinkCacheMatchesConfig(cfg, fingerprint)) {
      Logger::Info("Using cached GNSS link precompute: " + cfg.links_file.string(),
                   "LunarGnssODTS");
      return;
    }
    if (FileHasContent(cfg.links_file)) {
      Logger::Info(
          "GNSS link cache is stale for this config; regenerating " + cfg.links_file.string(),
          "LunarGnssODTS");
    }
    RemoveStaleDelayCache(cfg);

    Real t0_tdb = ConvertTime(GregorianToTime(cfg.start_epoch_utc), Time::UTC, Time::TDB);
    std::vector<ScheduledAppCall> receiver_app_calls = BuildReceiverAppSchedule(cfg);
    VecXd elapsed_s = ElapsedTimesFromSchedule(receiver_app_calls);
    VecXd ephem_elapsed_s = MakeEphemerisTimeVector(cfg);
    VecXd times_tdb = ShiftTimes(t0_tdb, elapsed_s);
    VecXd ephem_times_tdb = ShiftTimes(t0_tdb, ephem_elapsed_s);
    VecXd ephem_times_tai = ConvertTimeVector(ephem_times_tdb, Time::TDB, Time::TAI);

    std::vector<State> nominal_truth_states
        = BuildTruthStates(cfg, t0_tdb, elapsed_s, times_tdb, false, cfg.seed, "Truth trajectory");
    std::vector<RuntimeConstellation> constellations = BuildConstellations(cfg, ephem_times_tai);
    std::vector<Real> receive_times = ToRealVector(times_tdb);
    LunarGnssODTSConfig link_cfg = cfg;
    link_cfg.plasma.simulate_truth = false;
    std::vector<GNSSMeasurementsEpoch> precomputed = PrecomputeConstellationChannels(
        link_cfg, constellations, receive_times, nominal_truth_states);
    WriteLinksCsv(cfg, precomputed, nominal_truth_states);
    WriteLinkCacheMetadata(cfg, fingerprint);
  }

  void PrecomputeLunarGnssODTSLinksRange(const LunarGnssODTSConfig& input_cfg, int epoch_begin,
                                         int epoch_end) {
    LunarGnssODTSConfig cfg = ResolveRuntimeConfig(input_cfg);
    Real t0_tdb = ConvertTime(GregorianToTime(cfg.start_epoch_utc), Time::UTC, Time::TDB);
    std::vector<ScheduledAppCall> receiver_app_calls = BuildReceiverAppSchedule(cfg);
    const int n_epochs = static_cast<int>(receiver_app_calls.size());
    LUPNT_CHECK(epoch_begin >= 0 && epoch_begin <= epoch_end && epoch_end <= n_epochs,
                "Invalid GNSS link precompute epoch range [" + std::to_string(epoch_begin) + ", "
                    + std::to_string(epoch_end) + ") for " + std::to_string(n_epochs) + " epochs",
                "LunarGnssODTS");

    VecXd elapsed_s = ElapsedTimesFromSchedule(receiver_app_calls);
    VecXd ephem_elapsed_s = MakeEphemerisTimeVector(cfg);
    VecXd times_tdb = ShiftTimes(t0_tdb, elapsed_s);
    VecXd ephem_times_tdb = ShiftTimes(t0_tdb, ephem_elapsed_s);
    VecXd ephem_times_tai = ConvertTimeVector(ephem_times_tdb, Time::TDB, Time::TAI);

    std::vector<State> nominal_truth_states
        = BuildTruthStates(cfg, t0_tdb, elapsed_s, times_tdb, false, cfg.seed, "Truth trajectory");
    std::vector<RuntimeConstellation> constellations = BuildConstellations(cfg, ephem_times_tai);

    std::vector<Real> receive_times;
    std::vector<State> truth_states;
    receive_times.reserve(epoch_end - epoch_begin);
    truth_states.reserve(epoch_end - epoch_begin);
    for (int i = epoch_begin; i < epoch_end; ++i) {
      receive_times.push_back(times_tdb(i));
      truth_states.push_back(nominal_truth_states[i]);
    }

    LunarGnssODTSConfig link_cfg = cfg;
    link_cfg.plasma.simulate_truth = false;
    std::vector<GNSSMeasurementsEpoch> precomputed
        = PrecomputeConstellationChannels(link_cfg, constellations, receive_times, truth_states);
    WriteLinksCsv(cfg, precomputed, truth_states, epoch_begin);
  }

}  // namespace lupnt

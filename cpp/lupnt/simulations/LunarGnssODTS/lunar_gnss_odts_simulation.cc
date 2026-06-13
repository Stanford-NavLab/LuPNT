#include "lupnt/simulations/LunarGnssODTS/lunar_gnss_odts_simulation.h"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
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
      AppendFingerprintField(oss, "version", 14);
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
      AppendFingerprintField(oss, "apply_cn0_threshold", cfg.design.apply_cn0_threshold);
      AppendFingerprintField(oss, "cn0_threshold_dbhz", cfg.design.cn0_threshold_dbhz);
      AppendFingerprintField(oss, "use_pseudorange", cfg.use_pseudorange);
      AppendFingerprintField(oss, "use_doppler", cfg.use_doppler);
      AppendFingerprintField(oss, "use_tdcp", cfg.use_tdcp);
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

    Ptr<ClockDynamics> CreateClockDynamics(bool add_noise, int seed = 0) {
      auto clock = MakePtr<ClockDynamics>();
      clock->SetModel(ClockModel::OCXO);
      clock->SetClockBiasUnit(ClockBiasUnit::SECONDS);
      clock->SetAddNoise(add_noise);
      clock->SetSeed(seed);
      return clock;
    }

    Ptr<JointOrbitClockDynamics> CreateJointDynamics(Ptr<NumericalDynamics> orbit_dynamics,
                                                     bool add_clock_noise, int seed = 0) {
      auto joint = MakePtr<JointOrbitClockDynamics>();
      joint->SetOrbitDynamics(orbit_dynamics);
      joint->SetClockDynamics(CreateClockDynamics(add_clock_noise, seed));
      joint->SetUseClockRelativity(true);
      joint->SetRelativityCenterBody(BodyId::MOON);
      joint->SetAddClockNoise(add_clock_noise);
      joint->SetFrame(Frame::MOON_CI);
      joint->SetTimeStep(orbit_dynamics->GetTimeStep());
      joint->SetIntegrator(IntegratorType::RKF45);
      joint->SetIntegratorParams(IntegratorParams(20, 1.0e-12, 1.0e-12));
      return joint;
    }

    JointOrbitClockState MakeJointState(const Vec6& rv, double clock_bias_s,
                                        double clock_drift_sps) {
      Cart6 orbit(rv, Frame::MOON_CI);
      ClockState2 clock;
      clock.b() = clock_bias_s;
      clock.d() = clock_drift_sps;
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
      x_joint.SetUnits({"m", "m", "m", "m/s", "m/s", "m/s", "s", "s/s"});
      return JointOrbitClockState(x_joint);
    }

    State MakeAugmentedSrpState(const State& joint_state, double srp_coeff_m2_kg) {
      State x(9);
      x.head(8) = joint_state.head(8);
      x(8) = srp_coeff_m2_kg;
      x.SetFrame(joint_state.GetFrame());
      x.SetName("JointOrbitClockSrp");
      x.SetNames({"r_x", "r_y", "r_z", "v_x", "v_y", "v_z", "b", "d", "bcoeff_srp"});
      x.SetUnits({"m", "m", "m", "m/s", "m/s", "m/s", "s", "s/s", "m^2/kg"});
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

    ProcessNoiseFunction CreateProcessNoiseFunction(const LunarGnssODTSConfig& cfg) {
      return [cfg](const State& x, Real t0, Real tf) {
        const double dt = std::abs((tf - t0).val());
        MatXd Q = MatXd::Zero(x.size(), x.size());

        Mat3d Q_acc = std::pow(cfg.process_accel_sigma_mps2, 2) * Mat3d::Identity();
        Q.block(0, 0, 6, 6) = ProcessNoisePosVel(Q_acc, dt);

        Q.block(6, 6, 2, 2)
            = ClockDynamics::TwoStateNoise(ClockModel::OCXO, dt, ClockBiasUnit::SECONDS)
                  .cast<double>();
        if (EstimateSrp(cfg)) {
          Q(8, 8) = std::pow(cfg.process_srp_coeff_sigma_m2_kg_sqrt_s, 2) * dt;
        }
        return Q;
      };
    }

    MatXd InitialCovariance(const LunarGnssODTSConfig& cfg) {
      const int n = EstimateSrp(cfg) ? 9 : 8;
      MatXd P = MatXd::Zero(n, n);
      P.block(0, 0, 3, 3) = std::pow(cfg.initial_position_sigma_m, 2) * Mat3d::Identity();
      P.block(3, 3, 3, 3) = std::pow(cfg.initial_velocity_sigma_mps, 2) * Mat3d::Identity();
      P(6, 6) = std::pow(cfg.initial_clock_bias_sigma_s, 2);
      P(7, 7) = std::pow(cfg.initial_clock_drift_sigma_sps, 2);
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
      err(6) = SampleNormal(0.0, cfg.initial_clock_bias_sigma_s).val();
      err(7) = SampleNormal(0.0, cfg.initial_clock_drift_sigma_sps).val();
      if (EstimateSrp(cfg)) err(8) = SampleNormal(0.0, cfg.initial_srp_coeff_sigma_m2_kg).val();
      return err;
    }

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
      LUPNT_CHECK(!cfg.constellation.sp3_files.empty(), "GNSS filtering requires SP3 files",
                  "LunarGnssODTS");
      for (const auto& sp3 : cfg.constellation.sp3_files) {
        LUPNT_CHECK(std::filesystem::exists(sp3), "GNSS SP3 file not found: " + sp3.string(),
                    "LunarGnssODTS");
      }
      LUPNT_CHECK(std::filesystem::exists(cfg.constellation.antex_file),
                  "GNSS ANTEX file not found: " + cfg.constellation.antex_file.string(),
                  "LunarGnssODTS");

      const std::vector<int> gps_prns = cfg.constellation.gps_prns;
      out.push_back(BuildConstellationFromFiles(GnssConst::GPS, GnssFreq::L1, gps_prns, cfg,
                                                ephem_times_tai));

      if (cfg.constellation.include_galileo) {
        const std::vector<int> galileo_prns = cfg.constellation.galileo_prns;
        out.push_back(BuildConstellationFromFiles(GnssConst::GALILEO, GnssFreq::E1, galileo_prns,
                                                  cfg, ephem_times_tai));
      }
      return out;
    }

    std::vector<State> BuildTruthStates(const LunarGnssODTSConfig& cfg, Real t0_tdb,
                                        const VecXd& elapsed_s, const VecXd& times_tdb,
                                        bool add_clock_noise, int seed) {
      (void)elapsed_s;
      Ptr<NBodyDynamics> orbit_dynamics
          = CreateLunarDynamics(cfg.moon_gravity_degree_truth, cfg.moon_gravity_order_truth, false,
                                cfg, cfg.use_srp_truth, cfg.srp_coeff_truth_m2_kg);
      Ptr<JointOrbitClockDynamics> dynamics
          = CreateJointDynamics(orbit_dynamics, add_clock_noise, seed);
      Vec6 rv0 = ReceiverElfoInitialState(t0_tdb, cfg);
      State x = MakeJointState(rv0, cfg.clock_bias_s, cfg.clock_drift_sps);
      std::vector<State> states;
      states.reserve(times_tdb.size());
      for (int i = 0; i < times_tdb.size(); ++i) {
        if (i == 0) {
          states.push_back(x);
        } else {
          x = dynamics->Propagate(JointOrbitClockState(x), times_tdb(i - 1), times_tdb(i), nullptr);
          states.push_back(x);
        }
      }
      return states;
    }

    GnssMeasurementOptions MeasurementOptions(const LunarGnssODTSConfig& cfg, bool truth_model) {
      GnssMeasurementOptions options;
      options.observables = CurrentEpochObservables(cfg);
      options.clock_bias_unit = ClockBiasUnit::SECONDS;
      options.receive_time_scale = Time::TDB;
      options.ephemeris_time_scale = Time::TAI;
      options.frame = Frame::MOON_CI;
      options.solve_light_time = true;
      options.recompute_transmit_time_from_ephemeris = !truth_model;
      options.apply_transmitter_relativity = false;
      options.apply_shapiro_delay = truth_model;
      options.apply_visibility = true;
      options.apply_cn0_threshold = cfg.design.apply_cn0_threshold;
      options.cn0_threshold_dbhz = cfg.design.cn0_threshold_dbhz;
      options.apply_ionosphere_plasma_delay = truth_model && cfg.plasma.simulate_truth;
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
      for (auto& channel : channels) channel.shapiro_delay_m = 0.0;
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
          channel.sigma_pseudorange_m = std::hypot(channel.sigma_pseudorange_m.val(),
                                                   cfg.plasma.filter_pseudorange_noise_inflation_m);
          channel.sigma_doppler_hz = std::hypot(channel.sigma_doppler_hz.val(),
                                                cfg.plasma.filter_doppler_noise_inflation_hz);
        }
      }
    }

    VecXd ComputeMeasurementVector(const State& x, const std::vector<GnssChannel>& channels,
                                   const GnssMeasurementOptions& options, MatXd* H = nullptr) {
      const int n_obs = static_cast<int>(options.observables.size());
      VecXd y(n_obs * static_cast<int>(channels.size()));
      if (H != nullptr) H->setZero(y.size(), x.size());

      for (int i = 0; i < static_cast<int>(channels.size()); ++i) {
        GnssMeasurement measurement(channels[i]);
        MatXd H_i;
        VecXd y_i = measurement.Compute(x, H != nullptr ? &H_i : nullptr, options);
        y.segment(i * n_obs, n_obs) = y_i;
        if (H != nullptr) H->block(i * n_obs, 0, n_obs, x.size()) = H_i;
      }
      return y;
    }

    MatXd MeasurementCovariance(const std::vector<GnssChannel>& channels,
                                const GnssMeasurementOptions& options,
                                const LunarGnssODTSConfig& cfg) {
      (void)cfg;
      const int n_obs = static_cast<int>(options.observables.size());
      MatXd R = MatXd::Zero(n_obs * static_cast<int>(channels.size()),
                            n_obs * static_cast<int>(channels.size()));
      for (int i = 0; i < static_cast<int>(channels.size()); ++i) {
        for (int j = 0; j < n_obs; ++j) {
          const int row = i * n_obs + j;
          if (options.observables[j] == GnssObservable::PSEUDORANGE)
            R(row, row) = std::pow(channels[i].sigma_pseudorange_m.val(), 2);
          else if (options.observables[j] == GnssObservable::DOPPLER)
            R(row, row) = std::pow(channels[i].sigma_doppler_hz.val(), 2);
          else if (options.observables[j] == GnssObservable::CARRIER_PHASE)
            R(row, row) = std::pow(channels[i].sigma_carrier_phase_cycles.val(), 2);
        }
      }
      return R;
    }

    struct TdcpPair {
      GnssChannel current;
      GnssChannel previous;
    };

    std::vector<TdcpPair> MakeTdcpPairs(const std::vector<GnssChannel>& current,
                                        const std::vector<GnssChannel>& previous) {
      std::vector<TdcpPair> pairs;
      for (const auto& curr : current) {
        auto it = std::find_if(previous.begin(), previous.end(), [&](const GnssChannel& prev) {
          return curr.gnss_const == prev.gnss_const && curr.prn == prev.prn
                 && curr.frequency == prev.frequency;
        });
        if (it != previous.end()) pairs.push_back({curr, *it});
      }
      return pairs;
    }

    Real CarrierRangeMeters(const State& x, const GnssChannel& channel,
                            const GnssMeasurementOptions& options, MatXd* H = nullptr) {
      GnssMeasurement measurement(channel);
      MatXd H_cycles;
      VecXd y_cycles = measurement.Compute(x, H != nullptr ? &H_cycles : nullptr, options);
      const Real lambda = channel.Wavelength();
      if (H != nullptr) *H = lambda.val() * H_cycles;
      return lambda * y_cycles(0);
    }

    VecXd ComputeTdcpVector(const State& current_state, const State& previous_state,
                            const std::vector<TdcpPair>& pairs,
                            const GnssMeasurementOptions& options, MatXd* H = nullptr) {
      VecXd y(pairs.size());
      if (H != nullptr) H->setZero(pairs.size(), current_state.size() + previous_state.size());
      for (int i = 0; i < static_cast<int>(pairs.size()); ++i) {
        MatXd H_curr;
        MatXd H_prev;
        Real curr_m = CarrierRangeMeters(current_state, pairs[i].current, options,
                                         H != nullptr ? &H_curr : nullptr);
        Real prev_m = CarrierRangeMeters(previous_state, pairs[i].previous, options,
                                         H != nullptr ? &H_prev : nullptr);
        y(i) = curr_m - prev_m;
        if (H != nullptr) {
          H->block(i, 0, 1, current_state.size()) = H_curr;
          H->block(i, current_state.size(), 1, previous_state.size()) = -H_prev;
        }
      }
      return y;
    }

    MatXd TdcpCovariance(const std::vector<TdcpPair>& pairs, const LunarGnssODTSConfig& cfg) {
      MatXd R = MatXd::Zero(pairs.size(), pairs.size());
      for (int i = 0; i < static_cast<int>(pairs.size()); ++i) {
        const double curr_sigma_m = pairs[i].current.sigma_carrier_phase_cycles.val()
                                    * pairs[i].current.Wavelength().val();
        const double prev_sigma_m = pairs[i].previous.sigma_carrier_phase_cycles.val()
                                    * pairs[i].previous.Wavelength().val();
        R(i, i)
            = std::pow(curr_sigma_m, 2) + std::pow(prev_sigma_m, 2) + std::pow(cfg.tdcp_sigma_m, 2);
      }
      return R;
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

    FilterMeasurementFunction CreateCombinedMeasurementFunction(
        std::vector<GnssChannel> channels, std::vector<TdcpPair> tdcp_pairs,
        GnssMeasurementOptions current_options, GnssMeasurementOptions carrier_options,
        LunarGnssODTSConfig cfg) {
      return [channels = std::move(channels), tdcp_pairs = std::move(tdcp_pairs),
              current_options = std::move(current_options),
              carrier_options = std::move(carrier_options),
              cfg = std::move(cfg)](const State& x, MatXd* H, MatXd* R) {
        const int base_n
            = cfg.use_tdcp ? static_cast<int>(x.size() / 2) : static_cast<int>(x.size());
        State x_curr = cfg.use_tdcp ? State(x.head(base_n)) : x;
        State x_prev = cfg.use_tdcp ? State(x.tail(base_n)) : x;

        MatXd H_current_base;
        VecXd y_current = ComputeMeasurementVector(x_curr, channels, current_options,
                                                   H != nullptr ? &H_current_base : nullptr);
        MatXd H_tdcp;
        VecXd y_tdcp = cfg.use_tdcp ? ComputeTdcpVector(x_curr, x_prev, tdcp_pairs, carrier_options,
                                                        H != nullptr ? &H_tdcp : nullptr)
                                    : VecXd(0);

        VecXd y(y_current.size() + y_tdcp.size());
        y << y_current, y_tdcp;

        if (H != nullptr) {
          H->setZero(y.size(), x.size());
          H->block(0, 0, y_current.size(), base_n) = H_current_base;
          if (cfg.use_tdcp && y_tdcp.size() > 0) {
            H->block(y_current.size(), 0, y_tdcp.size(), 2 * base_n) = H_tdcp;
          }
        }

        if (R != nullptr) {
          MatXd R_current = MeasurementCovariance(channels, current_options, cfg);
          MatXd R_tdcp = cfg.use_tdcp ? TdcpCovariance(tdcp_pairs, cfg) : MatXd::Zero(0, 0);
          R->setZero(y.size(), y.size());
          if (R_current.size() > 0)
            R->topLeftCorner(R_current.rows(), R_current.cols()) = R_current;
          if (R_tdcp.size() > 0) {
            R->bottomRightCorner(R_tdcp.rows(), R_tdcp.cols()) = R_tdcp;
          }
        }
        return y;
      };
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
      MatXd R = TdcpCovariance(pairs, cfg);
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
      if (!GNSSMeasurements::ComputeVisibility(rx_eci, tx_eci, R_EARTH, Vec3::Zero())) return false;

      Vec3 tx_mci = ConvertFrame(t_tdb, tx_eci, Frame::ECI, Frame::MOON_CI);
      return GNSSMeasurements::ComputeVisibility(rx_mci, tx_mci, R_MOON, Vec3::Zero());
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

      for (const auto& runtime : constellations) {
        GNSSMeasurements gnss(runtime.constellation);
        gnss.SetFrequency(runtime.frequency);
        GnssMeasurementOptions candidate_options = MeasurementOptions(cfg, true);
        candidate_options.frame = Frame::ECI;
        candidate_options.apply_visibility = false;
        gnss.SetOptions(candidate_options);
        gnss.SetReceiverParams(cfg.design.receiver_params);
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
        }
        for (size_t i = 0; i < epochs.size(); ++i) {
          combined[i].channels.insert(combined[i].channels.end(), epochs[i].channels.begin(),
                                      epochs[i].channels.end());
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
      std::vector<std::string> columns = SplitCsvLine(header);
      std::unordered_map<std::string, int> idx;
      for (int i = 0; i < static_cast<int>(columns.size()); ++i) idx[columns[i]] = i;

      auto col = [&](const std::vector<std::string>& row, const std::string& name) -> std::string {
        auto it = idx.find(name);
        if (it == idx.end() || it->second >= static_cast<int>(row.size())) return "";
        return row[it->second];
      };

      std::string line;
      while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row = SplitCsvLine(line);
        LinkKey key;
        key.epoch_index = std::stoi(col(row, "epoch_index"));
        key.gnss_const = std::stoi(col(row, "gnss_const_id"));
        key.prn = std::stoi(col(row, "prn"));
        key.frequency = std::stoi(col(row, "frequency_id"));
        DelayRecord record;
        std::string delay = col(row, "ionosphere_plasma_delay_m");
        if (delay.empty()) delay = col(row, "total_delay_m");
        record.ionosphere_plasma_delay_m = delay.empty() ? 0.0 : std::stod(delay);
        delays[key] = record;
      }
      return delays;
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
                       const std::vector<State>& truth_states) {
      std::filesystem::create_directories(cfg.links_file.parent_path());
      std::ofstream out(cfg.links_file);
      out << std::scientific << std::setprecision(12);
      out << "link_id,epoch_index,t_tdb,t_utc,gnss_const_id,gnss_const,prn,frequency_id,frequency,"
             "rx_x_ecef_m,rx_y_ecef_m,rx_z_ecef_m,tx_x_ecef_m,tx_y_ecef_m,tx_z_ecef_m,"
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
          auto [tangent_radius_m, tangent_altitude_m]
              = EarthTangentRadiusAltitude(rx_ecef, tx_ecef);
          double tx_boresight_angle_deg = TransmitterBoresightAngleDeg(t_tdb, rx_mci, channel);
          out << link_id++ << "," << k << "," << t_tdb.val() << "," << t_utc.val() << ","
              << static_cast<int>(channel.gnss_const) << "," << GnssConstName(channel.gnss_const)
              << "," << channel.prn << "," << static_cast<int>(channel.frequency) << ","
              << GnssFreqName(channel.frequency) << "," << rx_ecef(0).val() << ","
              << rx_ecef(1).val() << "," << rx_ecef(2).val() << "," << tx_ecef(0).val() << ","
              << tx_ecef(1).val() << "," << tx_ecef(2).val() << "," << tangent_radius_m << ","
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
        Ptr<JointOrbitClockDynamics> dynamics_filter = CreateJointDynamics(
            orbit_dynamics_filter, false, context_.cfg.seed + 2000 + context_.mc_index);
        filter_->SetDynamicsFunction(CreateDynamicsFunction(dynamics_filter, context_.cfg));
        filter_->SetProcessNoiseFunction(CreateProcessNoiseFunction(context_.cfg));
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

        summary_.monte_carlo_index = context_.mc_index;
        summary_.num_epochs = static_cast<int>(context_.times_tdb.size());
        pos_err2_sum_ = 0.0;
        vel_err2_sum_ = 0.0;
        next_epoch_index_ = 0;
      }

      void Step(Real) override {
        const int k = next_epoch_index_++;
        LUPNT_CHECK(k < context_.times_tdb.size(), "Lunar ODTS app called after final epoch",
                    "LunarGnssODTS");
        std::vector<GnssChannel> truth_channels = context_.precomputed[k].channels;
        ApplyMeasurementSigmas(truth_channels, context_.cfg, false);
        std::vector<GnssChannel> filter_channels = MakeFilterChannels(truth_channels, context_.cfg);
        ApplyMeasurementSigmas(filter_channels, context_.cfg, true);

        if (k > 0) filter_->Predict(context_.times_tdb(k));

        VecXd y_true
            = ComputeMeasurementVector(context_.truth_states[k], truth_channels, truth_options_);
        VecXd y_obs = AddMeasurementNoise(y_true, truth_channels, truth_options_);
        std::vector<TdcpPair> truth_tdcp_pairs;
        std::vector<TdcpPair> filter_tdcp_pairs;
        if (context_.cfg.use_tdcp && k > 0) {
          truth_tdcp_pairs = MakeTdcpPairs(truth_channels, previous_truth_channels_);
          filter_tdcp_pairs = MakeTdcpPairs(filter_channels, previous_filter_channels_);
          if (!truth_tdcp_pairs.empty() && !filter_tdcp_pairs.empty()) {
            VecXd y_tdcp_true
                = ComputeTdcpVector(context_.truth_states[k], context_.truth_states[k - 1],
                                    truth_tdcp_pairs, truth_carrier_options_);
            VecXd y_tdcp_obs = AddTdcpNoise(y_tdcp_true, truth_tdcp_pairs, context_.cfg);
            VecXd y_combined(y_obs.size() + y_tdcp_obs.size());
            y_combined << y_obs, y_tdcp_obs;
            y_obs = y_combined;
          }
        }

        filter_->SetMeasurementFunction(
            CreateCombinedMeasurementFunction(filter_channels, filter_tdcp_pairs, filter_options_,
                                              filter_carrier_options_, context_.cfg));
        filter_->Update(y_obs);

        State x_est = CurrentFilterState(filter_->GetState(), context_.cfg);
        Vec3d dr = (context_.truth_states[k].head(3) - x_est.head(3)).cast<double>();
        Vec3d dv = (context_.truth_states[k].segment(3, 3) - x_est.segment(3, 3)).cast<double>();
        const double pos_err = dr.norm();
        const double vel_err = dv.norm();
        const double clk_b_err_m = C * (context_.truth_states[k](6) - x_est(6)).val();
        const double clk_d_err_mps = C * (context_.truth_states[k](7) - x_est(7)).val();
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
        const double sb = 3.0 * C * std::sqrt(std::max(0.0, P(6, 6)));
        const double sd = 3.0 * C * std::sqrt(std::max(0.0, P(7, 7)));
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

        previous_truth_channels_ = truth_channels;
        previous_filter_channels_ = filter_channels;

        if (k == context_.times_tdb.size() - 1) {
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
      std::vector<GnssChannel> previous_truth_channels_;
      std::vector<GnssChannel> previous_filter_channels_;
      LunarGnssODTSSummary summary_;
      double pos_err2_sum_ = 0.0;
      double vel_err2_sum_ = 0.0;
      int next_epoch_index_ = 0;
    };

    LunarGnssODTSSummary RunOneMonteCarlo(const LunarGnssODTSConfig& cfg, int mc_index, Real t0_tdb,
                                          const VecXd& elapsed_s, const VecXd& times_tdb,
                                          const VecXd& receiver_clock_s,
                                          std::vector<GNSSMeasurementsEpoch> precomputed) {
      std::vector<State> truth_states
          = BuildTruthStates(cfg, t0_tdb, elapsed_s, times_tdb, true, cfg.seed + 1000 + mc_index);

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
    }

    void ApplyDesignNode(const YAML::Node& node, LunarGnssODTSConfig& cfg) {
      if (!node) return;
      ReadReceiverParams(node["receiver_tracking"], cfg.design.receiver_params);

      const YAML::Node link = node["link_budget"];
      cfg.design.cn0_threshold_dbhz
          = ReadYaml(link, "cn0_threshold_dbhz", cfg.design.cn0_threshold_dbhz);
      cfg.design.apply_cn0_threshold
          = ReadYaml(link, "apply_cn0_threshold", cfg.design.apply_cn0_threshold);
      cfg.design.setup_transmitters
          = ReadYaml(link, "setup_transmitters", cfg.design.setup_transmitters);
      cfg.design.use_cn0_measurement_sigmas
          = ReadYaml(link, "use_cn0_measurement_sigmas", cfg.design.use_cn0_measurement_sigmas);

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

  LunarGnssODTSConfig LoadLunarGnssODTSConfig(const std::filesystem::path& path) {
    YAML::Node root = YAML::LoadFile(path.string());
    std::filesystem::path config_dir
        = path.has_parent_path() ? path.parent_path() : std::filesystem::current_path();
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

    const YAML::Node plasma = root["plasma"];
    cfg.plasma.simulate_truth = ReadYaml(plasma, "simulate_truth", cfg.plasma.simulate_truth);
    cfg.plasma.model_in_filter = ReadYaml(plasma, "model_in_filter", cfg.plasma.model_in_filter);
    const YAML::Node raytrace = plasma["raytrace"];
    cfg.plasma.raytrace_step_size_km
        = ReadYaml(raytrace, "step_size_km", cfg.plasma.raytrace_step_size_km);
    cfg.plasma.raytrace_correction
        = ReadYaml(raytrace, "correction", cfg.plasma.raytrace_correction);
    cfg.plasma.raytrace_fine_correction
        = ReadYaml(raytrace, "fine_correction", cfg.plasma.raytrace_fine_correction);
    cfg.plasma.raytrace_straight_ray
        = ReadYaml(raytrace, "straight_ray", cfg.plasma.raytrace_straight_ray);
    cfg.plasma.raytrace_compute_higher_order
        = ReadYaml(raytrace, "compute_higher_order", cfg.plasma.raytrace_compute_higher_order);
    cfg.plasma.raytrace_use_adaptive_step
        = ReadYaml(raytrace, "use_adaptive_step", cfg.plasma.raytrace_use_adaptive_step);
    cfg.plasma.raytrace_use_fortran_gcpm
        = ReadYaml(raytrace, "use_fortran_gcpm", cfg.plasma.raytrace_use_fortran_gcpm);
    cfg.plasma.raytrace_cutoff_radius_re
        = ReadYaml(raytrace, "cutoff_radius_re", cfg.plasma.raytrace_cutoff_radius_re);
    cfg.plasma.raytrace_gradient_step_km
        = ReadYaml(raytrace, "gradient_step_km", cfg.plasma.raytrace_gradient_step_km);
    cfg.plasma.raytrace_correction_tolerance_m
        = ReadYaml(raytrace, "correction_tolerance_m", cfg.plasma.raytrace_correction_tolerance_m);
    cfg.plasma.raytrace_kp = ReadYaml(raytrace, "kp", cfg.plasma.raytrace_kp);
    cfg.plasma.raytrace_integrator
        = ReadYaml(raytrace, "integrator", cfg.plasma.raytrace_integrator);
    cfg.plasma.raytrace_correction_method
        = ReadYaml(raytrace, "correction_method", cfg.plasma.raytrace_correction_method);
    cfg.plasma.filter_pseudorange_noise_inflation_m
        = ReadYaml(plasma, "filter_pseudorange_noise_inflation_m",
                   cfg.plasma.filter_pseudorange_noise_inflation_m);
    cfg.plasma.filter_doppler_noise_inflation_hz = ReadYaml(
        plasma, "filter_doppler_noise_inflation_hz", cfg.plasma.filter_doppler_noise_inflation_hz);

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
    cfg.process_clock_bias_sigma_s_sqrt_s = ReadYaml(filter, "process_clock_bias_sigma_s_sqrt_s",
                                                     cfg.process_clock_bias_sigma_s_sqrt_s);
    cfg.process_clock_drift_sigma_sps_sqrt_s = ReadYaml(
        filter, "process_clock_drift_sigma_sps_sqrt_s", cfg.process_clock_drift_sigma_sps_sqrt_s);
    cfg.process_srp_coeff_sigma_m2_kg_sqrt_s = ReadYaml(
        filter, "process_srp_coeff_sigma_m2_kg_sqrt_s", cfg.process_srp_coeff_sigma_m2_kg_sqrt_s);
    cfg.integration_step_s = ReadYaml(filter, "integration_step_s", cfg.integration_step_s);

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

    if (cfg.constellation.auto_select_sp3) {
      const auto [t_start_tai, t_end_tai] = EphemerisWindowTai(cfg);
      cfg.constellation.sp3_files
          = SelectSp3FilesForEpochWindow(cfg.constellation.sp3_directory, t_start_tai, t_end_tai);
    }
    return cfg;
  }

  std::vector<LunarGnssODTSSummary> RunLunarGnssODTSMonteCarlo(const LunarGnssODTSConfig& cfg) {
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

    std::vector<State> nominal_truth_states
        = BuildTruthStates(cfg, t0_tdb, elapsed_s, times_tdb, false, cfg.seed);
    std::vector<RuntimeConstellation> constellations = BuildConstellations(cfg, ephem_times_tai);
    std::vector<Real> receive_times = ToRealVector(times_tdb);
    LunarGnssODTSConfig link_cfg = cfg;
    link_cfg.plasma.simulate_truth = false;
    std::vector<GNSSMeasurementsEpoch> precomputed = PrecomputeConstellationChannels(
        link_cfg, constellations, receive_times, nominal_truth_states);
    if (cfg.plasma.simulate_truth) ApplyDelayFile(precomputed, cfg.delays_file);

    std::vector<LunarGnssODTSSummary> summaries;
    summaries.reserve(cfg.monte_carlo_runs);
    for (int mc = 0; mc < cfg.monte_carlo_runs; ++mc) {
      summaries.push_back(
          RunOneMonteCarlo(cfg, mc, t0_tdb, elapsed_s, times_tdb, receiver_clock_s, precomputed));
    }
    WriteSummary(cfg, summaries);
    return summaries;
  }

  void PrecomputeLunarGnssODTSLinks(const LunarGnssODTSConfig& cfg) {
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
        = BuildTruthStates(cfg, t0_tdb, elapsed_s, times_tdb, false, cfg.seed);
    std::vector<RuntimeConstellation> constellations = BuildConstellations(cfg, ephem_times_tai);
    std::vector<Real> receive_times = ToRealVector(times_tdb);
    LunarGnssODTSConfig link_cfg = cfg;
    link_cfg.plasma.simulate_truth = false;
    std::vector<GNSSMeasurementsEpoch> precomputed = PrecomputeConstellationChannels(
        link_cfg, constellations, receive_times, nominal_truth_states);
    WriteLinksCsv(cfg, precomputed, nominal_truth_states);
    WriteLinkCacheMetadata(cfg, fingerprint);
  }

  LunarGnssODTSSimulation::LunarGnssODTSSimulation(std::filesystem::path config_path)
      : config_path_(std::move(config_path)) {
    Setup();
  }

  LunarGnssODTSSimulation::LunarGnssODTSSimulation(LunarGnssODTSConfig config)
      : config_(std::move(config)) {
    Setup();
  }

  void LunarGnssODTSSimulation::Setup() {
    if (setup_complete_) return;
    if (!config_path_.empty()) config_ = LoadLunarGnssODTSConfig(config_path_);
    setup_complete_ = true;
  }

  void LunarGnssODTSSimulation::Precompute() {
    Setup();
    PrecomputeLunarGnssODTSLinks(config_);
  }

  void LunarGnssODTSSimulation::Run() {
    Setup();
    summaries_ = RunLunarGnssODTSMonteCarlo(config_);
  }

}  // namespace lupnt

#include "lupnt/applications/lunar_sat_odts/isl_odts_coordinator_app.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

#include "lupnt/agents/agent.h"
#include "lupnt/applications/ephemeris/lunanet_sat_app.h"
#include "lupnt/applications/lunar_sat_odts/isl_odts_app.h"
#include "lupnt/conversions/coordinate_conversions.h"
#include "lupnt/conversions/frame_converter.h"
#include "lupnt/core/asset_factory.h"
#include "lupnt/lupnt.h"
#include "lupnt/measurements/measurement.h"
#include "lupnt/measurements/measurement_utils.h"
#include "lupnt/numerics/filters/ekf.h"
#include "lupnt/numerics/filters/filter_utils.h"
#include "lupnt/simulations/simulation.h"

namespace lupnt {

  namespace {

    Ptr<NBodyDynamics> MakeOrbitDynamics(int moon_degree, int moon_order, bool use_autodiff,
                                         const IslOdtsConfig& cfg) {
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
      return dynamics;
    }

    Ptr<JointOrbitClockDynamics> MakeJointDynamics(Ptr<NumericalDynamics> orbit_dynamics,
                                                   bool add_clock_noise, int seed,
                                                   const IslOdtsConfig& cfg) {
      auto clock = MakePtr<ClockDynamics>();
      clock->SetModel(ClockModel::OCXO);
      clock->SetClockBiasUnit(ClockBiasUnit::SECONDS);
      clock->SetAddNoise(add_clock_noise);
      clock->SetSeed(seed);

      auto joint = MakePtr<JointOrbitClockDynamics>();
      joint->SetOrbitDynamics(orbit_dynamics);
      joint->SetClockDynamics(clock);
      joint->SetUseClockRelativity(true);
      joint->SetRelativityCenterBody(BodyId::MOON);
      joint->SetAddClockNoise(add_clock_noise);
      joint->SetFrame(Frame::MOON_CI);
      joint->SetTimeStep(cfg.integration_step_s);
      joint->SetIntegrator(IntegratorType::RKF45);
      joint->SetIntegratorParams(IntegratorParams(20, 1.0e-12, 1.0e-12));
      return joint;
    }

    State MakeInitialTruthState(const IslOdtsSatelliteConfig& sat) {
      Cart6 orbit0(sat.r0_m.cast<Real>(), sat.v0_mps.cast<Real>(), Frame::MOON_CI);
      ClockState2 clock0;
      clock0.b() = sat.clock_bias_s;
      clock0.d() = sat.clock_drift_sps;
      return State(JointOrbitClockState(orbit0, clock0));
    }

    VecXd SampleErrorPosVelClock(double sigma_r, double sigma_v, double sigma_b, double sigma_d) {
      VecXd e(8);
      e(0) = SampleNormal(0.0, sigma_r).val();
      e(1) = SampleNormal(0.0, sigma_r).val();
      e(2) = SampleNormal(0.0, sigma_r).val();
      e(3) = SampleNormal(0.0, sigma_v).val();
      e(4) = SampleNormal(0.0, sigma_v).val();
      e(5) = SampleNormal(0.0, sigma_v).val();
      e(6) = SampleNormal(0.0, sigma_b).val();
      e(7) = SampleNormal(0.0, sigma_d).val();
      return e;
    }

    double ComputeLinkBudgetCn0Dbhz(const IslLinkBudgetConfig& lb, double range_m) {
      if (!lb.enabled || !(range_m > 0.0)) return std::numeric_limits<double>::quiet_NaN();
      const double fspl_db = 20.0 * std::log10(4.0 * PI * range_m * lb.frequency_hz / C);
      const double received_power_dbw = lb.tx_power_dbw + lb.tx_gain_dbi + lb.rx_gain_dbi - fspl_db;
      constexpr double kBoltzmannDbwHzK = -228.6;  // 10*log10(1.380649e-23 J/K)
      const double noise_density_dbw_hz
          = kBoltzmannDbwHzK + 10.0 * std::log10(lb.system_noise_temp_k);
      return received_power_dbw - noise_density_dbw_hz;
    }

    Vec3d ParseVec3(const Config& node) {
      return Vec3d(node[0].as<double>(), node[1].as<double>(), node[2].as<double>());
    }

  }  // namespace

  // --- IslOdtsConfig <-> Config (YAML) translation -----------------------------

  IslOdtsConfig ConfigToIslOdtsConfig(Config& config) {
    IslOdtsConfig c;
    c.seed = config["seed"].as<int>(c.seed);
    c.start_epoch_utc = config["start_epoch_utc"].as<std::string>(c.start_epoch_utc);
    c.duration_s = config["duration_s"].as<double>(c.duration_s);
    c.dt_s = config["dt_s"].as<double>(c.dt_s);
    c.integration_step_s = config["integration_step_s"].as<double>(c.integration_step_s);

    LUPNT_CHECK(config["satellites"], "IslOdtsCoordinatorApp requires a `satellites` list",
                "IslOdtsCoordinatorApp");
    for (const auto& item : config["satellites"]) {
      Config s(item);
      IslOdtsSatelliteConfig sat;
      sat.name = s["name"].as<std::string>(std::string("SV"));
      sat.r0_m = ParseVec3(Config(s["r0_m"]));
      sat.v0_mps = ParseVec3(Config(s["v0_mps"]));
      sat.clock_bias_s = s["clock_bias_s"].as<double>(0.0);
      sat.clock_drift_sps = s["clock_drift_sps"].as<double>(0.0);
      c.satellites.push_back(sat);
    }

    if (config["surface_stations"]) {
      for (const auto& item : config["surface_stations"]) {
        Config s(item);
        IslSurfaceStationConfig st;
        st.enabled = s["enabled"].as<bool>(true);
        st.name = s["name"].as<std::string>(st.name);
        st.latitude_deg = s["latitude_deg"].as<double>(st.latitude_deg);
        st.longitude_deg = s["longitude_deg"].as<double>(st.longitude_deg);
        st.altitude_m = s["altitude_m"].as<double>(st.altitude_m);
        st.pseudorange_sigma_m = s["pseudorange_sigma_m"].as<double>(st.pseudorange_sigma_m);
        st.elevation_mask_deg = s["elevation_mask_deg"].as<double>(st.elevation_mask_deg);
        c.surface_stations.push_back(st);
      }
    }

    // Measurement suite / tuning.
    c.range_sigma_m = config["range_sigma_m"].as<double>(c.range_sigma_m);
    c.range_rate_sigma_mps = config["range_rate_sigma_mps"].as<double>(c.range_rate_sigma_mps);
    c.enable_two_way_time_transfer
        = config["enable_two_way_time_transfer"].as<bool>(c.enable_two_way_time_transfer);
    c.time_transfer_sigma_m = config["time_transfer_sigma_m"].as<double>(c.time_transfer_sigma_m);
    c.enable_two_way_frequency_transfer
        = config["enable_two_way_frequency_transfer"].as<bool>(c.enable_two_way_frequency_transfer);
    c.frequency_transfer_sigma_mps
        = config["frequency_transfer_sigma_mps"].as<double>(c.frequency_transfer_sigma_mps);
    c.enable_station_doppler = config["enable_station_doppler"].as<bool>(c.enable_station_doppler);
    c.station_doppler_sigma_mps
        = config["station_doppler_sigma_mps"].as<double>(c.station_doppler_sigma_mps);
    c.enable_centralized_ground_filter
        = config["enable_centralized_ground_filter"].as<bool>(c.enable_centralized_ground_filter);
    c.central_process_accel_sigma_mps2
        = config["central_process_accel_sigma_mps2"].as<double>(c.central_process_accel_sigma_mps2);
    c.central_outlier_threshold
        = config["central_outlier_threshold"].as<double>(c.central_outlier_threshold);
    c.consider_exchange_interval_s
        = config["consider_exchange_interval_s"].as<double>(c.consider_exchange_interval_s);
    c.exchange_use_covariance_intersection
        = config["exchange_use_covariance_intersection"].as<bool>(
            c.exchange_use_covariance_intersection);
    c.exchange_ci_weight = config["exchange_ci_weight"].as<double>(c.exchange_ci_weight);

    c.moon_gravity_degree_truth
        = config["moon_gravity_degree_truth"].as<int>(c.moon_gravity_degree_truth);
    c.moon_gravity_order_truth
        = config["moon_gravity_order_truth"].as<int>(c.moon_gravity_order_truth);
    c.moon_gravity_degree_filter
        = config["moon_gravity_degree_filter"].as<int>(c.moon_gravity_degree_filter);
    c.moon_gravity_order_filter
        = config["moon_gravity_order_filter"].as<int>(c.moon_gravity_order_filter);
    c.include_earth = config["include_earth"].as<bool>(c.include_earth);
    c.include_sun = config["include_sun"].as<bool>(c.include_sun);
    c.use_relativity = config["use_relativity"].as<bool>(c.use_relativity);

    c.initial_position_sigma_m
        = config["initial_position_sigma_m"].as<double>(c.initial_position_sigma_m);
    c.initial_velocity_sigma_mps
        = config["initial_velocity_sigma_mps"].as<double>(c.initial_velocity_sigma_mps);
    c.initial_clock_bias_sigma_s
        = config["initial_clock_bias_sigma_s"].as<double>(c.initial_clock_bias_sigma_s);
    c.initial_clock_drift_sigma_sps
        = config["initial_clock_drift_sigma_sps"].as<double>(c.initial_clock_drift_sigma_sps);
    c.consider_position_sigma_m
        = config["consider_position_sigma_m"].as<double>(c.consider_position_sigma_m);
    c.consider_velocity_sigma_mps
        = config["consider_velocity_sigma_mps"].as<double>(c.consider_velocity_sigma_mps);
    c.consider_clock_bias_sigma_s
        = config["consider_clock_bias_sigma_s"].as<double>(c.consider_clock_bias_sigma_s);
    c.consider_clock_drift_sigma_sps
        = config["consider_clock_drift_sigma_sps"].as<double>(c.consider_clock_drift_sigma_sps);
    c.process_accel_sigma_mps2
        = config["process_accel_sigma_mps2"].as<double>(c.process_accel_sigma_mps2);
    return c;
  }

  Config IslOdtsConfigToConfig(const IslOdtsConfig& cfg) {
    Config c;
    c["seed"] = cfg.seed;
    c["start_epoch_utc"] = cfg.start_epoch_utc;
    c["duration_s"] = cfg.duration_s;
    c["dt_s"] = cfg.dt_s;
    c["integration_step_s"] = cfg.integration_step_s;

    for (const IslOdtsSatelliteConfig& s : cfg.satellites) {
      Config node;
      node["name"] = s.name;
      Config r0, v0;
      for (int i = 0; i < 3; ++i) r0.push_back(s.r0_m(i));
      for (int i = 0; i < 3; ++i) v0.push_back(s.v0_mps(i));
      node["r0_m"] = r0;
      node["v0_mps"] = v0;
      node["clock_bias_s"] = s.clock_bias_s;
      node["clock_drift_sps"] = s.clock_drift_sps;
      c["satellites"].push_back(node);
    }

    std::vector<IslSurfaceStationConfig> stations = cfg.surface_stations;
    if (stations.empty() && cfg.surface_station.enabled) stations.push_back(cfg.surface_station);
    for (const IslSurfaceStationConfig& st : stations) {
      Config node;
      node["enabled"] = st.enabled;
      node["name"] = st.name;
      node["latitude_deg"] = st.latitude_deg;
      node["longitude_deg"] = st.longitude_deg;
      node["altitude_m"] = st.altitude_m;
      node["pseudorange_sigma_m"] = st.pseudorange_sigma_m;
      node["elevation_mask_deg"] = st.elevation_mask_deg;
      c["surface_stations"].push_back(node);
    }

    c["range_sigma_m"] = cfg.range_sigma_m;
    c["range_rate_sigma_mps"] = cfg.range_rate_sigma_mps;
    c["enable_two_way_time_transfer"] = cfg.enable_two_way_time_transfer;
    c["time_transfer_sigma_m"] = cfg.time_transfer_sigma_m;
    c["enable_two_way_frequency_transfer"] = cfg.enable_two_way_frequency_transfer;
    c["frequency_transfer_sigma_mps"] = cfg.frequency_transfer_sigma_mps;
    c["enable_station_doppler"] = cfg.enable_station_doppler;
    c["station_doppler_sigma_mps"] = cfg.station_doppler_sigma_mps;
    c["enable_centralized_ground_filter"] = cfg.enable_centralized_ground_filter;
    c["central_process_accel_sigma_mps2"] = cfg.central_process_accel_sigma_mps2;
    c["central_outlier_threshold"] = cfg.central_outlier_threshold;
    c["consider_exchange_interval_s"] = cfg.consider_exchange_interval_s;
    c["exchange_use_covariance_intersection"] = cfg.exchange_use_covariance_intersection;
    c["exchange_ci_weight"] = cfg.exchange_ci_weight;
    c["moon_gravity_degree_truth"] = cfg.moon_gravity_degree_truth;
    c["moon_gravity_order_truth"] = cfg.moon_gravity_order_truth;
    c["moon_gravity_degree_filter"] = cfg.moon_gravity_degree_filter;
    c["moon_gravity_order_filter"] = cfg.moon_gravity_order_filter;
    c["include_earth"] = cfg.include_earth;
    c["include_sun"] = cfg.include_sun;
    c["use_relativity"] = cfg.use_relativity;
    c["initial_position_sigma_m"] = cfg.initial_position_sigma_m;
    c["initial_velocity_sigma_mps"] = cfg.initial_velocity_sigma_mps;
    c["initial_clock_bias_sigma_s"] = cfg.initial_clock_bias_sigma_s;
    c["initial_clock_drift_sigma_sps"] = cfg.initial_clock_drift_sigma_sps;
    c["consider_position_sigma_m"] = cfg.consider_position_sigma_m;
    c["consider_velocity_sigma_mps"] = cfg.consider_velocity_sigma_mps;
    c["consider_clock_bias_sigma_s"] = cfg.consider_clock_bias_sigma_s;
    c["consider_clock_drift_sigma_sps"] = cfg.consider_clock_drift_sigma_sps;
    c["process_accel_sigma_mps2"] = cfg.process_accel_sigma_mps2;
    return c;
  }

  // --- IslOdtsCoordinatorApp ----------------------------------------------------

  IslOdtsCoordinatorApp::IslOdtsCoordinatorApp(Config& config) : Application(config) {
    cfg_ = ConfigToIslOdtsConfig(config);
  }

  IslOdtsCoordinatorApp::IslOdtsCoordinatorApp(IslOdtsConfig config) : cfg_(std::move(config)) {}

  void IslOdtsCoordinatorApp::Setup() {
    LUPNT_CHECK(agent_, "Agent not set", "IslOdtsCoordinatorApp");
    LUPNT_CHECK(cfg_.satellites.size() >= 2,
                "IslOdtsConfig.satellites must contain at least 2 satellites",
                "IslOdtsCoordinatorApp");
    // Defer the heavy precompute to the first Step (so a programmatically-set config is
    // honored). Here we only schedule the per-epoch Steps at the measurement cadence.
    Simulation* sim = agent_->GetSimulation();
    sim->Schedule(0.0, [this](Real t) { Step(t); }, 1.0 / cfg_.dt_s, Event::Priority::APPLICATION);
  }

  void IslOdtsCoordinatorApp::Initialize() {
    const IslOdtsConfig& cfg = cfg_;
    LUPNT_CHECK(cfg.duration_s > 0.0, "duration_s must be positive", "IslOdtsCoordinatorApp");
    LUPNT_CHECK(cfg.dt_s > 0.0, "dt_s must be positive", "IslOdtsCoordinatorApp");
    LUPNT_CHECK(cfg.integration_step_s > 0.0, "integration_step_s must be positive",
                "IslOdtsCoordinatorApp");
    LUPNT_CHECK(cfg.range_sigma_m > 0.0, "range_sigma_m must be positive", "IslOdtsCoordinatorApp");
    LUPNT_CHECK(cfg.range_rate_sigma_mps > 0.0, "range_rate_sigma_mps must be positive",
                "IslOdtsCoordinatorApp");

    n_sat_ = static_cast<int>(cfg.satellites.size());
    n_links_ = n_sat_ - 1;
    n_state_ = kSubStateSize * n_sat_;
    const int n_sat = n_sat_;
    const int n_links = n_links_;
    const int n_state = n_state_;
    const double nan = std::numeric_limits<double>::quiet_NaN();

    RandomEngine::SetSeed(static_cast<unsigned int>(cfg.seed));

    // The coordinator time-keeps in ABSOLUTE TDB seconds past J2000 (t0_tdb_ + elapsed),
    // passing those directly to the orbit dynamics and frame converters -- exactly as the
    // former monolithic driver did. Orbit propagation adds `GetLupntEpoch()` to its time
    // argument (numerical_orbit_dynamics.cc: `t_tdb = t + GetLupntEpoch()`), so we reset the
    // global epoch to 0 here (the base `Simulation` sets it from the scenario `epoch:`) to
    // avoid double-counting the epoch and to reproduce the monolith's numerics bit-for-bit.
    SetLupntEpoch(0.0);

    t0_tdb_ = ConvertTime(GregorianToTime(cfg.start_epoch_utc), Time::UTC, Time::TDB);
    N_ = static_cast<int>(std::floor(cfg.duration_s / cfg.dt_s + 1.0e-9)) + 1;
    const int N = N_;

    VecXd t_s(N);
    for (int k = 0; k < N; ++k) t_s(k) = k * cfg.dt_s;

    // --- Truth dynamics: independent per satellite so stochastic clock-noise draws
    // don't interleave between satellites.
    dyn_truth_.assign(n_sat, nullptr);
    for (int j = 0; j < n_sat; ++j) {
      Ptr<NBodyDynamics> orbit_truth = MakeOrbitDynamics(cfg.moon_gravity_degree_truth,
                                                         cfg.moon_gravity_order_truth, false, cfg);
      dyn_truth_[j] = MakeJointDynamics(orbit_truth, true, cfg.seed + 1000 * (j + 1), cfg);
    }

    x_truth_.assign(n_sat, State());
    x_truth0_.assign(n_sat, VecXd());
    for (int j = 0; j < n_sat; ++j) {
      x_truth_[j] = MakeInitialTruthState(cfg.satellites[j]);
      x_truth0_[j] = x_truth_[j].cast<double>();
    }

    // --- Per-filter state ordering.
    order_.assign(n_sat, {});
    pos_.assign(n_sat, std::vector<int>(n_sat, -1));
    for (int j = 0; j < n_sat; ++j) {
      order_[j].push_back(j);
      for (int g = 0; g < n_sat; ++g)
        if (g != j) order_[j].push_back(g);
      for (int b = 0; b < n_sat; ++b) pos_[j][order_[j][b]] = b;
    }

    // --- N parallel onboard Schmidt-EKF filters, one per satellite, each hosted like
    // flight software inside its own LunaNetSatApp.
    apps_.assign(n_sat, nullptr);
    sat_apps_.assign(n_sat, nullptr);
    for (int j = 0; j < n_sat; ++j) {
      VecXd x0j(n_state);
      MatXd P0j = MatXd::Zero(n_state, n_state);
      for (int b = 0; b < n_sat; ++b) {
        const int g = order_[j][b];
        const double sr = (b == 0) ? cfg.initial_position_sigma_m : cfg.consider_position_sigma_m;
        const double sv
            = (b == 0) ? cfg.initial_velocity_sigma_mps : cfg.consider_velocity_sigma_mps;
        const double sb
            = (b == 0) ? cfg.initial_clock_bias_sigma_s : cfg.consider_clock_bias_sigma_s;
        const double sd
            = (b == 0) ? cfg.initial_clock_drift_sigma_sps : cfg.consider_clock_drift_sigma_sps;
        x0j.segment(kSubStateSize * b, kSubStateSize)
            = x_truth0_[g] + SampleErrorPosVelClock(sr, sv, sb, sd);
        P0j.block(kSubStateSize * b, kSubStateSize * b, kSubStateSize, kSubStateSize)
            = InitialCovariancePosVelClock(sr, sv, sb, sd);
      }

      IslOdtsAppParams app_params;
      app_params.n_sat = n_sat;
      app_params.range_sigma_m = cfg.range_sigma_m;
      app_params.range_rate_sigma_mps = cfg.range_rate_sigma_mps;
      app_params.pseudorange_sigma_m = cfg.surface_stations.empty()
                                           ? cfg.surface_station.pseudorange_sigma_m
                                           : cfg.surface_stations.front().pseudorange_sigma_m;
      app_params.process_accel_sigma_mps2 = cfg.process_accel_sigma_mps2;
      app_params.include_time_transfer = cfg.enable_two_way_time_transfer;
      app_params.time_transfer_sigma_m = cfg.time_transfer_sigma_m;
      app_params.include_frequency_transfer = cfg.enable_two_way_frequency_transfer;
      app_params.frequency_transfer_sigma_mps = cfg.frequency_transfer_sigma_mps;
      app_params.include_anchor_doppler = cfg.enable_station_doppler;
      app_params.anchor_doppler_sigma_mps = cfg.station_doppler_sigma_mps;
      apps_[j] = MakePtr<IslOdtsApp>(app_params);

      Ptr<NBodyDynamics> orbit_filter = MakeOrbitDynamics(cfg.moon_gravity_degree_filter,
                                                          cfg.moon_gravity_order_filter, true, cfg);
      Ptr<JointOrbitClockDynamics> dyn_filter = MakeJointDynamics(orbit_filter, false, 0, cfg);
      apps_[j]->Configure(t0_tdb_, x0j, P0j, dyn_filter);

      sat_apps_[j] = MakePtr<LunaNetSatApp>();
      sat_apps_[j]->SetName("lunanet_sat_" + cfg.satellites[j].name);
      sat_apps_[j]->AddSubApp(apps_[j]);
      sat_apps_[j]->Setup();
    }

    // --- Lunar surface station beacon network (fixed in the Moon principal-axis frame).
    stations_ = cfg.surface_stations;
    if (stations_.empty() && cfg.surface_station.enabled) stations_.push_back(cfg.surface_station);
    n_stn_ = static_cast<int>(stations_.size());
    const int n_stn = n_stn_;

    station_bf_.assign(n_stn, Vec3::Zero());
    for (int s = 0; s < n_stn; ++s) {
      State lla(3);
      VecX lla_v(3);
      lla_v << stations_[s].latitude_deg, stations_[s].longitude_deg, stations_[s].altitude_m;
      lla = lla_v;
      State station_cart = LatLonAltToCart(lla, R_MOON, 0.0);
      station_bf_[s] = Vec3(station_cart.head(3));
    }

    // --- Result buffers.
    results_ = IslOdtsResults{};
    results_.t_s = t_s;
    results_.satellite_names.resize(n_sat);
    for (int j = 0; j < n_sat; ++j) results_.satellite_names[j] = cfg.satellites[j].name;
    results_.truth_states.assign(n_sat, MatXd(N, kSubStateSize));
    results_.est.assign(n_sat, MatXd::Zero(N, n_state));
    results_.cov_diag.assign(n_sat, MatXd::Zero(N, n_state));
    results_.cov_own_full.assign(n_sat, MatXd::Zero(N, kSubStateSize * kSubStateSize));
    results_.range_true_m.resize(N, n_links);
    results_.range_rate_true_mps.resize(N, n_links);
    results_.range_obs_m = MatXd::Constant(N, n_links, nan);
    results_.range_rate_obs_mps = MatXd::Constant(N, n_links, nan);
    results_.range_resid_m.assign(n_sat, MatXd::Constant(N, n_links, nan));
    results_.cn0_dbhz.resize(N, n_links);
    results_.time_transfer_true_m = MatXd::Constant(N, n_links, nan);
    results_.time_transfer_obs_m = MatXd::Constant(N, n_links, nan);
    results_.station_pos_mci.assign(n_stn, MatXd::Constant(N, 3, nan));
    results_.station_visible = MatXd::Zero(N, n_sat);
    results_.station_pr_true_m = MatXd::Constant(N, n_sat, nan);
    results_.station_pr_obs_m = MatXd::Constant(N, n_sat, nan);
    results_.station_pr_resid_m = MatXd::Constant(N, n_sat, nan);
    results_.est_central = MatXd::Zero(N, n_state);
    results_.cov_central_full.assign(n_sat, MatXd::Zero(N, kSubStateSize * kSubStateSize));

    // --- Centralized ground filter.
    central_ = nullptr;
    if (cfg.enable_centralized_ground_filter && n_stn > 0) {
      Ptr<NBodyDynamics> orbit_c = MakeOrbitDynamics(cfg.moon_gravity_degree_filter,
                                                     cfg.moon_gravity_order_filter, true, cfg);
      Ptr<JointOrbitClockDynamics> dyn_c = MakeJointDynamics(orbit_c, false, 0, cfg);
      central_ = MakePtr<EKF>();
      VecXd x0c(n_state);
      MatXd P0c = MatXd::Zero(n_state, n_state);
      for (int j = 0; j < n_sat; ++j) {
        x0c.segment(kSubStateSize * j, kSubStateSize)
            = x_truth0_[j]
              + SampleErrorPosVelClock(cfg.initial_position_sigma_m, cfg.initial_velocity_sigma_mps,
                                       cfg.initial_clock_bias_sigma_s,
                                       cfg.initial_clock_drift_sigma_sps);
        P0c.block(kSubStateSize * j, kSubStateSize * j, kSubStateSize, kSubStateSize)
            = InitialCovariancePosVelClock(
                cfg.initial_position_sigma_m, cfg.initial_velocity_sigma_mps,
                cfg.initial_clock_bias_sigma_s, cfg.initial_clock_drift_sigma_sps);
      }
      State x0c_state(n_state);
      x0c_state = x0c.cast<Real>();
      x0c_state.SetFrame(Frame::MOON_CI);
      central_->SetTime(t0_tdb_);
      central_->SetState(x0c_state);
      central_->SetCovariance(P0c);
      central_->SetDynamicsFunction(
          [dyn_c, n_sat](const State& x, Real t0, Real tf, const State*, MatXd* F) -> State {
            const int n = n_sat * kSubStateSize;
            State xf(n);
            xf.SetFrame(x.GetFrame());
            if (F != nullptr) F->setZero(n, n);
            for (int j = 0; j < n_sat; ++j) {
              const int off = j * kSubStateSize;
              State sub(kSubStateSize);
              sub = x.segment(off, kSubStateSize);
              sub.SetFrame(x.GetFrame());
              JointOrbitClockState xj(sub);
              if (F != nullptr) {
                MatXd Fj;
                State xfj = dyn_c->Propagate(xj, t0, tf, nullptr, &Fj);
                xf.segment(off, kSubStateSize) = xfj;
                F->block(off, off, kSubStateSize, kSubStateSize) = Fj;
              } else {
                State xfj = dyn_c->Propagate(xj, t0, tf, nullptr);
                xf.segment(off, kSubStateSize) = xfj;
              }
            }
            return xf;
          });
      const double sigma_a = cfg.central_process_accel_sigma_mps2;
      central_->SetProcessNoiseFunction([n_sat, sigma_a](const State& x, Real t0, Real tf) {
        const double dt = std::abs((tf - t0).val());
        MatXd Qsub = MatXd::Zero(kSubStateSize, kSubStateSize);
        Mat3d Qacc = std::pow(sigma_a, 2) * Mat3d::Identity();
        Qsub.block(0, 0, 6, 6) = ProcessNoisePosVel(Qacc, dt);
        Qsub.block(6, 6, 2, 2) = ClockDynamics::TwoStateNoise(ClockModel::OCXO, dt).cast<double>();
        MatXd Q = MatXd::Zero(x.size(), x.size());
        for (int j = 0; j < n_sat; ++j)
          Q.block(j * kSubStateSize, j * kSubStateSize, kSubStateSize, kSubStateSize) = Qsub;
        return Q;
      });
      central_->SetOutlierThreshold(cfg.central_outlier_threshold);
    }

    // --- Epoch 0: record truth + initial estimates (no measurement).
    for (int j = 0; j < n_sat; ++j)
      results_.truth_states[j].row(0) = x_truth_[j].cast<double>().transpose();
    {
      MatXd pr = MatXd::Zero(n_sat, n_sat), prr = MatXd::Zero(n_sat, n_sat);
      for (int a = 0; a < n_sat; ++a)
        for (int b = a + 1; b < n_sat; ++b) {
          Vec2 y
              = RangeAndRangeRate(VecX(x_truth_[a].head(3)), VecX(x_truth_[b].head(3)),
                                  VecX(x_truth_[a].segment(3, 3)), VecX(x_truth_[b].segment(3, 3)));
          pr(a, b) = pr(b, a) = y(0).val();
          prr(a, b) = prr(b, a) = y(1).val();
        }
      for (int i = 0; i < n_links; ++i) {
        results_.range_true_m(0, i) = pr(0, i + 1);
        results_.range_rate_true_mps(0, i) = prr(0, i + 1);
        results_.cn0_dbhz(0, i) = ComputeLinkBudgetCn0Dbhz(cfg.link_budget, pr(0, i + 1));
      }
    }
    RecordFilters(0);
    RecordCentral(0);

    next_exchange_s_ = (cfg.consider_exchange_interval_s > 0.0)
                           ? cfg.consider_exchange_interval_s
                           : std::numeric_limits<double>::infinity();
  }

  void IslOdtsCoordinatorApp::RecordFilters(int k) {
    const int n_sat = n_sat_;
    for (int j = 0; j < n_sat; ++j) {
      State est_state = apps_[j]->GetEstimate();
      VecXd est_d = est_state.cast<double>();
      MatXd P = apps_[j]->GetCovariance();
      VecXd cov_d = P.diagonal();
      for (int g = 0; g < n_sat; ++g) {
        const int b = pos_[j][g];
        results_.est[j].block(k, kSubStateSize * g, 1, kSubStateSize)
            = est_d.segment(kSubStateSize * b, kSubStateSize).transpose();
        results_.cov_diag[j].block(k, kSubStateSize * g, 1, kSubStateSize)
            = cov_d.segment(kSubStateSize * b, kSubStateSize).transpose();
      }
      MatXd P_own = P.block(0, 0, kSubStateSize, kSubStateSize);
      for (int r = 0; r < kSubStateSize; ++r)
        for (int c = 0; c < kSubStateSize; ++c)
          results_.cov_own_full[j](k, kSubStateSize * r + c) = P_own(r, c);
    }
  }

  void IslOdtsCoordinatorApp::RecordCentral(int k) {
    if (!central_) return;
    const int n_sat = n_sat_;
    VecXd xc = central_->GetState().cast<double>();
    results_.est_central.row(k) = xc.transpose();
    MatXd Pc = central_->GetCovariance();
    for (int j = 0; j < n_sat; ++j) {
      MatXd Pj = Pc.block(kSubStateSize * j, kSubStateSize * j, kSubStateSize, kSubStateSize);
      for (int r = 0; r < kSubStateSize; ++r)
        for (int c = 0; c < kSubStateSize; ++c)
          results_.cov_central_full[j](k, kSubStateSize * r + c) = Pj(r, c);
    }
  }

  void IslOdtsCoordinatorApp::Step(Real t) {
    if (!initialized_) {
      Initialize();
      initialized_ = true;
    }
    int k = static_cast<int>(std::lround(t.val() / cfg_.dt_s));
    if (k < 1 || k >= N_) return;
    RunEpoch(k);
  }

  void IslOdtsCoordinatorApp::RunEpoch(int k) {
    const IslOdtsConfig& cfg = cfg_;
    const int n_sat = n_sat_;
    const int n_links = n_links_;
    const int n_state = n_state_;
    const int n_stn = n_stn_;
    const VecXd& t_s = results_.t_s;

    const Real tkm1 = t0_tdb_ + t_s(k - 1);
    const Real tk = t0_tdb_ + t_s(k);

    // 1. Propagate all truth states.
    for (int j = 0; j < n_sat; ++j)
      x_truth_[j] = dyn_truth_[j]->Propagate(JointOrbitClockState(x_truth_[j]), tkm1, tk, nullptr);

    // 2. Pairwise crosslink truth + one noisy observation per unordered pair.
    const bool tt_on = cfg.enable_two_way_time_transfer;
    const bool ft_on = cfg.enable_two_way_frequency_transfer;
    MatXd pr_true = MatXd::Zero(n_sat, n_sat), prr_true = MatXd::Zero(n_sat, n_sat);
    for (int a = 0; a < n_sat; ++a)
      for (int b = a + 1; b < n_sat; ++b) {
        Vec2 y
            = RangeAndRangeRate(VecX(x_truth_[a].head(3)), VecX(x_truth_[b].head(3)),
                                VecX(x_truth_[a].segment(3, 3)), VecX(x_truth_[b].segment(3, 3)));
        pr_true(a, b) = pr_true(b, a) = y(0).val();
        prr_true(a, b) = prr_true(b, a) = y(1).val();
      }
    MatXd pr_obs = MatXd::Zero(n_sat, n_sat), prr_obs = MatXd::Zero(n_sat, n_sat);
    MatXd tt_obs = MatXd::Zero(n_sat, n_sat);
    MatXd ft_obs = MatXd::Zero(n_sat, n_sat);
    for (int a = 0; a < n_sat; ++a) {
      for (int b = a + 1; b < n_sat; ++b) {
        const double r = pr_true(a, b) + SampleNormal(0.0, cfg.range_sigma_m).val();
        const double rr = prr_true(a, b) + SampleNormal(0.0, cfg.range_rate_sigma_mps).val();
        pr_obs(a, b) = pr_obs(b, a) = r;
        prr_obs(a, b) = prr_obs(b, a) = rr;
        if (tt_on) {
          const double db = C * (x_truth_[a].cast<double>()(6) - x_truth_[b].cast<double>()(6));
          const double tt = db + SampleNormal(0.0, cfg.time_transfer_sigma_m).val();
          tt_obs(a, b) = tt;
          tt_obs(b, a) = -tt;
        }
        if (ft_on) {
          const double dd = C * (x_truth_[a].cast<double>()(7) - x_truth_[b].cast<double>()(7));
          const double ft = dd + SampleNormal(0.0, cfg.frequency_transfer_sigma_mps).val();
          ft_obs(a, b) = ft;
          ft_obs(b, a) = -ft;
        }
      }
    }

    // 3. Surface-station beacon(s).
    const bool sd_on = cfg.enable_station_doppler;
    std::vector<std::vector<Vec3d>> sat_anchor_pos(n_sat);
    std::vector<std::vector<Vec3d>> sat_anchor_vel(n_sat);
    std::vector<std::vector<double>> sat_anchor_obs(n_sat);
    std::vector<std::vector<double>> sat_anchor_dopp_obs(n_sat);
    std::vector<std::tuple<Vec3d, Vec3d, int>> central_meas;
    std::vector<double> central_pr_obs, central_dp_obs;
    for (int s = 0; s < n_stn; ++s) {
      Vec6 st6;
      st6 << station_bf_[s], Vec3::Zero();
      Vec6 st_mci = ConvertFrame(tk, st6, Frame::MOON_PA, Frame::MOON_CI);
      const Vec3d rst = Vec3(st_mci.head(3)).cast<double>();
      const Vec3d vst = Vec3(st_mci.tail(3)).cast<double>();
      results_.station_pos_mci[s].row(k) = rst.transpose();

      for (int j = 0; j < n_sat; ++j) {
        Vec6 xt;
        xt << x_truth_[j].head(3), x_truth_[j].segment(3, 3);
        Vec6 xt_bf = ConvertFrame(tk, xt, Frame::MOON_CI, Frame::MOON_PA);
        Cart3 r_sat_bf(Vec3(xt_bf.head(3)), Frame::MOON_PA);
        Cart3 r_gs_bf(station_bf_[s], Frame::MOON_PA);
        State aer = CartToAzElRange(r_sat_bf, r_gs_bf);
        if ((aer(1) * DEG).val() <= stations_[s].elevation_mask_deg) continue;

        const Vec3d r_j = x_truth_[j].cast<double>().head(3);
        const Vec3d v_j = x_truth_[j].cast<double>().segment(3, 3);
        const double b_j = x_truth_[j].cast<double>()(6);
        const double d_j = x_truth_[j].cast<double>()(7);
        const double pr_t = (r_j - rst).norm() + C * b_j;
        const double pr_o = pr_t + SampleNormal(0.0, stations_[s].pseudorange_sigma_m).val();
        const Vec3d u = (r_j - rst) / (r_j - rst).norm();
        const double dp_t = u.dot(v_j - vst) + C * d_j;
        const double dp_o = dp_t + SampleNormal(0.0, cfg.station_doppler_sigma_mps).val();

        sat_anchor_pos[j].push_back(rst);
        sat_anchor_vel[j].push_back(vst);
        sat_anchor_obs[j].push_back(pr_o);
        sat_anchor_dopp_obs[j].push_back(dp_o);
        central_meas.emplace_back(rst, vst, j);
        central_pr_obs.push_back(pr_o);
        central_dp_obs.push_back(dp_o);
        results_.station_visible(k, j) += 1.0;
        if (std::isnan(results_.station_pr_true_m(k, j))) {
          results_.station_pr_true_m(k, j) = pr_t;
          results_.station_pr_obs_m(k, j) = pr_o;
        }
      }
    }

    // 4a. Distributed onboard filters.
    for (int j = 0; j < n_sat; ++j) {
      IslOdtsMeasurementEpoch meas;
      meas.crosslink_range_m.resize(n_links);
      meas.crosslink_range_rate_mps.resize(n_links);
      if (tt_on) meas.crosslink_time_transfer_m.resize(n_links);
      if (ft_on) meas.crosslink_frequency_transfer_mps.resize(n_links);
      for (int i = 0; i < n_links; ++i) {
        const int g = order_[j][i + 1];
        meas.crosslink_range_m(i) = pr_obs(j, g);
        meas.crosslink_range_rate_mps(i) = prr_obs(j, g);
        if (tt_on) meas.crosslink_time_transfer_m(i) = tt_obs(j, g);
        if (ft_on) meas.crosslink_frequency_transfer_mps(i) = ft_obs(j, g);
      }
      const int n_anchor_j = static_cast<int>(sat_anchor_pos[j].size());
      if (n_anchor_j > 0) {
        meas.anchor_pos_mci = sat_anchor_pos[j];
        meas.anchor_pseudorange_m.resize(n_anchor_j);
        for (int a = 0; a < n_anchor_j; ++a) meas.anchor_pseudorange_m(a) = sat_anchor_obs[j][a];
        if (sd_on) {
          meas.anchor_vel_mci = sat_anchor_vel[j];
          meas.anchor_doppler_mps.resize(n_anchor_j);
          for (int a = 0; a < n_anchor_j; ++a)
            meas.anchor_doppler_mps(a) = sat_anchor_dopp_obs[j][a];
        }
      }
      apps_[j]->StageMeasurements(meas);
      sat_apps_[j]->Step(tk);

      VecXd resid = apps_[j]->GetPrefitResidual();
      for (int i = 0; i < n_links; ++i) results_.range_resid_m[j](k, i) = resid(2 * i);
      if (n_anchor_j > 0) results_.station_pr_resid_m(k, j) = resid(2 * n_links);
    }

    // 4b. Centralized ground filter.
    if (central_) {
      central_->Predict(tk);
      const int np = static_cast<int>(central_meas.size());
      if (np > 0) {
        const double sigma_pr = cfg.surface_stations.empty()
                                    ? cfg.surface_station.pseudorange_sigma_m
                                    : cfg.surface_stations.front().pseudorange_sigma_m;
        const double sigma_dp = cfg.station_doppler_sigma_mps;
        const bool dopp_on = sd_on;
        central_->SetMeasurementFunction([central_meas, sigma_pr, sigma_dp, dopp_on](
                                             const State& x, MatXd* H, MatXd* R) -> VecXd {
          const int nn = static_cast<int>(central_meas.size());
          const int mm = dopp_on ? 2 * nn : nn;
          const VecXd xd = x.cast<double>();
          VecXd y(mm);
          if (H != nullptr) H->setZero(mm, x.size());
          if (R != nullptr) *R = MatXd::Zero(mm, mm);
          for (int i = 0; i < nn; ++i) {
            const Vec3d& rst = std::get<0>(central_meas[i]);
            const Vec3d& vst = std::get<1>(central_meas[i]);
            const int j = std::get<2>(central_meas[i]);
            const Vec3d r_j = xd.segment(kSubStateSize * j, 3);
            const Vec3d v_j = xd.segment(kSubStateSize * j + 3, 3);
            const double b_j = xd(kSubStateSize * j + 6);
            const double d_j = xd(kSubStateSize * j + 7);
            const Vec3d dr = r_j - rst;
            const double rng = dr.norm();
            const Vec3d u = dr / rng;
            y(i) = rng + C * b_j;
            if (H != nullptr) {
              H->block(i, kSubStateSize * j, 1, 3) = u.transpose();
              (*H)(i, kSubStateSize* j + 6) = C;
            }
            if (R != nullptr) (*R)(i, i) = sigma_pr * sigma_pr;
            if (dopp_on) {
              const Vec3d dv = v_j - vst;
              y(nn + i) = u.dot(dv) + C * d_j;
              if (H != nullptr) {
                H->block(nn + i, kSubStateSize * j, 1, 3)
                    = ((dv - u * u.dot(dv)) / rng).transpose();
                H->block(nn + i, kSubStateSize * j + 3, 1, 3) = u.transpose();
                (*H)(nn + i, kSubStateSize * j + 7) = C;
              }
              if (R != nullptr) (*R)(nn + i, nn + i) = sigma_dp * sigma_dp;
            }
          }
          return y;
        });
        const int m = dopp_on ? 2 * np : np;
        VecXd cy(m);
        for (int i = 0; i < np; ++i) {
          cy(i) = central_pr_obs[i];
          if (dopp_on) cy(np + i) = central_dp_obs[i];
        }
        central_->Update(cy);
      }
    }

    // 5. Consider-state exchange.
    if (t_s(k) + 1.0e-9 >= next_exchange_s_) {
      std::vector<VecXd> own_mean(n_sat);
      std::vector<MatXd> own_cov(n_sat);
      for (int i = 0; i < n_sat; ++i) {
        own_mean[i] = apps_[i]->GetEstimate().cast<double>().segment(0, kSubStateSize);
        own_cov[i] = apps_[i]->GetCovariance().block(0, 0, kSubStateSize, kSubStateSize);
      }
      for (int j = 0; j < n_sat; ++j) {
        State xj = apps_[j]->GetEstimate();
        MatXd Pj = apps_[j]->GetCovariance();

        if (cfg.exchange_use_covariance_intersection) {
          const VecXd xjd = xj.cast<double>();
          VecXd xj_new = xjd;
          for (int b = 1; b < n_sat; ++b) {
            const int g = order_[j][b];
            const MatXd Pbb
                = Pj.block(kSubStateSize * b, kSubStateSize * b, kSubStateSize, kSubStateSize);
            const VecXd mb = xjd.segment(kSubStateSize * b, kSubStateSize);

            const VecXd s = Pbb.diagonal().cwiseMax(1e-300).cwiseSqrt();
            const VecXd sinv = s.cwiseInverse();
            const MatXd Sinv = sinv.asDiagonal();
            const MatXd Yloc = (Sinv * Pbb * Sinv).inverse();
            const MatXd Ybc = (Sinv * own_cov[g] * Sinv).inverse();
            const VecXd yloc = Yloc * sinv.cwiseProduct(mb);
            const VecXd ybc = Ybc * sinv.cwiseProduct(own_mean[g]);

            auto fuse = [&](double w, MatXd& Pbb_new, VecXd& mb_new) {
              MatXd Yf = w * Yloc + (1.0 - w) * Ybc;
              VecXd yf = w * yloc + (1.0 - w) * ybc;
              MatXd Pt = Yf.inverse();
              Pt = 0.5 * (Pt + Pt.transpose());
              Pbb_new = s.asDiagonal() * Pt * s.asDiagonal();
              mb_new = s.cwiseProduct(Pt * yf);
            };

            MatXd Pbb_best;
            VecXd mb_best;
            if (cfg.exchange_ci_weight >= 0.0 && cfg.exchange_ci_weight <= 1.0) {
              fuse(cfg.exchange_ci_weight, Pbb_best, mb_best);
            } else {
              const int n_grid = 19;
              double best_obj = std::numeric_limits<double>::infinity();
              for (int gi = 1; gi <= n_grid; ++gi) {
                const double w = static_cast<double>(gi) / (n_grid + 1);
                MatXd Pbb_w;
                VecXd mb_w;
                fuse(w, Pbb_w, mb_w);
                const double obj = (Sinv * Pbb_w * Sinv).trace();
                if (obj < best_obj) {
                  best_obj = obj;
                  Pbb_best = Pbb_w;
                  mb_best = mb_w;
                }
              }
            }

            xj_new.segment(kSubStateSize * b, kSubStateSize) = mb_best;
            Pj.block(kSubStateSize * b, 0, kSubStateSize, n_state).setZero();
            Pj.block(0, kSubStateSize * b, n_state, kSubStateSize).setZero();
            Pj.block(kSubStateSize * b, kSubStateSize * b, kSubStateSize, kSubStateSize) = Pbb_best;
          }
          xj = xj_new.cast<Real>();
        } else {
          for (int b = 1; b < n_sat; ++b) {
            const int g = order_[j][b];
            xj.segment(kSubStateSize * b, kSubStateSize) = own_mean[g].cast<Real>();
            Pj.block(kSubStateSize * b, 0, kSubStateSize, n_state).setZero();
            Pj.block(0, kSubStateSize * b, n_state, kSubStateSize).setZero();
            Pj.block(kSubStateSize * b, kSubStateSize * b, kSubStateSize, kSubStateSize)
                = own_cov[g];
          }
        }

        apps_[j]->GetFilter()->SetState(xj);
        apps_[j]->GetFilter()->SetCovariance(Pj);
      }
      next_exchange_s_ += cfg.consider_exchange_interval_s;
    }

    // 6. Record truth, estimates, crosslink-geometry + station diagnostics.
    for (int j = 0; j < n_sat; ++j)
      results_.truth_states[j].row(k) = x_truth_[j].cast<double>().transpose();
    for (int i = 0; i < n_links; ++i) {
      results_.range_true_m(k, i) = pr_true(0, i + 1);
      results_.range_rate_true_mps(k, i) = prr_true(0, i + 1);
      results_.range_obs_m(k, i) = pr_obs(0, i + 1);
      results_.range_rate_obs_mps(k, i) = prr_obs(0, i + 1);
      results_.cn0_dbhz(k, i) = ComputeLinkBudgetCn0Dbhz(cfg.link_budget, pr_true(0, i + 1));
      if (tt_on) {
        results_.time_transfer_true_m(k, i)
            = C * (x_truth_[0].cast<double>()(6) - x_truth_[i + 1].cast<double>()(6));
        results_.time_transfer_obs_m(k, i) = tt_obs(0, i + 1);
      }
    }
    RecordFilters(k);
    RecordCentral(k);
  }

  REGISTER_FACTORY_CLASS(Application, IslOdtsCoordinatorApp)

}  // namespace lupnt

#include "lupnt/applications/rover/surface_rover_nav_app.h"

#include <cmath>

#include "lupnt/agents/agent.h"
#include "lupnt/core/constants.h"
#include "lupnt/core/error.h"
#include "lupnt/lupnt.h"
#include "lupnt/simulations/simulation.h"
#include "lupnt/simulations/world.h"

namespace lupnt {

  namespace {
    constexpr int I_DR = 0;    // position error
    constexpr int I_DV = 3;    // velocity error
    constexpr int I_DTH = 6;   // attitude error (nav frame)
    constexpr int I_DBA = 9;   // accel bias error
    constexpr int I_DBG = 12;  // gyro bias error
    constexpr int I_DCB = 15;  // clock-bias error
    constexpr int I_DCD = 16;  // clock-drift error
    constexpr int NX = kSurfaceNavErrorStateSize;

    Mat3d Skew3d(const Vec3d& v) {
      Mat3d S;
      S << 0.0, -v(2), v(1), v(2), 0.0, -v(0), -v(1), v(0), 0.0;
      return S;
    }

    // Exponential map SO3: rotation matrix of the rotation-vector `phi` (Rodrigues).
    Mat3d ExpSO3(const Vec3d& phi) {
      double th = phi.norm();
      Mat3d W = Skew3d(phi);
      if (th < 1e-9) return Mat3d::Identity() + W;
      double s = std::sin(th) / th;
      double c = (1.0 - std::cos(th)) / (th * th);
      return Mat3d::Identity() + s * W + c * (W * W);
    }

    // Rotation-vector (log) of a small rotation matrix, in the same convention as ExpSO3.
    Vec3d LogSO3(const Mat3d& R) {
      double c = std::clamp(0.5 * (R.trace() - 1.0), -1.0, 1.0);
      double th = std::acos(c);
      Vec3d w(R(2, 1) - R(1, 2), R(0, 2) - R(2, 0), R(1, 0) - R(0, 1));
      if (th < 1e-9) return 0.5 * w;
      return (th / (2.0 * std::sin(th))) * w;
    }
  }  // namespace

  namespace {
    Vec3d ParseVec3(const Config& node) {
      return Vec3d(node[0].as<double>(), node[1].as<double>(), node[2].as<double>());
    }
  }  // namespace

  SurfaceRoverNavApp::SurfaceRoverNavApp(Config& config) : Application(config) {
    self_driving_ = true;
    SurfaceNavConfig c;
    c.seed = config["seed"].as<int>(c.seed);
    c.start_epoch_utc = config["start_epoch_utc"].as<std::string>(c.start_epoch_utc);
    c.duration_s = config["duration_s"].as<double>(c.duration_s);
    c.dt_s = config["dt_s"].as<double>(c.dt_s);
    c.dem_max_res_m = config["dem_max_res_m"].as<double>(c.dem_max_res_m);
    c.rover_start_east_m = config["rover_start_east_m"].as<double>(c.rover_start_east_m);
    c.rover_start_north_m = config["rover_start_north_m"].as<double>(c.rover_start_north_m);
    c.rover_speed_mps = config["rover_speed_mps"].as<double>(c.rover_speed_mps);
    c.rover_heading_deg = config["rover_heading_deg"].as<double>(c.rover_heading_deg);
    c.rover_turn_rate_dps = config["rover_turn_rate_dps"].as<double>(c.rover_turn_rate_dps);
    c.rover_clock_bias_s = config["rover_clock_bias_s"].as<double>(c.rover_clock_bias_s);
    c.rover_clock_drift_sps = config["rover_clock_drift_sps"].as<double>(c.rover_clock_drift_sps);
    c.accel_noise_density = config["accel_noise_density"].as<double>(c.accel_noise_density);
    c.accel_bias_rw = config["accel_bias_rw"].as<double>(c.accel_bias_rw);
    c.gyro_noise_density = config["gyro_noise_density"].as<double>(c.gyro_noise_density);
    c.gyro_bias_rw = config["gyro_bias_rw"].as<double>(c.gyro_bias_rw);
    c.accel_bias0 = config["accel_bias0"].as<double>(c.accel_bias0);
    c.gyro_bias0 = config["gyro_bias0"].as<double>(c.gyro_bias0);
    c.elevation_mask_deg = config["elevation_mask_deg"].as<double>(c.elevation_mask_deg);
    c.pseudorange_sigma_m = config["pseudorange_sigma_m"].as<double>(c.pseudorange_sigma_m);
    c.sise_m = config["sise_m"].as<double>(c.sise_m);
    c.init_pos_sigma_m = config["init_pos_sigma_m"].as<double>(c.init_pos_sigma_m);
    c.init_vel_sigma_mps = config["init_vel_sigma_mps"].as<double>(c.init_vel_sigma_mps);
    c.init_att_sigma_deg = config["init_att_sigma_deg"].as<double>(c.init_att_sigma_deg);
    c.init_accel_bias_sigma = config["init_accel_bias_sigma"].as<double>(c.init_accel_bias_sigma);
    c.init_gyro_bias_sigma = config["init_gyro_bias_sigma"].as<double>(c.init_gyro_bias_sigma);
    c.init_clock_bias_sigma_s
        = config["init_clock_bias_sigma_s"].as<double>(c.init_clock_bias_sigma_s);
    c.init_clock_drift_sigma_sps
        = config["init_clock_drift_sigma_sps"].as<double>(c.init_clock_drift_sigma_sps);
    c.dem_sigma_m = config["dem_sigma_m"].as<double>(c.dem_sigma_m);
    c.filter_bias_rw_scale = config["filter_bias_rw_scale"].as<double>(c.filter_bias_rw_scale);
    c.enable_dem_constraint = config["enable_dem_constraint"].as<bool>(c.enable_dem_constraint);

    LUPNT_CHECK(config["satellites"], "SurfaceRoverNavApp requires a `satellites` list",
                "SurfaceRoverNavApp");
    for (const auto& item : config["satellites"]) {
      Config s(item);
      LcrnsSatConfig sat;
      sat.name = s["name"].as<std::string>(std::string("SV"));
      sat.r0_m = ParseVec3(Config(s["r0_m"]));
      sat.v0_mps = ParseVec3(Config(s["v0_mps"]));
      c.satellites.push_back(sat);
    }
    cfg_ = c;
  }

  void SurfaceRoverNavApp::Configure(double t0, const Vec3d& r0, const Vec3d& v0, const Mat3d& R0,
                                     const Vec3d& ba0, const Vec3d& bg0, double cb0, double cd0,
                                     const MatXd& P0) {
    LUPNT_CHECK(P0.rows() == NX && P0.cols() == NX, "P0 must be kSurfaceNavErrorStateSize square",
                "SurfaceRoverNavApp::Configure");
    t_ = t0;
    r_ = r0;
    v_ = v0;
    R_b2n_ = R0;
    ba_ = ba0;
    bg_ = bg0;
    cb_ = cb0;
    cd_ = cd0;
    P_ = P0;
  }

  MatXd SurfaceRoverNavApp::ProcessNoise(double dt) const {
    MatXd Q = MatXd::Zero(NX, NX);
    // Velocity <- accel white noise; attitude <- gyro white noise (Kalibr densities).
    double qv = params_.accel_noise_density * params_.accel_noise_density * dt;
    double qth = params_.gyro_noise_density * params_.gyro_noise_density * dt;
    // Accel/gyro bias random walks.
    double qba = params_.accel_bias_rw * params_.accel_bias_rw * dt;
    double qbg = params_.gyro_bias_rw * params_.gyro_bias_rw * dt;
    for (int i = 0; i < 3; ++i) {
      Q(I_DV + i, I_DV + i) = qv;
      Q(I_DTH + i, I_DTH + i) = qth;
      Q(I_DBA + i, I_DBA + i) = qba;
      Q(I_DBG + i, I_DBG + i) = qbg;
    }
    // Clock two-state random walk.
    double sb2 = params_.clock_bias_process_sigma * params_.clock_bias_process_sigma;
    double sd2 = params_.clock_drift_process_sigma * params_.clock_drift_process_sigma;
    Q(I_DCB, I_DCB) = sb2 * dt + sd2 * dt * dt * dt / 3.0;
    Q(I_DCB, I_DCD) = sd2 * dt * dt / 2.0;
    Q(I_DCD, I_DCB) = sd2 * dt * dt / 2.0;
    Q(I_DCD, I_DCD) = sd2 * dt;
    return Q;
  }

  void SurfaceRoverNavApp::Predict(const SurfaceImuMeasurement& imu, double dt) {
    LUPNT_CHECK(P_.rows() == NX, "filter not configured", "SurfaceRoverNavApp::Predict");

    // Bias-corrected IMU.
    Vec3d w_corr = imu.gyro - bg_;   // body angular rate
    Vec3d a_corr = imu.accel - ba_;  // body specific force
    Vec3d f_n = R_b2n_ * a_corr;     // specific force in the nav frame (gravity excluded)

    // Nominal mechanization: attitude, velocity, position, clock.
    Mat3d R_new = R_b2n_ * ExpSO3(w_corr * dt);
    World* world = agent_ ? agent_->GetWorld() : nullptr;
    double rn = r_.norm();
    Vec3d g = world ? world->Gravity(r_) : Vec3d(-GM_MOON / (rn * rn * rn) * r_);
    Vec3d a_nav = f_n + g;
    r_ = r_ + v_ * dt + 0.5 * a_nav * dt * dt;
    v_ = v_ + a_nav * dt;
    R_b2n_ = R_new;
    cb_ = cb_ + cd_ * dt;

    // Error-state transition (continuous INS model, nav frame ~ inertial over one step).
    MatXd F = MatXd::Zero(NX, NX);
    for (int i = 0; i < 3; ++i) {
      F(I_DR + i, I_DV + i) = 1.0;  // dr' = dv
    }
    F.block<3, 3>(I_DV, I_DTH) = -Skew3d(f_n);  // dv' <- attitude error
    F.block<3, 3>(I_DV, I_DBA) = -R_b2n_;       // dv' <- accel bias
    F.block<3, 3>(I_DTH, I_DBG) = -R_b2n_;      // dtheta' <- gyro bias
    F(I_DCB, I_DCD) = 1.0;                      // clock

    MatXd Phi = MatXd::Identity(NX, NX) + F * dt;
    P_ = Phi * P_ * Phi.transpose() + ProcessNoise(dt);
    P_ = 0.5 * (P_ + P_.transpose());
    t_ += dt;
  }

  void SurfaceRoverNavApp::InjectErrorState(const VecXd& dx) {
    r_ += dx.segment<3>(I_DR);
    v_ += dx.segment<3>(I_DV);
    // Attitude: nav-frame small-angle correction R_true = Exp(dtheta) R_est.
    R_b2n_ = ExpSO3(dx.segment<3>(I_DTH)) * R_b2n_;
    ba_ += dx.segment<3>(I_DBA);
    bg_ += dx.segment<3>(I_DBG);
    cb_ += dx(I_DCB);
    cd_ += dx(I_DCD);
  }

  void SurfaceRoverNavApp::UpdateScalar(const VecXd& H, double z_pred, double z_meas,
                                        double variance) {
    LUPNT_CHECK(H.size() == NX, "H must be size kSurfaceNavErrorStateSize",
                "SurfaceRoverNavApp::UpdateScalar");
    VecXd PHt = P_ * H;
    double S = H.dot(PHt) + variance;
    VecXd K = PHt / S;
    double dz = z_meas - z_pred;
    InjectErrorState(K * dz);
    // Joseph-form covariance update.
    MatXd I = MatXd::Identity(NX, NX);
    MatXd IKH = I - K * H.transpose();
    P_ = IKH * P_ * IKH.transpose() + (K * K.transpose()) * variance;
    P_ = 0.5 * (P_ + P_.transpose());
  }

  void SurfaceRoverNavApp::UpdateLans(const std::vector<SurfaceLansMeasurement>& meas) {
    for (const auto& m : meas) {
      // Rebuild the nominal context each iteration: the previous scalar update injected its
      // correction into r_/cb_, and the model linearizes about the *current* nominal state.
      NavErrorContext nom;
      nom.r = r_;
      nom.clock_bias_s = cb_;
      nom.error_state_size = NX;
      MatXd H;
      MeasData md = m.Compute(nom, &H);
      UpdateScalar(H.row(0).transpose(), md.value(0), m.pseudorange_m, md.covariance(0, 0));
    }
  }

  // --- Self-driving scenario (agent-based Simulation) --------------------------

  void SurfaceRoverNavApp::Setup() {
    if (!self_driving_) return;  // params-struct driver runs the filter externally
    LUPNT_CHECK(agent_, "Agent not set", "SurfaceRoverNavApp::Setup");
    // Defer the heavy precompute to the first Step so any config set programmatically between
    // Simulation construction and run() is honored; here we only schedule the per-epoch Steps.
    dt_ = cfg_.dt_s;
    Simulation* sim = agent_->GetSimulation();
    sim->Schedule(0.0, [this](Real t) { Step(t); }, 1.0 / dt_, Event::Priority::APPLICATION);
  }

  void SurfaceRoverNavApp::InitScenario() {
    World* world = agent_->GetWorld();
    LUPNT_CHECK(world && world->HasTerrain(),
                "SurfaceRoverNavApp requires a World with a `dem:` terrain block",
                "SurfaceRoverNavApp::Setup");

    rng_.seed(cfg_.seed);
    nd_.reset();

    N_ = std::max(2, static_cast<int>(std::round(cfg_.duration_s / cfg_.dt_s)) + 1);
    dt_ = cfg_.dt_s;
    n_sat_ = static_cast<int>(cfg_.satellites.size());
    const double dt = dt_;
    const double sqrt_dt = std::sqrt(dt);

    // ---- 1. Local ENU terrain frame (from the shared World) -------------------
    const LunarDem& dem = world->GetDem();
    R_enu2pa_ = world->REnuToWorld();
    R_pa2enu_ = R_enu2pa_.transpose();
    r_center_pa_ = world->SiteCenterWorld();
    up_hat_pa_ = R_enu2pa_.col(2);
    ds_ = std::max(cfg_.dem_max_res_m, 1.0);

    auto EnuToPa
        = [&](double E, double Nn, double U) -> Vec3d { return world->EnuToWorld(E, Nn, U); };
    auto Elevation = [&](double E, double Nn) -> double { return world->GetElevation(E, Nn); };

    // Body(vehicle)-to-nav rotation for a rover heading `hdg` [rad]: x=forward, y=left, z=up.
    auto BodyToPa = [&](double hdg) -> Mat3d {
      Vec3d xb(std::cos(hdg), std::sin(hdg), 0.0);
      Vec3d zb(0.0, 0.0, 1.0);
      Vec3d yb = zb.cross(xb);
      Mat3d R_body2enu;
      R_body2enu.col(0) = xb;
      R_body2enu.col(1) = yb;
      R_body2enu.col(2) = zb;
      return R_enu2pa_ * R_body2enu;
    };

    // ---- 2. Rover truth: ENU path, Moon-fixed position/velocity/accel, attitude/rate ----
    Ee_.assign(N_, 0.0);
    Nn_.assign(N_, 0.0);
    Uu_.assign(N_, 0.0);
    std::vector<double> Hdg(N_, 0.0);
    r_truth_.assign(N_, Vec3d::Zero());
    v_truth_.assign(N_, Vec3d::Zero());
    R_truth_.assign(N_, Mat3d::Identity());
    for (int k = 0; k < N_; ++k) {
      double t = k * dt;
      Hdg[k] = (cfg_.rover_heading_deg + cfg_.rover_turn_rate_dps * t) * RAD;
      if (k == 0) {
        Ee_[0] = cfg_.rover_start_east_m;
        Nn_[0] = cfg_.rover_start_north_m;
      } else {
        Ee_[k] = Ee_[k - 1] + cfg_.rover_speed_mps * std::cos(Hdg[k - 1]) * dt;
        Nn_[k] = Nn_[k - 1] + cfg_.rover_speed_mps * std::sin(Hdg[k - 1]) * dt;
      }
      Uu_[k] = Elevation(Ee_[k], Nn_[k]);
      R_truth_[k] = BodyToPa(Hdg[k]);
    }
    for (int k = 0; k < N_; ++k) r_truth_[k] = EnuToPa(Ee_[k], Nn_[k], Uu_[k]);
    for (int k = 0; k < N_; ++k) {
      int kp = std::min(k + 1, N_ - 1), km = std::max(k - 1, 0);
      double span = (kp - km) * dt;
      v_truth_[k]
          = (span > 0.0) ? Vec3d((r_truth_[kp] - r_truth_[km]) / span) : Vec3d(Vec3d::Zero());
    }

    f_body_truth_.assign(N_, Vec3d::Zero());
    w_body_truth_.assign(N_, Vec3d::Zero());
    for (int k = 0; k < N_; ++k) {
      int kp = std::min(k + 1, N_ - 1), km = std::max(k - 1, 0);
      double span = (kp - km) * dt;
      Vec3d a_total
          = (span > 0.0) ? Vec3d((v_truth_[kp] - v_truth_[km]) / span) : Vec3d(Vec3d::Zero());
      f_body_truth_[k] = R_truth_[k].transpose() * (a_total - world->Gravity(r_truth_[k]));
      Mat3d Rdot
          = (span > 0.0) ? Mat3d((R_truth_[kp] - R_truth_[km]) / span) : Mat3d(Mat3d::Zero());
      Mat3d Wx = R_truth_[k].transpose() * Rdot;  // [w x] in body
      w_body_truth_[k] = Vec3d(Wx(2, 1), Wx(0, 2), Wx(1, 0));
    }

    // Truth IMU biases: constant offset + random walk (Kalibr).
    Vec3d ba_truth = Gauss3(cfg_.accel_bias0);
    Vec3d bg_truth = Gauss3(cfg_.gyro_bias0);
    ba_truth_k_.assign(N_, Vec3d::Zero());
    bg_truth_k_.assign(N_, Vec3d::Zero());
    for (int k = 0; k < N_; ++k) {
      if (k > 0) {
        ba_truth += Gauss3(cfg_.accel_bias_rw * sqrt_dt);
        bg_truth += Gauss3(cfg_.gyro_bias_rw * sqrt_dt);
      }
      ba_truth_k_[k] = ba_truth;
      bg_truth_k_[k] = bg_truth;
    }

    // ---- 3. Propagate LCRNS truth orbits (Keplerian) --------------------------
    const Real t0_tdb = ConvertTime(GregorianToTime(cfg_.start_epoch_utc), Time::UTC, Time::TDB);
    CartesianTwoBodyDynamics sat_dyn(GM_MOON);
    sat_dyn.SetTimeStep(dt);
    std::vector<Vec6> sat_ci(n_sat_);
    for (int j = 0; j < n_sat_; ++j)
      sat_ci[j]
          = (Vec6() << cfg_.satellites[j].r0_m.cast<Real>(), cfg_.satellites[j].v0_mps.cast<Real>())
                .finished();
    sat_pa_.assign(N_, std::vector<Vec3d>(n_sat_));
    for (int k = 0; k < N_; ++k) {
      Real tk = t0_tdb + k * dt;
      if (k > 0) {
        Real tkm1 = t0_tdb + (k - 1) * dt;
        for (int j = 0; j < n_sat_; ++j) {
          Cart6 st = sat_dyn.Propagate(Cart6(sat_ci[j], Frame::MOON_CI), tkm1, tk, nullptr);
          sat_ci[j] = st.head(6);
        }
      }
      for (int j = 0; j < n_sat_; ++j) {
        Vec6 rv_pa = ConvertFrame(tk, sat_ci[j], Frame::MOON_CI, Frame::MOON_PA);
        sat_pa_[k][j] = rv_pa.head(3).cast<double>();
      }
    }

    sise_bias_.assign(n_sat_, 0.0);
    for (int j = 0; j < n_sat_; ++j) sise_bias_[j] = Gauss(cfg_.sise_m);

    // ---- 4. Filter tuning + perturbed initial estimate ------------------------
    params_.accel_noise_density = cfg_.accel_noise_density;
    params_.accel_bias_rw = cfg_.filter_bias_rw_scale * cfg_.accel_bias_rw;
    params_.gyro_noise_density = cfg_.gyro_noise_density;
    params_.gyro_bias_rw = cfg_.filter_bias_rw_scale * cfg_.gyro_bias_rw;
    params_.pseudorange_sigma_m = cfg_.pseudorange_sigma_m;
    params_.sise_m = cfg_.sise_m;
    params_.dem_sigma_m = cfg_.dem_sigma_m;

    Vec3d r0 = r_truth_[0] + Gauss3(cfg_.init_pos_sigma_m);
    Vec3d v0 = v_truth_[0] + Gauss3(cfg_.init_vel_sigma_mps);
    Vec3d att_err0 = Gauss3(cfg_.init_att_sigma_deg * RAD);
    // R_true = Exp(dtheta) R_est  =>  R_est = Exp(-dtheta) R_true.
    Mat3d dR0 = Mat3d::Identity() - Skew3d(att_err0);  // small-angle Exp(-att_err0)
    Mat3d R0 = dR0 * R_truth_[0];
    double cb0 = ClockBiasTruth(0) + Gauss(cfg_.init_clock_bias_sigma_s);
    double cd0 = cfg_.rover_clock_drift_sps + Gauss(cfg_.init_clock_drift_sigma_sps);

    VecXd p0(kSurfaceNavErrorStateSize);
    p0 << Vec3d::Constant(cfg_.init_pos_sigma_m * cfg_.init_pos_sigma_m),
        Vec3d::Constant(cfg_.init_vel_sigma_mps * cfg_.init_vel_sigma_mps),
        Vec3d::Constant(std::pow(cfg_.init_att_sigma_deg * RAD, 2)),
        Vec3d::Constant(cfg_.init_accel_bias_sigma * cfg_.init_accel_bias_sigma),
        Vec3d::Constant(cfg_.init_gyro_bias_sigma * cfg_.init_gyro_bias_sigma),
        cfg_.init_clock_bias_sigma_s * cfg_.init_clock_bias_sigma_s,
        cfg_.init_clock_drift_sigma_sps * cfg_.init_clock_drift_sigma_sps;
    MatXd P0 = p0.asDiagonal();
    Configure(0.0, r0, v0, R0, Vec3d::Zero(), Vec3d::Zero(), cb0, cd0, P0);

    // ---- 5. Results allocation ------------------------------------------------
    res_ = SurfaceNavResults();
    res_.site_id = dem.site().id;
    res_.site_name = dem.site().name;
    res_.site_lat_deg = dem.site().lat_deg;
    res_.site_lon_deg = dem.site().lon_deg;
    res_.dem_x = dem.x();
    res_.dem_y = dem.y();
    res_.dem_elevation = dem.elevation();
    res_.dem_center_x = dem.center_x();
    res_.dem_center_y = dem.center_y();
    res_.time_s = VecXd::Zero(N_);
    res_.pos_err_enu = MatXd::Zero(N_, 3);
    res_.pos_sigma_enu = MatXd::Zero(N_, 3);
    res_.pos_err_norm = VecXd::Zero(N_);
    res_.clock_bias_err = VecXd::Zero(N_);
    res_.clock_bias_sigma = VecXd::Zero(N_);
    res_.n_visible = VecXi::Zero(N_);
    res_.accel_bias_err = MatXd::Zero(N_, 3);
    res_.accel_bias_sigma = MatXd::Zero(N_, 3);
    res_.gyro_bias_err = MatXd::Zero(N_, 3);
    res_.gyro_bias_sigma = MatXd::Zero(N_, 3);
    res_.att_err_deg = MatXd::Zero(N_, 3);
    res_.att_sigma_deg = MatXd::Zero(N_, 3);
    res_.rover_track_enu_truth = MatXd::Zero(N_, 2);
    res_.rover_track_enu_est = MatXd::Zero(N_, 2);
    res_.rover_alt_truth = VecXd::Zero(N_);
    for (int j = 0; j < n_sat_; ++j) res_.satellite_names.push_back(cfg_.satellites[j].name);

    // Epoch 0: seed the host agent's truth state and log the initial estimate/covariance.
    auto* rover = dynamic_cast<AgentWithDynamics*>(agent_);
    if (rover) {
      Vec6 rv;
      rv << r_truth_[0].cast<Real>(), v_truth_[0].cast<Real>();
      rover->SetTime(0.0);
      rover->SetState(Cart6(rv, Frame::MOON_PA));
    }
    LogEpoch(0);
  }

  void SurfaceRoverNavApp::Step(Real t) {
    if (!self_driving_) return;
    if (!initialized_) {
      InitScenario();  // precompute truth/measurements/initial estimate on first Step
      initialized_ = true;
    }
    int k = static_cast<int>(std::lround(t.val() / dt_));
    if (k < 1 || k >= N_) return;
    const double dt = dt_;
    const double sqrt_dt = std::sqrt(dt);
    World* world = agent_->GetWorld();

    // Advance the rover agent's truth state (kinematic surface trajectory).
    auto* rover = dynamic_cast<AgentWithDynamics*>(agent_);
    if (rover) {
      Vec6 rv;
      rv << r_truth_[k].cast<Real>(), v_truth_[k].cast<Real>();
      rover->SetTime(k * dt);
      rover->SetState(Cart6(rv, Frame::MOON_PA));
    }

    // Simulate the IMU (Kalibr: bias + white noise) over [k-1, k] and predict.
    SurfaceImuMeasurement imu;
    imu.accel
        = f_body_truth_[k - 1] + ba_truth_k_[k - 1] + Gauss3(cfg_.accel_noise_density / sqrt_dt);
    imu.gyro
        = w_body_truth_[k - 1] + bg_truth_k_[k - 1] + Gauss3(cfg_.gyro_noise_density / sqrt_dt);
    Predict(imu, dt);

    // LANS pseudorange updates for visible satellites.
    std::vector<SurfaceLansMeasurement> lans;
    for (int j = 0; j < n_sat_; ++j) {
      Vec3d los = sat_pa_[k][j] - r_truth_[k];
      double range = los.norm();
      double sin_el = (range > 0.0) ? (los.dot(up_hat_pa_) / range) : -1.0;
      if (sin_el < std::sin(cfg_.elevation_mask_deg * RAD)) continue;
      SurfaceLansMeasurement m;
      m.r_sat = sat_pa_[k][j];
      m.sigma_m = cfg_.pseudorange_sigma_m;
      m.sise_m = cfg_.sise_m;
      m.pseudorange_m
          = range + C * ClockBiasTruth(k) + sise_bias_[j] + Gauss(cfg_.pseudorange_sigma_m);
      lans.push_back(m);
    }
    res_.n_visible(k) = static_cast<int>(lans.size());
    UpdateLans(lans);

    // DEM altitude constraint: enforce U(estimate) == DEM(East, North).
    if (cfg_.enable_dem_constraint) {
      Vec3d enu = R_pa2enu_ * (r_ - r_center_pa_);
      double E = enu(0), Ne = enu(1), U = enu(2);
      double dem_elev = world->GetElevation(E, Ne);
      double sE = (world->GetElevation(E + ds_, Ne) - world->GetElevation(E - ds_, Ne)) / (2 * ds_);
      double sN = (world->GetElevation(E, Ne + ds_) - world->GetElevation(E, Ne - ds_)) / (2 * ds_);
      VecXd H = VecXd::Zero(kSurfaceNavErrorStateSize);
      H.head(3) = (R_pa2enu_.row(2) - sE * R_pa2enu_.row(0) - sN * R_pa2enu_.row(1)).transpose();
      UpdateScalar(H, U - dem_elev, 0.0, cfg_.dem_sigma_m * cfg_.dem_sigma_m);
    }

    LogEpoch(k);
  }

  void SurfaceRoverNavApp::LogEpoch(int k) {
    const MatXd& P = P_;
    Vec3d r_est = r_;
    Vec3d err_pa = r_truth_[k] - r_est;
    Vec3d err_enu = R_pa2enu_ * err_pa;
    res_.time_s(k) = k * dt_;
    res_.pos_err_enu.row(k) = err_enu.transpose();
    res_.pos_err_norm(k) = err_pa.norm();
    Mat3d P_enu = R_pa2enu_ * P.block<3, 3>(0, 0) * R_enu2pa_;
    for (int i = 0; i < 3; ++i) res_.pos_sigma_enu(k, i) = std::sqrt(std::max(0.0, P_enu(i, i)));
    res_.clock_bias_err(k) = ClockBiasTruth(k) - cb_;
    res_.clock_bias_sigma(k) = std::sqrt(std::max(0.0, P(15, 15)));

    Vec3d ba_err = ba_truth_k_[k] - ba_;
    Vec3d bg_err = bg_truth_k_[k] - bg_;
    res_.accel_bias_err.row(k) = ba_err.transpose();
    res_.gyro_bias_err.row(k) = bg_err.transpose();
    Vec3d att_err = LogSO3(R_truth_[k] * R_b2n_.transpose());  // nav-frame error
    res_.att_err_deg.row(k) = (att_err / RAD).transpose();
    for (int i = 0; i < 3; ++i) {
      res_.accel_bias_sigma(k, i) = std::sqrt(std::max(0.0, P(9 + i, 9 + i)));
      res_.gyro_bias_sigma(k, i) = std::sqrt(std::max(0.0, P(12 + i, 12 + i)));
      res_.att_sigma_deg(k, i) = std::sqrt(std::max(0.0, P(6 + i, 6 + i))) / RAD;
    }

    Vec3d enu_est = R_pa2enu_ * (r_est - r_center_pa_);
    res_.rover_track_enu_est.row(k) = Vec2d(enu_est(0), enu_est(1)).transpose();
    res_.rover_track_enu_truth.row(k) = Vec2d(Ee_[k], Nn_[k]).transpose();
    res_.rover_alt_truth(k) = Uu_[k];
  }

  void SurfaceRoverNavApp::Log(Real /*t*/) {}

  REGISTER_FACTORY_CLASS(Application, SurfaceRoverNavApp)

}  // namespace lupnt

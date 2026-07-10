#include "lupnt/applications/lander/lander_nav_app.h"

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
    constexpr int NX = kLanderNavErrorStateSize;

    Mat3d Skew3d(const Vec3d& v) {
      Mat3d S;
      S << 0.0, -v(2), v(1), v(2), 0.0, -v(0), -v(1), v(0), 0.0;
      return S;
    }

    // Unit quaternion of the rotation-vector `phi` (SO3 exponential map on the quaternion
    // manifold): the small-rotation "delta quaternion" dq(phi) = [cos|phi|/2,
    // sin(|phi|/2) phi/|phi|], used to compose gyro increments and inject the MEKF
    // multiplicative attitude error.
    Eigen::Quaterniond QuatFromRotVec(const Vec3d& phi) {
      double th = phi.norm();
      if (th < 1e-9) {
        return Eigen::Quaterniond(1.0, 0.5 * phi.x(), 0.5 * phi.y(), 0.5 * phi.z()).normalized();
      }
      double half = 0.5 * th;
      double s = std::sin(half) / th;
      return Eigen::Quaterniond(std::cos(half), s * phi.x(), s * phi.y(), s * phi.z());
    }

    Vec3d Gravity(const Vec3d& r) {
      double rn = r.norm();
      return (rn > 0.0) ? Vec3d(-GM_MOON / (rn * rn * rn) * r) : Vec3d(Vec3d::Zero());
    }

    // Rotation-vector (log) of a small rotation matrix (nav-frame attitude error metric).
    Vec3d LogSO3(const Mat3d& R) {
      double c = std::clamp(0.5 * (R.trace() - 1.0), -1.0, 1.0);
      double th = std::acos(c);
      Vec3d w(R(2, 1) - R(1, 2), R(0, 2) - R(2, 0), R(1, 0) - R(0, 1));
      if (th < 1e-9) return 0.5 * w;
      return (th / (2.0 * std::sin(th))) * w;
    }

    // Smoothstep S(tau) = tau^2 (3 - 2 tau), clamped to [0, 1].
    double SmoothStep(double tau) {
      tau = std::clamp(tau, 0.0, 1.0);
      return tau * tau * (3.0 - 2.0 * tau);
    }

    Vec3d ParseVec3(const Config& node) {
      return Vec3d(node[0].as<double>(), node[1].as<double>(), node[2].as<double>());
    }
  }  // namespace

  void LanderNavApp::Configure(double t0, const Vec3d& r0, const Vec3d& v0, const Mat3d& R0,
                               const Vec3d& ba0, const Vec3d& bg0, double cb0, double cd0,
                               const MatXd& P0) {
    LUPNT_CHECK(P0.rows() == NX && P0.cols() == NX, "P0 must be kLanderNavErrorStateSize square",
                "LanderNavApp::Configure");
    t_ = t0;
    r_ = r0;
    v_ = v0;
    q_b2n_ = Eigen::Quaterniond(R0).normalized();  // seed the MEKF quaternion from the DCM
    ba_ = ba0;
    bg_ = bg0;
    cb_ = cb0;
    cd_ = cd0;
    P_ = P0;
  }

  MatXd LanderNavApp::ProcessNoise(double dt) const {
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

  void LanderNavApp::Predict(const SurfaceImuMeasurement& imu, double dt) {
    LUPNT_CHECK(P_.rows() == NX, "filter not configured", "LanderNavApp::Predict");

    // Bias-corrected IMU.
    Vec3d w_corr = imu.gyro - bg_;            // body angular rate
    Vec3d a_corr = imu.accel - ba_;           // body specific force
    Mat3d R_old = q_b2n_.toRotationMatrix();  // body-to-nav DCM from the nominal quaternion
    Vec3d f_n = R_old * a_corr;               // specific force in the nav frame (gravity excluded)

    // Nominal mechanization: attitude quaternion (gyro increment), velocity, position, clock.
    // q_b2n <- q_b2n (x) dq(w_corr dt): the body-frame rotation increment right-multiplies.
    q_b2n_ = (q_b2n_ * QuatFromRotVec(w_corr * dt)).normalized();
    Mat3d R_new = q_b2n_.toRotationMatrix();
    World* world = agent_ ? agent_->GetWorld() : nullptr;
    Vec3d a_nav = f_n + (world ? world->Gravity(r_) : Gravity(r_));
    r_ = r_ + v_ * dt + 0.5 * a_nav * dt * dt;
    v_ = v_ + a_nav * dt;
    cb_ = cb_ + cd_ * dt;

    // Error-state transition (continuous INS model, nav frame ~ inertial over one step).
    MatXd F = MatXd::Zero(NX, NX);
    for (int i = 0; i < 3; ++i) {
      F(I_DR + i, I_DV + i) = 1.0;  // dr' = dv
    }
    F.block<3, 3>(I_DV, I_DTH) = -Skew3d(f_n);  // dv' <- attitude error
    F.block<3, 3>(I_DV, I_DBA) = -R_new;        // dv' <- accel bias
    F.block<3, 3>(I_DTH, I_DBG) = -R_new;       // dtheta' <- gyro bias
    F(I_DCB, I_DCD) = 1.0;                      // clock

    MatXd Phi = MatXd::Identity(NX, NX) + F * dt;
    P_ = Phi * P_ * Phi.transpose() + ProcessNoise(dt);
    P_ = 0.5 * (P_ + P_.transpose());
    t_ += dt;
  }

  void LanderNavApp::InjectErrorState(const VecXd& dx) {
    r_ += dx.segment<3>(I_DR);
    v_ += dx.segment<3>(I_DV);
    // MEKF multiplicative attitude reset: inject the nav-frame error rotation into the nominal
    // quaternion, q_b2n <- dq(dtheta) (x) q_b2n (equiv. R_true = Exp(dtheta) R_est), then the
    // 3-parameter error is implicitly zeroed (absorbed by the quaternion).
    q_b2n_ = (QuatFromRotVec(dx.segment<3>(I_DTH)) * q_b2n_).normalized();
    ba_ += dx.segment<3>(I_DBA);
    bg_ += dx.segment<3>(I_DBG);
    cb_ += dx(I_DCB);
    cd_ += dx(I_DCD);
  }

  void LanderNavApp::UpdateScalar(const VecXd& H, double z_pred, double z_meas, double variance) {
    LUPNT_CHECK(H.size() == NX, "H must be size kLanderNavErrorStateSize",
                "LanderNavApp::UpdateScalar");
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

  void LanderNavApp::UpdateVector(const MatXd& H, const VecXd& y, const MatXd& R) {
    LUPNT_CHECK(H.cols() == NX, "H must have kLanderNavErrorStateSize columns",
                "LanderNavApp::UpdateVector");
    MatXd PHt = P_ * H.transpose();  // NX x m
    MatXd S = H * PHt + R;           // m x m
    MatXd K = PHt * S.inverse();     // NX x m
    InjectErrorState(K * y);
    // Joseph-form covariance update.
    MatXd I = MatXd::Identity(NX, NX);
    MatXd IKH = I - K * H;
    P_ = IKH * P_ * IKH.transpose() + K * R * K.transpose();
    P_ = 0.5 * (P_ + P_.transpose());
  }

  void LanderNavApp::UpdateLans(const std::vector<SurfaceLansMeasurement>& meas) {
    for (const auto& m : meas) {
      // Rebuild the nominal context each iteration (prior update mutated r_/cb_).
      NavErrorContext nom;
      nom.r = r_;
      nom.clock_bias_s = cb_;
      nom.error_state_size = NX;
      MatXd H;
      MeasData md = m.Compute(nom, &H);
      UpdateScalar(H.row(0).transpose(), md.value(0), m.pseudorange_m, md.covariance(0, 0));
    }
  }

  void LanderNavApp::UpdateAltimeter(const LanderAltimeterMeasurement& meas) {
    // `meas.predicted_altitude_m` / `meas.h_pos` are the DEM-derived nominal prediction and
    // position partial, supplied by the driving simulation.
    NavErrorContext nom;
    nom.error_state_size = NX;
    MatXd H;
    MeasData md = meas.Compute(nom, &H);
    UpdateScalar(H.row(0).transpose(), md.value(0), meas.altitude_m, md.covariance(0, 0));
  }

  void LanderNavApp::UpdateCrater(const LanderCraterMeasurement& meas) {
    NavErrorContext nom;
    nom.r = r_;
    nom.R_b2n = q_b2n_.toRotationMatrix();
    nom.error_state_size = NX;
    MatXd H;
    MeasData md = meas.Compute(nom, &H);
    if (md.value.size() == 0) return;           // degenerate geometry
    VecXd y = meas.los_body - Vec3d(md.value);  // residual (unit-vector, ~ tangent plane)
    UpdateVector(H, y, md.covariance);
  }

  void LanderNavApp::UpdateCraters(const std::vector<LanderCraterMeasurement>& meas) {
    for (const auto& m : meas) UpdateCrater(m);
  }

  // --- Config constructor + self-driving scenario -----------------------------

  LanderNavApp::LanderNavApp(Config& config) : Application(config) {
    self_driving_ = true;
    LanderNavConfig c;
    c.seed = config["seed"].as<int>(c.seed);
    c.start_epoch_utc = config["start_epoch_utc"].as<std::string>(c.start_epoch_utc);
    c.duration_s = config["duration_s"].as<double>(c.duration_s);
    c.dt_s = config["dt_s"].as<double>(c.dt_s);
    c.dem_max_res_m = config["dem_max_res_m"].as<double>(c.dem_max_res_m);
    c.descent_start_east_m = config["descent_start_east_m"].as<double>(c.descent_start_east_m);
    c.descent_start_north_m = config["descent_start_north_m"].as<double>(c.descent_start_north_m);
    c.descent_end_east_m = config["descent_end_east_m"].as<double>(c.descent_end_east_m);
    c.descent_end_north_m = config["descent_end_north_m"].as<double>(c.descent_end_north_m);
    c.descent_start_alt_m = config["descent_start_alt_m"].as<double>(c.descent_start_alt_m);
    c.descent_end_alt_m = config["descent_end_alt_m"].as<double>(c.descent_end_alt_m);
    c.descent_heading_deg = config["descent_heading_deg"].as<double>(c.descent_heading_deg);
    c.lander_clock_bias_s = config["lander_clock_bias_s"].as<double>(c.lander_clock_bias_s);
    c.lander_clock_drift_sps
        = config["lander_clock_drift_sps"].as<double>(c.lander_clock_drift_sps);
    c.accel_noise_density = config["accel_noise_density"].as<double>(c.accel_noise_density);
    c.accel_bias_rw = config["accel_bias_rw"].as<double>(c.accel_bias_rw);
    c.gyro_noise_density = config["gyro_noise_density"].as<double>(c.gyro_noise_density);
    c.gyro_bias_rw = config["gyro_bias_rw"].as<double>(c.gyro_bias_rw);
    c.accel_bias0 = config["accel_bias0"].as<double>(c.accel_bias0);
    c.gyro_bias0 = config["gyro_bias0"].as<double>(c.gyro_bias0);
    c.enable_altimeter = config["enable_altimeter"].as<bool>(c.enable_altimeter);
    c.altimeter_sigma_m = config["altimeter_sigma_m"].as<double>(c.altimeter_sigma_m);
    c.altimeter_max_range_m = config["altimeter_max_range_m"].as<double>(c.altimeter_max_range_m);
    c.enable_craters = config["enable_craters"].as<bool>(c.enable_craters);
    c.n_craters = config["n_craters"].as<int>(c.n_craters);
    c.crater_field_radius_m = config["crater_field_radius_m"].as<double>(c.crater_field_radius_m);
    c.camera_fov_deg = config["camera_fov_deg"].as<double>(c.camera_fov_deg);
    c.max_craters_per_epoch = config["max_craters_per_epoch"].as<int>(c.max_craters_per_epoch);
    c.crater_sigma_arcsec = config["crater_sigma_arcsec"].as<double>(c.crater_sigma_arcsec);
    c.enable_lunanet = config["enable_lunanet"].as<bool>(c.enable_lunanet);
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
    c.filter_bias_rw_scale = config["filter_bias_rw_scale"].as<double>(c.filter_bias_rw_scale);
    if (config["satellites"]) {
      for (const auto& item : config["satellites"]) {
        Config s(item);
        LcrnsSatConfig sat;
        sat.name = s["name"].as<std::string>(std::string("SV"));
        sat.r0_m = ParseVec3(Config(s["r0_m"]));
        sat.v0_mps = ParseVec3(Config(s["v0_mps"]));
        c.satellites.push_back(sat);
      }
    }
    cfg_ = c;
  }

  void LanderNavApp::Setup() {
    if (!self_driving_) return;
    LUPNT_CHECK(agent_, "Agent not set", "LanderNavApp::Setup");
    // Defer the heavy precompute to the first Step so a reference trajectory set via
    // SetReferenceTrajectoryEnu (after construction, before run()) is honored.
    dt_ = cfg_.dt_s;
    Simulation* sim = agent_->GetSimulation();
    sim->Schedule(0.0, [this](Real t) { Step(t); }, 1.0 / dt_, Event::Priority::APPLICATION);
  }

  void LanderNavApp::InitScenario() {
    World* world = agent_->GetWorld();
    LUPNT_CHECK(world && world->HasTerrain(),
                "LanderNavApp requires a World with a `dem:` terrain block", "LanderNavApp::Setup");

    rng_.seed(cfg_.seed);
    nd_.reset();
    ud_.reset();

    dt_ = cfg_.dt_s;
    n_sat_ = static_cast<int>(cfg_.satellites.size());
    crater_sigma_rad_ = cfg_.crater_sigma_arcsec * RAD / 3600.0;
    const double dt = dt_;
    const double sqrt_dt = std::sqrt(dt);

    const bool use_ref_traj = cfg_.ref_traj_enu.rows() > 0;
    LUPNT_CHECK(!use_ref_traj || cfg_.ref_traj_enu.cols() == 3,
                "ref_traj_enu must have 3 columns (East, North, Up)", "LanderNavApp::Setup");
    N_ = use_ref_traj ? static_cast<int>(cfg_.ref_traj_enu.rows())
                      : std::max(2, static_cast<int>(std::round(cfg_.duration_s / cfg_.dt_s)) + 1);

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

    // ---- 2. Synthetic crater-landmark map (scattered on the terrain) ----------
    n_crat_ = cfg_.enable_craters ? std::max(0, cfg_.n_craters) : 0;
    crater_pa_.assign(n_crat_, Vec3d::Zero());
    res_.crater_enu = MatXd::Zero(n_crat_, 2);
    for (int j = 0; j < n_crat_; ++j) {
      double ang = 2.0 * M_PI * ud_(rng_);
      double rad = cfg_.crater_field_radius_m * std::sqrt(ud_(rng_));
      double E = rad * std::cos(ang), Nn = rad * std::sin(ang);
      double U = Elevation(E, Nn);
      crater_pa_[j] = EnuToPa(E, Nn, U);
      res_.crater_enu(j, 0) = E;
      res_.crater_enu(j, 1) = Nn;
    }

    // ---- 3. Descent truth: ENU path, Moon-fixed pos/vel/accel, attitude/rate ---
    const double hdg = cfg_.descent_heading_deg * RAD;
    Ee_.assign(N_, 0.0);
    Nn_.assign(N_, 0.0);
    Uu_.assign(N_, 0.0);
    Alt_.assign(N_, 0.0);
    r_truth_.assign(N_, Vec3d::Zero());
    v_truth_.assign(N_, Vec3d::Zero());
    R_truth_.assign(N_, Mat3d::Identity());
    for (int k = 0; k < N_; ++k) {
      if (use_ref_traj) {
        Ee_[k] = cfg_.ref_traj_enu(k, 0);
        Nn_[k] = cfg_.ref_traj_enu(k, 1);
        Uu_[k] = cfg_.ref_traj_enu(k, 2);
        Alt_[k] = Uu_[k] - Elevation(Ee_[k], Nn_[k]);
      } else {
        double tau = (cfg_.duration_s > 0.0) ? (k * dt / cfg_.duration_s) : 1.0;
        double s = SmoothStep(tau);
        Ee_[k]
            = cfg_.descent_start_east_m + (cfg_.descent_end_east_m - cfg_.descent_start_east_m) * s;
        Nn_[k] = cfg_.descent_start_north_m
                 + (cfg_.descent_end_north_m - cfg_.descent_start_north_m) * s;
        Alt_[k]
            = cfg_.descent_start_alt_m + (cfg_.descent_end_alt_m - cfg_.descent_start_alt_m) * s;
        Uu_[k] = Elevation(Ee_[k], Nn_[k]) + Alt_[k];
      }
      R_truth_[k] = BodyToPa(hdg);
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
      Mat3d Wx = R_truth_[k].transpose() * Rdot;
      w_body_truth_[k] = Vec3d(Wx(2, 1), Wx(0, 2), Wx(1, 0));
    }

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

    // ---- 4. Propagate LunaNet truth orbits (Keplerian) ------------------------
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

    // ---- 5. Filter tuning + perturbed initial estimate ------------------------
    params_.accel_noise_density = cfg_.accel_noise_density;
    params_.accel_bias_rw = cfg_.filter_bias_rw_scale * cfg_.accel_bias_rw;
    params_.gyro_noise_density = cfg_.gyro_noise_density;
    params_.gyro_bias_rw = cfg_.filter_bias_rw_scale * cfg_.gyro_bias_rw;
    params_.pseudorange_sigma_m = cfg_.pseudorange_sigma_m;
    params_.sise_m = cfg_.sise_m;
    params_.altimeter_sigma_m = cfg_.altimeter_sigma_m;
    params_.crater_sigma_rad = crater_sigma_rad_;

    Vec3d r0 = r_truth_[0] + Gauss3(cfg_.init_pos_sigma_m);
    Vec3d v0 = v_truth_[0] + Gauss3(cfg_.init_vel_sigma_mps);
    Vec3d att_err0 = Gauss3(cfg_.init_att_sigma_deg * RAD);
    Mat3d dR0 = Mat3d::Identity() - Skew3d(att_err0);
    Mat3d R0 = dR0 * R_truth_[0];
    double cb0 = ClockBiasTruth(0) + Gauss(cfg_.init_clock_bias_sigma_s);
    double cd0 = cfg_.lander_clock_drift_sps + Gauss(cfg_.init_clock_drift_sigma_sps);

    VecXd p0(kLanderNavErrorStateSize);
    p0 << Vec3d::Constant(cfg_.init_pos_sigma_m * cfg_.init_pos_sigma_m),
        Vec3d::Constant(cfg_.init_vel_sigma_mps * cfg_.init_vel_sigma_mps),
        Vec3d::Constant(std::pow(cfg_.init_att_sigma_deg * RAD, 2)),
        Vec3d::Constant(cfg_.init_accel_bias_sigma * cfg_.init_accel_bias_sigma),
        Vec3d::Constant(cfg_.init_gyro_bias_sigma * cfg_.init_gyro_bias_sigma),
        cfg_.init_clock_bias_sigma_s * cfg_.init_clock_bias_sigma_s,
        cfg_.init_clock_drift_sigma_sps * cfg_.init_clock_drift_sigma_sps;
    MatXd P0 = p0.asDiagonal();
    Configure(0.0, r0, v0, R0, Vec3d::Zero(), Vec3d::Zero(), cb0, cd0, P0);

    // ---- 6. Results allocation ------------------------------------------------
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
    res_.vel_err_norm = VecXd::Zero(N_);
    res_.clock_bias_err = VecXd::Zero(N_);
    res_.clock_bias_sigma = VecXd::Zero(N_);
    res_.n_visible_sat = VecXi::Zero(N_);
    res_.n_craters = VecXi::Zero(N_);
    res_.accel_bias_err = MatXd::Zero(N_, 3);
    res_.accel_bias_sigma = MatXd::Zero(N_, 3);
    res_.gyro_bias_err = MatXd::Zero(N_, 3);
    res_.gyro_bias_sigma = MatXd::Zero(N_, 3);
    res_.att_err_deg = MatXd::Zero(N_, 3);
    res_.att_sigma_deg = MatXd::Zero(N_, 3);
    res_.att_quat_est = MatXd::Zero(N_, 4);
    res_.traj_enu_truth = MatXd::Zero(N_, 3);
    res_.traj_enu_est = MatXd::Zero(N_, 3);
    res_.alt_truth = VecXd::Zero(N_);
    res_.alt_est = VecXd::Zero(N_);
    for (int j = 0; j < n_sat_; ++j) res_.satellite_names.push_back(cfg_.satellites[j].name);

    auto* lander = dynamic_cast<AgentWithDynamics*>(agent_);
    if (lander) {
      Vec6 rv;
      rv << r_truth_[0].cast<Real>(), v_truth_[0].cast<Real>();
      lander->SetTime(0.0);
      lander->SetState(Cart6(rv, Frame::MOON_PA));
    }
    LogEpoch(0);
  }

  void LanderNavApp::Step(Real t) {
    if (!self_driving_) return;
    if (!initialized_) {
      InitScenario();  // precompute descent truth/measurements/initial estimate on first Step
      initialized_ = true;
    }
    int k = static_cast<int>(std::lround(t.val() / dt_));
    if (k < 1 || k >= N_) return;
    const double dt = dt_;
    const double sqrt_dt = std::sqrt(dt);

    auto* lander = dynamic_cast<AgentWithDynamics*>(agent_);
    if (lander) {
      Vec6 rv;
      rv << r_truth_[k].cast<Real>(), v_truth_[k].cast<Real>();
      lander->SetTime(k * dt);
      lander->SetState(Cart6(rv, Frame::MOON_PA));
    }

    // Simulate the IMU (Kalibr: bias + white noise) over [k-1, k] and predict.
    SurfaceImuMeasurement imu;
    imu.accel
        = f_body_truth_[k - 1] + ba_truth_k_[k - 1] + Gauss3(cfg_.accel_noise_density / sqrt_dt);
    imu.gyro
        = w_body_truth_[k - 1] + bg_truth_k_[k - 1] + Gauss3(cfg_.gyro_noise_density / sqrt_dt);
    Predict(imu, dt);

    // LunaNet (LANS) pseudorange updates for visible satellites.
    if (cfg_.enable_lunanet) {
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
      res_.n_visible_sat(k) = static_cast<int>(lans.size());
      UpdateLans(lans);
    }

    // Radar altimeter: height above the terrain directly below (nadir), if within range.
    World* world = agent_->GetWorld();
    if (cfg_.enable_altimeter && Alt_[k] <= cfg_.altimeter_max_range_m) {
      LanderAltimeterMeasurement m;
      m.sigma_m = cfg_.altimeter_sigma_m;
      m.altitude_m = Alt_[k] + Gauss(cfg_.altimeter_sigma_m);
      Vec3d enu = R_pa2enu_ * (r_ - r_center_pa_);
      double E = enu(0), Ne = enu(1), U = enu(2);
      double alt_pred = U - world->GetElevation(E, Ne);
      double sE = (world->GetElevation(E + ds_, Ne) - world->GetElevation(E - ds_, Ne)) / (2 * ds_);
      double sN = (world->GetElevation(E, Ne + ds_) - world->GetElevation(E, Ne - ds_)) / (2 * ds_);
      Vec3d H_pos = (R_pa2enu_.row(2) - sE * R_pa2enu_.row(0) - sN * R_pa2enu_.row(1)).transpose();
      m.predicted_altitude_m = alt_pred;
      m.h_pos = H_pos;
      UpdateAltimeter(m);
    }

    // Crater-landmark bearings: down-looking camera within a nadir FOV cone.
    if (cfg_.enable_craters && n_crat_ > 0) {
      double cos_fov = std::cos(cfg_.camera_fov_deg * RAD);
      std::vector<std::pair<double, int>> cands;
      for (int j = 0; j < n_crat_; ++j) {
        Vec3d los = crater_pa_[j] - r_truth_[k];
        double range = los.norm();
        if (range <= 0.0) continue;
        Vec3d u_n = los / range;
        double cos_nadir = -u_n.dot(up_hat_pa_);
        if (cos_nadir < cos_fov) continue;
        cands.emplace_back(range, j);
      }
      std::sort(cands.begin(), cands.end());
      int n_use = std::min(static_cast<int>(cands.size()), cfg_.max_craters_per_epoch);
      std::vector<LanderCraterMeasurement> craters;
      for (int i = 0; i < n_use; ++i) {
        int j = cands[i].second;
        Vec3d los = crater_pa_[j] - r_truth_[k];
        Vec3d u_n = los / los.norm();
        Vec3d u_body = R_truth_[k].transpose() * u_n + Gauss3(crater_sigma_rad_);
        LanderCraterMeasurement m;
        m.r_crater = crater_pa_[j];
        m.los_body = u_body.normalized();
        m.sigma_rad = crater_sigma_rad_;
        craters.push_back(m);
      }
      res_.n_craters(k) = static_cast<int>(craters.size());
      UpdateCraters(craters);
    }

    LogEpoch(k);
  }

  void LanderNavApp::LogEpoch(int k) {
    World* world = agent_->GetWorld();
    const MatXd& P = P_;
    Vec3d r_est = r_;
    Vec3d v_est = v_;
    Vec3d err_pa = r_truth_[k] - r_est;
    Vec3d err_enu = R_pa2enu_ * err_pa;
    res_.time_s(k) = k * dt_;
    res_.pos_err_enu.row(k) = err_enu.transpose();
    res_.pos_err_norm(k) = err_pa.norm();
    res_.vel_err_norm(k) = (v_truth_[k] - v_est).norm();
    Mat3d P_enu = R_pa2enu_ * P.block<3, 3>(0, 0) * R_enu2pa_;
    for (int i = 0; i < 3; ++i) res_.pos_sigma_enu(k, i) = std::sqrt(std::max(0.0, P_enu(i, i)));
    res_.clock_bias_err(k) = ClockBiasTruth(k) - cb_;
    res_.clock_bias_sigma(k) = std::sqrt(std::max(0.0, P(15, 15)));

    Vec3d ba_err = ba_truth_k_[k] - ba_;
    Vec3d bg_err = bg_truth_k_[k] - bg_;
    res_.accel_bias_err.row(k) = ba_err.transpose();
    res_.gyro_bias_err.row(k) = bg_err.transpose();
    Vec3d att_err = LogSO3(R_truth_[k] * q_b2n_.toRotationMatrix().transpose());
    res_.att_err_deg.row(k) = (att_err / RAD).transpose();
    res_.att_quat_est.row(k) = quaternion().transpose();
    for (int i = 0; i < 3; ++i) {
      res_.accel_bias_sigma(k, i) = std::sqrt(std::max(0.0, P(9 + i, 9 + i)));
      res_.gyro_bias_sigma(k, i) = std::sqrt(std::max(0.0, P(12 + i, 12 + i)));
      res_.att_sigma_deg(k, i) = std::sqrt(std::max(0.0, P(6 + i, 6 + i))) / RAD;
    }

    Vec3d enu_est = R_pa2enu_ * (r_est - r_center_pa_);
    res_.traj_enu_truth.row(k) = Vec3d(Ee_[k], Nn_[k], Uu_[k]).transpose();
    res_.traj_enu_est.row(k) = enu_est.transpose();
    res_.alt_truth(k) = Alt_[k];
    res_.alt_est(k) = enu_est(2) - world->GetElevation(enu_est(0), enu_est(1));
  }

  void LanderNavApp::Log(Real /*t*/) {}

  REGISTER_FACTORY_CLASS(Application, LanderNavApp)

}  // namespace lupnt

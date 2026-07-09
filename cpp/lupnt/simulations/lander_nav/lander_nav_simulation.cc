#include "lupnt/simulations/lander_nav/lander_nav_simulation.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "lupnt/agents/lander.h"
#include "lupnt/applications/lander_nav_app.h"
#include "lupnt/interfaces/lola_dem.h"
#include "lupnt/lupnt.h"
#include "lupnt/measurements/lander_measurements.h"
#include "lupnt/measurements/surface_measurements.h"

namespace lupnt {

  namespace {
    Vec3d Gravity(const Vec3d& r) {
      double rn = r.norm();
      return (rn > 0.0) ? Vec3d(-GM_MOON / (rn * rn * rn) * r) : Vec3d(Vec3d::Zero());
    }

    Mat3d Skew3d(const Vec3d& v) {
      Mat3d S;
      S << 0.0, -v(2), v(1), v(2), 0.0, -v(0), -v(1), v(0), 0.0;
      return S;
    }

    // Rotation-vector (log) of a small rotation matrix, in the same convention as ExpSO3.
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
  }  // namespace

  LanderNavResults RunLanderNav(const LanderNavConfig& cfg) {
    std::mt19937 rng(cfg.seed);
    std::normal_distribution<double> nd(0.0, 1.0);
    std::uniform_real_distribution<double> ud(0.0, 1.0);
    auto gauss = [&](double sigma) { return sigma * nd(rng); };
    auto gauss3 = [&](double sigma) { return Vec3d(gauss(sigma), gauss(sigma), gauss(sigma)); };

    const int n_sat = static_cast<int>(cfg.satellites.size());
    // A supplied reference trajectory (N x 3 ENU) overrides the built-in smoothstep and sets N.
    const bool use_ref_traj = cfg.ref_traj_enu.rows() > 0;
    LUPNT_CHECK(!use_ref_traj || cfg.ref_traj_enu.cols() == 3,
                "ref_traj_enu must have 3 columns (East, North, Up)", "RunLanderNav");
    const int N = use_ref_traj
                      ? static_cast<int>(cfg.ref_traj_enu.rows())
                      : std::max(2, static_cast<int>(std::round(cfg.duration_s / cfg.dt_s)) + 1);
    const double dt = cfg.dt_s;
    const double sqrt_dt = std::sqrt(dt);
    const double crater_sigma_rad = cfg.crater_sigma_arcsec * RAD / 3600.0;

    // ---- 1. Load the DEM and build the local ENU terrain frame ----------------
    LunarDem dem
        = LoadLolaDem(cfg.site_lat_deg, cfg.site_lon_deg, cfg.dem_half_width_m, cfg.dem_max_res_m);
    const double cx = dem.center_x();
    const double cy = dem.center_y();

    Vec3 r_center_real = LatLonAltToCart(Vec3(cfg.site_lat_deg, cfg.site_lon_deg, 0.0), R_MOON);
    Cart3 r_center_state(r_center_real, Frame::MOON_PA);
    Mat3d R_enu2pa = RotEastNorthUpToCart(r_center_state, R_MOON).cast<double>();
    Mat3d R_pa2enu = R_enu2pa.transpose();
    Vec3d r_center_pa = r_center_real.cast<double>();
    Vec3d up_hat_pa = R_enu2pa.col(2);

    auto EnuToPa = [&](double E, double Nn, double U) -> Vec3d {
      return r_center_pa + R_enu2pa * Vec3d(E, Nn, U);
    };
    auto Elevation
        = [&](double E, double Nn) -> double { return dem.GetElevation(cx + E, cy + Nn); };

    // Body(vehicle)-to-nav rotation for heading `hdg` [rad]: x=forward, y=left, z=up.
    auto BodyToPa = [&](double hdg) -> Mat3d {
      Vec3d xb(std::cos(hdg), std::sin(hdg), 0.0);
      Vec3d zb(0.0, 0.0, 1.0);
      Vec3d yb = zb.cross(xb);
      Mat3d R_body2enu;
      R_body2enu.col(0) = xb;
      R_body2enu.col(1) = yb;
      R_body2enu.col(2) = zb;
      return R_enu2pa * R_body2enu;
    };

    // ---- 2. Synthetic crater-landmark map (scattered on the terrain) ----------
    const int n_crat = cfg.enable_craters ? std::max(0, cfg.n_craters) : 0;
    std::vector<Vec3d> crater_pa(n_crat);
    MatXd crater_enu = MatXd::Zero(n_crat, 2);
    for (int j = 0; j < n_crat; ++j) {
      double ang = 2.0 * M_PI * ud(rng);
      double rad = cfg.crater_field_radius_m * std::sqrt(ud(rng));
      double E = rad * std::cos(ang), Nn = rad * std::sin(ang);
      double U = Elevation(E, Nn);
      crater_pa[j] = EnuToPa(E, Nn, U);
      crater_enu(j, 0) = E;
      crater_enu(j, 1) = Nn;
    }

    // ---- 3. Descent truth: ENU path, Moon-fixed position/velocity/accel, attitude/rate ----
    const double hdg = cfg.descent_heading_deg * RAD;
    std::vector<double> Ee(N), Nn(N), Uu(N), Alt(N);
    std::vector<Vec3d> r_truth(N), v_truth(N);
    std::vector<Mat3d> R_truth(N);
    for (int k = 0; k < N; ++k) {
      if (use_ref_traj) {
        // Externally-supplied ENU reference path (U is height above the site datum).
        Ee[k] = cfg.ref_traj_enu(k, 0);
        Nn[k] = cfg.ref_traj_enu(k, 1);
        Uu[k] = cfg.ref_traj_enu(k, 2);
        Alt[k] = Uu[k] - Elevation(Ee[k], Nn[k]);  // height above the terrain below
      } else {
        // Built-in smoothstep descent (terrain-following altitude).
        double tau = (cfg.duration_s > 0.0) ? (k * dt / cfg.duration_s) : 1.0;
        double s = SmoothStep(tau);
        Ee[k] = cfg.descent_start_east_m + (cfg.descent_end_east_m - cfg.descent_start_east_m) * s;
        Nn[k]
            = cfg.descent_start_north_m + (cfg.descent_end_north_m - cfg.descent_start_north_m) * s;
        Alt[k] = cfg.descent_start_alt_m + (cfg.descent_end_alt_m - cfg.descent_start_alt_m) * s;
        Uu[k] = Elevation(Ee[k], Nn[k]) + Alt[k];
      }
      R_truth[k] = BodyToPa(hdg);
    }
    for (int k = 0; k < N; ++k) r_truth[k] = EnuToPa(Ee[k], Nn[k], Uu[k]);
    for (int k = 0; k < N; ++k) {
      int kp = std::min(k + 1, N - 1), km = std::max(k - 1, 0);
      double span = (kp - km) * dt;
      v_truth[k] = (span > 0.0) ? Vec3d((r_truth[kp] - r_truth[km]) / span) : Vec3d(Vec3d::Zero());
    }

    // Truth specific force (body) and angular rate (body) at each epoch.
    std::vector<Vec3d> f_body_truth(N), w_body_truth(N);
    for (int k = 0; k < N; ++k) {
      int kp = std::min(k + 1, N - 1), km = std::max(k - 1, 0);
      double span = (kp - km) * dt;
      Vec3d a_total
          = (span > 0.0) ? Vec3d((v_truth[kp] - v_truth[km]) / span) : Vec3d(Vec3d::Zero());
      f_body_truth[k] = R_truth[k].transpose() * (a_total - Gravity(r_truth[k]));
      Mat3d Rdot = (span > 0.0) ? Mat3d((R_truth[kp] - R_truth[km]) / span) : Mat3d(Mat3d::Zero());
      Mat3d Wx = R_truth[k].transpose() * Rdot;  // [w x] in body
      w_body_truth[k] = Vec3d(Wx(2, 1), Wx(0, 2), Wx(1, 0));
    }

    // Truth IMU biases: constant offset + random walk (Kalibr).
    Vec3d ba_truth = gauss3(cfg.accel_bias0);
    Vec3d bg_truth = gauss3(cfg.gyro_bias0);
    std::vector<Vec3d> ba_truth_k(N), bg_truth_k(N);
    for (int k = 0; k < N; ++k) {
      if (k > 0) {
        ba_truth += gauss3(cfg.accel_bias_rw * sqrt_dt);
        bg_truth += gauss3(cfg.gyro_bias_rw * sqrt_dt);
      }
      ba_truth_k[k] = ba_truth;
      bg_truth_k[k] = bg_truth;
    }

    // ---- 4. Propagate LunaNet truth orbits (Keplerian) ------------------------
    const Real t0_tdb = ConvertTime(GregorianToTime(cfg.start_epoch_utc), Time::UTC, Time::TDB);
    CartesianTwoBodyDynamics sat_dyn(GM_MOON);
    sat_dyn.SetTimeStep(dt);
    std::vector<Vec6> sat_ci(n_sat);
    for (int j = 0; j < n_sat; ++j)
      sat_ci[j]
          = (Vec6() << cfg.satellites[j].r0_m.cast<Real>(), cfg.satellites[j].v0_mps.cast<Real>())
                .finished();
    std::vector<std::vector<Vec3d>> sat_pa(N, std::vector<Vec3d>(n_sat));
    for (int k = 0; k < N; ++k) {
      Real tk = t0_tdb + k * dt;
      if (k > 0) {
        Real tkm1 = t0_tdb + (k - 1) * dt;
        for (int j = 0; j < n_sat; ++j) {
          Cart6 st = sat_dyn.Propagate(Cart6(sat_ci[j], Frame::MOON_CI), tkm1, tk, nullptr);
          sat_ci[j] = st.head(6);
        }
      }
      for (int j = 0; j < n_sat; ++j) {
        Vec6 rv_pa = ConvertFrame(tk, sat_ci[j], Frame::MOON_CI, Frame::MOON_PA);
        sat_pa[k][j] = rv_pa.head(3).cast<double>();
      }
    }

    std::vector<double> sise_bias(n_sat);
    for (int j = 0; j < n_sat; ++j) sise_bias[j] = gauss(cfg.sise_m);

    // ---- 5. Build the Lander agent and attach the navigation application ------
    LanderNavAppParams params;
    params.accel_noise_density = cfg.accel_noise_density;
    params.accel_bias_rw = cfg.filter_bias_rw_scale * cfg.accel_bias_rw;
    params.gyro_noise_density = cfg.gyro_noise_density;
    params.gyro_bias_rw = cfg.filter_bias_rw_scale * cfg.gyro_bias_rw;
    params.pseudorange_sigma_m = cfg.pseudorange_sigma_m;
    params.sise_m = cfg.sise_m;
    params.altimeter_sigma_m = cfg.altimeter_sigma_m;
    params.crater_sigma_rad = crater_sigma_rad;

    Lander lander;
    lander.SetName("Lander");
    auto nav = std::make_shared<LanderNavApp>(params);
    lander.SetApplication(nav);  // wires nav->GetAgent() == &lander (app hosted on the agent)

    auto ClockBiasTruth
        = [&](int k) { return cfg.lander_clock_bias_s + cfg.lander_clock_drift_sps * (k * dt); };

    // Initial estimate: truth + perturbations; IMU biases start at zero (unknown).
    Vec3d r0 = r_truth[0] + gauss3(cfg.init_pos_sigma_m);
    Vec3d v0 = v_truth[0] + gauss3(cfg.init_vel_sigma_mps);
    Vec3d att_err0 = gauss3(cfg.init_att_sigma_deg * RAD);
    // R_true = Exp(dtheta) R_est  =>  R_est = Exp(-dtheta) R_true.
    Mat3d dR0 = Mat3d::Identity() - Skew3d(att_err0);  // small-angle Exp(-att_err0)
    Mat3d R0 = dR0 * R_truth[0];
    double cb0 = ClockBiasTruth(0) + gauss(cfg.init_clock_bias_sigma_s);
    double cd0 = cfg.lander_clock_drift_sps + gauss(cfg.init_clock_drift_sigma_sps);

    VecXd p0(kLanderNavErrorStateSize);
    p0 << Vec3d::Constant(cfg.init_pos_sigma_m * cfg.init_pos_sigma_m),
        Vec3d::Constant(cfg.init_vel_sigma_mps * cfg.init_vel_sigma_mps),
        Vec3d::Constant(std::pow(cfg.init_att_sigma_deg * RAD, 2)),
        Vec3d::Constant(cfg.init_accel_bias_sigma * cfg.init_accel_bias_sigma),
        Vec3d::Constant(cfg.init_gyro_bias_sigma * cfg.init_gyro_bias_sigma),
        cfg.init_clock_bias_sigma_s * cfg.init_clock_bias_sigma_s,
        cfg.init_clock_drift_sigma_sps * cfg.init_clock_drift_sigma_sps;
    MatXd P0 = p0.asDiagonal();
    nav->Configure(0.0, r0, v0, R0, Vec3d::Zero(), Vec3d::Zero(), cb0, cd0, P0);

    // ---- 6. Results allocation ------------------------------------------------
    LanderNavResults res;
    res.site_id = dem.site().id;
    res.site_name = dem.site().name;
    res.site_lat_deg = dem.site().lat_deg;
    res.site_lon_deg = dem.site().lon_deg;
    res.dem_x = dem.x();
    res.dem_y = dem.y();
    res.dem_elevation = dem.elevation();
    res.dem_center_x = cx;
    res.dem_center_y = cy;
    res.crater_enu = crater_enu;
    res.time_s = VecXd::Zero(N);
    res.pos_err_enu = MatXd::Zero(N, 3);
    res.pos_sigma_enu = MatXd::Zero(N, 3);
    res.pos_err_norm = VecXd::Zero(N);
    res.vel_err_norm = VecXd::Zero(N);
    res.clock_bias_err = VecXd::Zero(N);
    res.clock_bias_sigma = VecXd::Zero(N);
    res.n_visible_sat = VecXi::Zero(N);
    res.n_craters = VecXi::Zero(N);
    res.accel_bias_err = MatXd::Zero(N, 3);
    res.accel_bias_sigma = MatXd::Zero(N, 3);
    res.gyro_bias_err = MatXd::Zero(N, 3);
    res.gyro_bias_sigma = MatXd::Zero(N, 3);
    res.att_err_deg = MatXd::Zero(N, 3);
    res.att_sigma_deg = MatXd::Zero(N, 3);
    res.att_quat_est = MatXd::Zero(N, 4);
    res.traj_enu_truth = MatXd::Zero(N, 3);
    res.traj_enu_est = MatXd::Zero(N, 3);
    res.alt_truth = VecXd::Zero(N);
    res.alt_est = VecXd::Zero(N);
    for (int j = 0; j < n_sat; ++j) res.satellite_names.push_back(cfg.satellites[j].name);

    const double ds = std::max(cfg.dem_max_res_m, 1.0);

    auto LogEpoch = [&](int k) {
      const MatXd& P = nav->covariance();
      Vec3d r_est = nav->position();
      Vec3d v_est = nav->velocity();
      Vec3d err_pa = r_truth[k] - r_est;
      Vec3d err_enu = R_pa2enu * err_pa;
      res.time_s(k) = k * dt;
      res.pos_err_enu.row(k) = err_enu.transpose();
      res.pos_err_norm(k) = err_pa.norm();
      res.vel_err_norm(k) = (v_truth[k] - v_est).norm();
      Mat3d P_enu = R_pa2enu * P.block<3, 3>(0, 0) * R_enu2pa;
      for (int i = 0; i < 3; ++i) res.pos_sigma_enu(k, i) = std::sqrt(std::max(0.0, P_enu(i, i)));
      res.clock_bias_err(k) = ClockBiasTruth(k) - nav->clock_bias();
      res.clock_bias_sigma(k) = std::sqrt(std::max(0.0, P(15, 15)));

      Vec3d ba_err = ba_truth_k[k] - nav->accel_bias();
      Vec3d bg_err = bg_truth_k[k] - nav->gyro_bias();
      res.accel_bias_err.row(k) = ba_err.transpose();
      res.gyro_bias_err.row(k) = bg_err.transpose();
      Vec3d att_err = LogSO3(R_truth[k] * nav->attitude().transpose());  // nav-frame error
      res.att_err_deg.row(k) = (att_err / RAD).transpose();
      res.att_quat_est.row(k) = nav->quaternion().transpose();  // MEKF nominal quaternion [w,x,y,z]
      for (int i = 0; i < 3; ++i) {
        res.accel_bias_sigma(k, i) = std::sqrt(std::max(0.0, P(9 + i, 9 + i)));
        res.gyro_bias_sigma(k, i) = std::sqrt(std::max(0.0, P(12 + i, 12 + i)));
        res.att_sigma_deg(k, i) = std::sqrt(std::max(0.0, P(6 + i, 6 + i))) / RAD;
      }

      Vec3d enu_est = R_pa2enu * (r_est - r_center_pa);
      res.traj_enu_truth.row(k) = Vec3d(Ee[k], Nn[k], Uu[k]).transpose();
      res.traj_enu_est.row(k) = enu_est.transpose();
      res.alt_truth(k) = Alt[k];
      res.alt_est(k) = enu_est(2) - Elevation(enu_est(0), enu_est(1));
    };

    // ---- 7. Filter loop -------------------------------------------------------
    LogEpoch(0);
    for (int k = 1; k < N; ++k) {
      // Advance the lander agent's truth state and drive its hosted nav app's IMU predict.
      Vec6 rv;
      rv << r_truth[k].cast<Real>(), v_truth[k].cast<Real>();
      lander.SetTime(k * dt);
      lander.SetState(Cart6(rv, Frame::MOON_PA));

      // Simulate the IMU (Kalibr: bias + white noise) over [k-1, k] and predict.
      SurfaceImuMeasurement imu;
      imu.accel
          = f_body_truth[k - 1] + ba_truth_k[k - 1] + gauss3(cfg.accel_noise_density / sqrt_dt);
      imu.gyro = w_body_truth[k - 1] + bg_truth_k[k - 1] + gauss3(cfg.gyro_noise_density / sqrt_dt);
      nav->Predict(imu, dt);

      // LunaNet (LANS) pseudorange updates for visible satellites.
      if (cfg.enable_lunanet) {
        std::vector<SurfaceLansMeasurement> lans;
        for (int j = 0; j < n_sat; ++j) {
          Vec3d los = sat_pa[k][j] - r_truth[k];
          double range = los.norm();
          double sin_el = (range > 0.0) ? (los.dot(up_hat_pa) / range) : -1.0;
          if (sin_el < std::sin(cfg.elevation_mask_deg * RAD)) continue;
          SurfaceLansMeasurement m;
          m.r_sat = sat_pa[k][j];
          m.sigma_m = cfg.pseudorange_sigma_m;
          m.sise_m = cfg.sise_m;
          m.pseudorange_m
              = range + C * ClockBiasTruth(k) + sise_bias[j] + gauss(cfg.pseudorange_sigma_m);
          lans.push_back(m);
        }
        res.n_visible_sat(k) = static_cast<int>(lans.size());
        nav->UpdateLans(lans);
      }

      // Radar altimeter: height above the terrain directly below (nadir), if within range.
      if (cfg.enable_altimeter && Alt[k] <= cfg.altimeter_max_range_m) {
        LanderAltimeterMeasurement m;
        m.sigma_m = cfg.altimeter_sigma_m;
        m.altitude_m = Alt[k] + gauss(cfg.altimeter_sigma_m);
        // Predicted altitude and Jacobian at the current estimate.
        Vec3d enu = R_pa2enu * (nav->position() - r_center_pa);
        double E = enu(0), Ne = enu(1), U = enu(2);
        double alt_pred = U - Elevation(E, Ne);
        double sE = (Elevation(E + ds, Ne) - Elevation(E - ds, Ne)) / (2 * ds);
        double sN = (Elevation(E, Ne + ds) - Elevation(E, Ne - ds)) / (2 * ds);
        Vec3d H_pos = (R_pa2enu.row(2) - sE * R_pa2enu.row(0) - sN * R_pa2enu.row(1)).transpose();
        m.predicted_altitude_m = alt_pred;
        m.h_pos = H_pos;
        nav->UpdateAltimeter(m);
      }

      // Crater-landmark bearings: down-looking camera within a nadir FOV cone.
      if (cfg.enable_craters && n_crat > 0) {
        double cos_fov = std::cos(cfg.camera_fov_deg * RAD);
        std::vector<std::pair<double, int>> cands;  // (range, crater index)
        for (int j = 0; j < n_crat; ++j) {
          Vec3d los = crater_pa[j] - r_truth[k];
          double range = los.norm();
          if (range <= 0.0) continue;
          Vec3d u_n = los / range;
          double cos_nadir = -u_n.dot(up_hat_pa);  // angle from nadir boresight
          if (cos_nadir < cos_fov) continue;       // outside the camera cone
          cands.emplace_back(range, j);
        }
        std::sort(cands.begin(), cands.end());
        int n_use = std::min(static_cast<int>(cands.size()), cfg.max_craters_per_epoch);
        std::vector<LanderCraterMeasurement> craters;
        for (int i = 0; i < n_use; ++i) {
          int j = cands[i].second;
          Vec3d los = crater_pa[j] - r_truth[k];
          Vec3d u_n = los / los.norm();
          Vec3d u_body = R_truth[k].transpose() * u_n + gauss3(crater_sigma_rad);
          LanderCraterMeasurement m;
          m.r_crater = crater_pa[j];
          m.los_body = u_body.normalized();
          m.sigma_rad = crater_sigma_rad;
          craters.push_back(m);
        }
        res.n_craters(k) = static_cast<int>(craters.size());
        nav->UpdateCraters(craters);
      }

      LogEpoch(k);
    }

    return res;
  }

}  // namespace lupnt

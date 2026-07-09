#include "lupnt/applications/lander_nav_app.h"

#include <cmath>

#include "lupnt/core/constants.h"
#include "lupnt/core/error.h"

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
    Vec3d a_nav = f_n + Gravity(r_);
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

}  // namespace lupnt

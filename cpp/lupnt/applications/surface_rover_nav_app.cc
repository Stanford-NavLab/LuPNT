#include "lupnt/applications/surface_rover_nav_app.h"

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

    Vec3d Gravity(const Vec3d& r) {
      double rn = r.norm();
      return (rn > 0.0) ? Vec3d(-GM_MOON / (rn * rn * rn) * r) : Vec3d(Vec3d::Zero());
    }
  }  // namespace

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
    Vec3d a_nav = f_n + Gravity(r_);
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

}  // namespace lupnt

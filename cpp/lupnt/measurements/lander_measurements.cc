#include "lupnt/measurements/lander_measurements.h"

namespace lupnt {

  namespace {
    Mat3d Skew3d(const Vec3d& v) {
      Mat3d S;
      S << 0.0, -v(2), v(1), v(2), 0.0, -v(0), -v(1), v(0), 0.0;
      return S;
    }
  }  // namespace

  MeasData LanderAltimeterMeasurement::Compute(const NavErrorContext& nom, MatXd* H) const {
    MeasData md;
    md.timestamp = timestamp;
    md.value = VecXd::Constant(1, predicted_altitude_m);
    md.covariance = MatXd::Constant(1, 1, NoiseVariance());
    if (H != nullptr) {
      *H = MatXd::Zero(1, nom.error_state_size);
      H->block(0, config.i_dr, 1, 3) = h_pos.transpose();  // d(altitude)/d(r), Moon-fixed frame
    }
    return md;
  }

  Vec3d LanderCraterMeasurement::PredictedLosBody(const Vec3d& r, const Mat3d& R_b2n) const {
    Vec3d los = r_crater - r;
    double range = los.norm();
    if (range <= 0.0) return Vec3d::Zero();
    return R_b2n.transpose() * (los / range);
  }

  MeasData LanderCraterMeasurement::Compute(const NavErrorContext& nom, MatXd* H) const {
    MeasData md;
    md.timestamp = timestamp;
    Vec3d los = r_crater - nom.r;
    double range = los.norm();
    if (range <= 0.0) return md;             // degenerate geometry: empty value -> caller skips
    Vec3d u_n = los / range;                 // nav-frame unit line-of-sight
    md.value = nom.R_b2n.transpose() * u_n;  // predicted body line-of-sight
    md.covariance = (sigma_rad * sigma_rad) * MatXd::Identity(3, 3);
    if (H != nullptr) {
      // d(u_n)/d(r) = -(1/range)(I - u_n u_n^T); attitude term uses [u_n]_x (nav frame).
      Mat3d dun_dr = -(1.0 / range) * (Mat3d::Identity() - u_n * u_n.transpose());
      *H = MatXd::Zero(3, nom.error_state_size);
      H->block(0, config.i_dr, 3, 3) = nom.R_b2n.transpose() * dun_dr;
      H->block(0, config.i_dth, 3, 3) = nom.R_b2n.transpose() * Skew3d(u_n);
    }
    return md;
  }

}  // namespace lupnt

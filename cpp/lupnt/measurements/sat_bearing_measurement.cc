#include "lupnt/measurements/sat_bearing_measurement.h"

namespace lupnt {

  VecX SatBearingMeasurement::Model(const VecX& x) const {
    const Vec3 r_obs = x.segment(config_.idx_position, 3);
    const Vec3 target = config_.target_pos_mci.cast<Real>();
    const Vec3 dr = target - r_obs;
    const Vec3 u = dr / dr.norm();
    VecX y(3);
    y = u;
    return y;
  }

  MatXd SatBearingMeasurement::Covariance() const {
    return (config_.sigma_rad * config_.sigma_rad) * MatXd::Identity(NumRows(), NumRows());
  }

  MeasData SatBearingMeasurement::Compute(const State& x, MatXd* H) const {
    MeasData md;
    VecX y;
    if (H != nullptr) {
      // Reseed a fresh autodiff state and differentiate the raw-VecX model (same pattern as
      // `IslCrosslinkMeasurement`), keeping the derivative seeds intact for the EKF update.
      VecX x_tmp = x.cast<double>();
      auto f = [this](const VecX& xx) -> VecX { return Model(xx); };
      jacobian(f, wrt(x_tmp), at(x_tmp), y, *H);
    } else {
      y = Model(x);
    }
    md.value = y.cast<double>();
    md.covariance = Covariance();
    return md;
  }

}  // namespace lupnt

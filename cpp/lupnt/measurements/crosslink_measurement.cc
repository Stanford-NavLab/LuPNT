#include "lupnt/measurements/crosslink_measurement.h"

#include "lupnt/core/constants.h"
#include "lupnt/measurements/measurement.h"
#include "lupnt/measurements/measurement_utils.h"

namespace lupnt {

  int IslCrosslinkMeasurement::NumRows() const {
    return 2 * config_.n_links + static_cast<int>(config_.anchor_pos_mci.size());
  }

  VecX IslCrosslinkMeasurement::Model(const VecX& x) const {
    const int n_links = config_.n_links;
    const int sub = config_.sub_state_size;
    const int n_anchor = static_cast<int>(config_.anchor_pos_mci.size());
    const int m = 2 * n_links + n_anchor;

    VecX y(m);
    VecX r0 = x.segment(config_.idx_position, 3);
    VecX v0 = x.segment(config_.idx_velocity, 3);
    for (int i = 0; i < n_links; ++i) {
      const int off = sub * (i + 1);
      Vec2 yi = RangeAndRangeRate(r0, VecX(x.segment(off, 3)), v0, VecX(x.segment(off + 3, 3)));
      y(2 * i) = yi(0);
      y(2 * i + 1) = yi(1);
    }

    const Vec3 r_hub = x.segment(config_.idx_position, 3);
    const Real b_hub = x(config_.idx_clock_bias);  // hub clock bias [s]
    for (int ad = 0; ad < n_anchor; ++ad) {
      const Vec3 ra = config_.anchor_pos_mci[ad].cast<Real>();
      y(2 * n_links + ad) = (r_hub - ra).norm() + C * b_hub;
    }
    return y;
  }

  MatXd IslCrosslinkMeasurement::Covariance() const {
    const int n_links = config_.n_links;
    const int n_anchor = static_cast<int>(config_.anchor_pos_mci.size());
    MatXd R = MatXd::Zero(NumRows(), NumRows());
    for (int i = 0; i < n_links; ++i) {
      R(2 * i, 2 * i) = config_.sigma_range_m * config_.sigma_range_m;
      R(2 * i + 1, 2 * i + 1) = config_.sigma_range_rate_mps * config_.sigma_range_rate_mps;
    }
    for (int ad = 0; ad < n_anchor; ++ad) {
      R(2 * n_links + ad, 2 * n_links + ad)
          = config_.sigma_pseudorange_m * config_.sigma_pseudorange_m;
    }
    return R;
  }

  MeasData IslCrosslinkMeasurement::Compute(const State& x, MatXd* H) const {
    MeasData md;
    VecX y;
    if (H != nullptr) {
      // Reseed a fresh autodiff state (drops incoming derivative seeds) and differentiate
      // the raw-VecX model, matching the original inline behavior.
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

#include "lupnt/measurements/crosslink_measurement.h"

#include "lupnt/core/constants.h"
#include "lupnt/measurements/measurement.h"
#include "lupnt/measurements/measurement_utils.h"

namespace lupnt {

  int IslCrosslinkMeasurement::NumRows() const {
    const int n_anchor = static_cast<int>(config_.anchor_pos_mci.size());
    return 2 * config_.n_links + n_anchor + (config_.include_anchor_doppler ? n_anchor : 0)
           + (config_.include_time_transfer ? config_.n_links : 0)
           + (config_.include_frequency_transfer ? config_.n_links : 0);
  }

  VecX IslCrosslinkMeasurement::Model(const VecX& x) const {
    const int n_links = config_.n_links;
    const int sub = config_.sub_state_size;
    const int n_anchor = static_cast<int>(config_.anchor_pos_mci.size());
    const int n_ad = config_.include_anchor_doppler ? n_anchor : 0;
    const int n_tt = config_.include_time_transfer ? n_links : 0;
    const int n_ft = config_.include_frequency_transfer ? n_links : 0;
    const int m = 2 * n_links + n_anchor + n_ad + n_tt + n_ft;

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
    const Vec3 v_hub = x.segment(config_.idx_velocity, 3);
    const Real b_hub = x(config_.idx_clock_bias);   // hub clock bias [s]
    const Real d_hub = x(config_.idx_clock_drift);  // hub clock drift [s/s]
    for (int ad = 0; ad < n_anchor; ++ad) {
      const Vec3 ra = config_.anchor_pos_mci[ad].cast<Real>();
      y(2 * n_links + ad) = (r_hub - ra).norm() + C * b_hub;
    }

    // One-way anchor Doppler (pseudorange-rate) rows: u . (v_hub - v_anchor) + C * d_hub.
    for (int ad = 0; ad < n_ad; ++ad) {
      const Vec3 ra = config_.anchor_pos_mci[ad].cast<Real>();
      const Vec3 va = config_.anchor_vel_mci[ad].cast<Real>();
      const Vec3 dr = r_hub - ra;
      const Vec3 u = dr / dr.norm();
      y(2 * n_links + n_anchor + ad) = u.dot(v_hub - va) + C * d_hub;
    }

    // Two-way time-transfer rows: range-equivalent clock-bias difference between the hub
    // and each linked satellite, C * (b_hub - b_link_i).
    for (int i = 0; i < n_tt; ++i) {
      const Real b_link = x(sub * (i + 1) + config_.idx_clock_bias);
      y(2 * n_links + n_anchor + n_ad + i) = C * (b_hub - b_link);
    }

    // Two-way frequency-transfer rows: range-rate-equivalent clock-drift difference,
    // C * (d_hub - d_link_i).
    for (int i = 0; i < n_ft; ++i) {
      const Real d_link = x(sub * (i + 1) + config_.idx_clock_drift);
      y(2 * n_links + n_anchor + n_ad + n_tt + i) = C * (d_hub - d_link);
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
    const int n_ad = config_.include_anchor_doppler ? n_anchor : 0;
    for (int ad = 0; ad < n_ad; ++ad) {
      const int r = 2 * n_links + n_anchor + ad;
      R(r, r) = config_.sigma_anchor_doppler_mps * config_.sigma_anchor_doppler_mps;
    }
    const int n_tt = config_.include_time_transfer ? n_links : 0;
    for (int i = 0; i < n_tt; ++i) {
      const int r = 2 * n_links + n_anchor + n_ad + i;
      R(r, r) = config_.sigma_time_transfer_m * config_.sigma_time_transfer_m;
    }
    const int n_ft = config_.include_frequency_transfer ? n_links : 0;
    for (int i = 0; i < n_ft; ++i) {
      const int r = 2 * n_links + n_anchor + n_ad + n_tt + i;
      R(r, r) = config_.sigma_frequency_transfer_mps * config_.sigma_frequency_transfer_mps;
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

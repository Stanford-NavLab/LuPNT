#include "lupnt/measurements/ground_range_measurement.h"

namespace lupnt {

  MeasData GroundStationRangeMeasurement::Compute(const State& x, MatXd* H) const {
    const int n_rows = NumRows();
    Vec6d xs = x.head(6).cast<double>();
    const Vec6d& xgs = config_.reference_state;
    Vec3d dr = xs.head(3) - xgs.head(3);
    Vec3d dv = xs.tail(3) - xgs.tail(3);
    double rho = dr.norm();
    Vec3d u = dr / rho;
    double rho_dot = dr.dot(dv) / rho;

    VecXd y(n_rows);
    MatXd h_obs = MatXd::Zero(n_rows, 6);
    MatXd R = MatXd::Zero(n_rows, n_rows);
    int row = 0;
    if (config_.use_range) {
      h_obs.block(row, 0, 1, 3) = u.transpose();
      y(row) = rho;
      R(row, row) = config_.range_sigma_m * config_.range_sigma_m;
      row++;
    }
    if (config_.use_range_rate) {
      h_obs.block(row, 0, 1, 3) = ((dv - rho_dot * u) / rho).transpose();
      h_obs.block(row, 3, 1, 3) = u.transpose();
      y(row) = rho_dot;
      R(row, row) = config_.range_rate_sigma_mps * config_.range_rate_sigma_mps;
      row++;
    }

    MeasData md;
    md.value = y;
    md.covariance = R;
    if (H != nullptr) *H = h_obs;
    return md;
  }

}  // namespace lupnt

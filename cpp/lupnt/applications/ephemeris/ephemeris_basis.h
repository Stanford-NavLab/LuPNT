#pragma once
/**
 * @file ephemeris_basis.h
 * @brief Internal numerics helpers shared by `CartesianEphemeris` (ephemeris.cc)
 * and `Almanac` (almanac.cc): a Chebyshev polynomial basis (and its time
 * derivative), angle unwrapping, cumulative-trapezoid integration, the
 * Radial/Transverse/Normal (RTN) rotation matrix, and the RTN-decomposed
 * RMS/95th-percentile fit-error statistics both classes report from `EvalError`.
 *
 * Not part of the public API (no Python bindings) -- purely a private
 * implementation detail of ephemeris.cc and almanac.cc in lupnt/applications.
 */

#include <algorithm>
#include <cmath>

#include "lupnt/applications/ephemeris/lunanet_ephemeris.h"
#include "lupnt/core/constants.h"
#include "lupnt/core/definitions.h"

namespace lupnt {

  /// @brief Chebyshev polynomial basis (T_n basis, constant term first),
  /// evaluated at normalized time z_k = 2*t_k/t_fit.
  /// @return An (N x (order+1)) matrix whose column j is T_j(z_k).
  inline MatXd ChebyshevBasis(const VecXd& t_k, double t_fit, int order) {
    const int n = static_cast<int>(t_k.size());
    MatXd T(n, order + 1);
    T.col(0).setOnes();
    if (order == 0) return T;
    VecXd z = (2.0 / t_fit) * t_k;
    T.col(1) = z;
    for (int j = 2; j <= order; ++j) {
      T.col(j) = 2.0 * z.array() * T.col(j - 1).array() - T.col(j - 2).array();
    }
    return T;
  }

  /// @brief Time derivative (d/dt) of `ChebyshevBasis(t_k, t_fit, order)`.
  inline MatXd ChebyshevBasisDt(const VecXd& t_k, double t_fit, int order) {
    const int n = static_cast<int>(t_k.size());
    MatXd T = ChebyshevBasis(t_k, t_fit, order);
    MatXd Tdot(n, order + 1);
    Tdot.col(0).setZero();
    if (order == 0) return Tdot;
    const double dzdt = 2.0 / t_fit;
    VecXd z = (2.0 / t_fit) * t_k;
    Tdot.col(1).setConstant(dzdt);
    for (int j = 2; j <= order; ++j) {
      Tdot.col(j) = 2.0 * z.array() * Tdot.col(j - 1).array() + 2.0 * T.col(j - 1).array() * dzdt
                    - Tdot.col(j - 2).array();
    }
    return Tdot;
  }

  /// @brief Unwraps a sequence of angles [rad] so consecutive samples differ by
  /// less than pi, adding multiples of 2*pi (like `numpy.unwrap`). Used to make an
  /// angular orbital-element time series (e.g. mean anomaly) safe to fit with a
  /// polynomial.
  inline VecXd UnwrapAngles(const VecXd& angles) {
    VecXd out = angles;
    for (int i = 1; i < out.size(); ++i) {
      while (out(i) - out(i - 1) > PI) out(i) -= 2.0 * PI;
      while (out(i) - out(i - 1) < -PI) out(i) += 2.0 * PI;
    }
    return out;
  }

  /// @brief Cumulative trapezoidal integral of samples `y(t)`: `out(0) = 0`,
  /// `out(k)` = integral of y from `t(0)` to `t(k)`.
  inline VecXd CumTrapz(const VecXd& y, const VecXd& t) {
    VecXd out = VecXd::Zero(y.size());
    for (int i = 1; i < y.size(); ++i) {
      out(i) = out(i - 1) + 0.5 * (y(i) + y(i - 1)) * (t(i) - t(i - 1));
    }
    return out;
  }

  /// @brief Radial-Transverse-Normal (RTN) rotation matrix at Cartesian state
  /// `rv` (position `rv.head(3)`, velocity `rv.tail(3)`): rows are the unit
  /// vectors R (radial, along position), T (transverse, completes the right-handed
  /// triad), N (orbit-normal, along r x v). `v_rtn = RtnMatrix(rv) * v_inertial`.
  inline Mat3d RtnMatrix(const Vec6d& rv) {
    const Vec3d r = rv.head<3>();
    const Vec3d v = rv.tail<3>();
    const Vec3d R = r.normalized();
    const Vec3d N = (r.cross(v)).normalized();
    const Vec3d T = N.cross(R);
    Mat3d M;
    M.row(0) = R.transpose();
    M.row(1) = T.transpose();
    M.row(2) = N.transpose();
    return M;
  }

  /// @brief The 95th percentile of `|v|`'s entries (nearest-rank, matching the
  /// diagnostic percentiles reported by `EphemerisFitErrorStats`).
  inline double Percentile95(VecXd v) {
    std::sort(v.data(), v.data() + v.size());
    int idx = std::max(0, static_cast<int>(std::ceil(0.95 * v.size())) - 1);
    idx = std::min<int>(idx, static_cast<int>(v.size()) - 1);
    return v(idx);
  }

  /// @brief RTN-decomposed RMS/95th-percentile position and velocity fit-error
  /// statistics between a fitted trajectory `rv_fit` and a reference trajectory
  /// `rv_ref` (both [N x 6], same epochs), used by `CartesianEphemeris::EvalError`
  /// and `Almanac::EvalError`.
  inline EphemerisFitErrorStats ComputeFitErrorStats(const MatXd& rv_fit, const MatXd& rv_ref) {
    const int n = static_cast<int>(rv_ref.rows());
    MatXd diff = rv_fit - rv_ref;
    MatXd diff_rtn(n, 6);
    for (int i = 0; i < n; ++i) {
      const Vec6d rv_row = rv_ref.row(i).transpose();
      const Mat3d M = RtnMatrix(rv_row);
      diff_rtn.row(i).head<3>() = (M * diff.row(i).head<3>().transpose()).transpose();
      diff_rtn.row(i).tail<3>() = (M * diff.row(i).tail<3>().transpose()).transpose();
    }

    EphemerisFitErrorStats stats;
    for (int d = 0; d < 3; ++d) {
      stats.rms_pos_m(d) = std::sqrt(diff_rtn.col(d).array().square().mean());
      stats.rms_vel_mps(d) = std::sqrt(diff_rtn.col(3 + d).array().square().mean());
      stats.p95_pos_m(d) = Percentile95(diff_rtn.col(d).array().abs().matrix());
      stats.p95_vel_mps(d) = Percentile95(diff_rtn.col(3 + d).array().abs().matrix());
    }
    stats.rms_pos_m(3) = std::sqrt(diff_rtn.leftCols(3).array().square().sum() / n);
    stats.rms_vel_mps(3) = std::sqrt(diff_rtn.rightCols(3).array().square().sum() / n);

    VecXd pos_norm(n), vel_norm(n);
    for (int i = 0; i < n; ++i) {
      pos_norm(i) = diff_rtn.row(i).head<3>().norm();
      vel_norm(i) = diff_rtn.row(i).tail<3>().norm();
    }
    stats.p95_pos_m(3) = Percentile95(pos_norm);
    stats.p95_vel_mps(3) = Percentile95(vel_norm);

    return stats;
  }

}  // namespace lupnt

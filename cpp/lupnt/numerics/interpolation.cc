#include "lupnt/numerics/interpolation.h"

#include <string>

namespace lupnt {

  /// @brief Linear interpolation in 1D
  /// @param x Vector of x values
  /// @param z Vector of z values
  /// @param xi x interpolation point
  double LinearInterp1d(const VecXd& x, const VecXd& z, double xi) {
    if (x.size() != z.size()) throw std::runtime_error("Invalid size");
    if (xi < x(0) - EPS || xi > x(x.size() - 1) + EPS) throw std::runtime_error("Out of range");

    int i = 0;
    for (i = 0; i < x.size() - 1; ++i)
      if (xi <= x(i + 1) + EPS) break;

    double dx0 = (xi - x(i)) / (x(i + 1) - x(i));
    double dx1 = (x(i + 1) - xi) / (x(i + 1) - x(i));
    double result = z(i) * dx1 + z(i + 1) * dx0;
    return result;
  }

  /// @brief Linear interpolation in 2D
  /// @param x Vector of x values
  /// @param y Vector of y values
  /// @param z Matrix of z values
  /// @param xi x interpolation point
  double LinearInterp2d(const VecXd& x, const VecXd& y, const MatXd& z, double xi, double yi) {
    if (x.size() != z.rows() || y.size() != z.cols()) throw std::runtime_error("Invalid size");
    if (xi < x(0) - EPS || xi > x(x.size() - 1) + EPS || yi < y(0) - EPS
        || yi > y(y.size() - 1) + EPS)
      throw std::runtime_error("Out of range");

    int i = 0, j = 0;
    for (i = 0; i < x.size() - 1; ++i)
      if (xi <= x(i + 1) + EPS) break;
    for (j = 0; j < y.size() - 1; ++j)
      if (yi <= y(j + 1) + EPS) break;

    double dx0 = (xi - x(i)) / (x(i + 1) - x(i));
    double dx1 = (x(i + 1) - xi) / (x(i + 1) - x(i));
    double dy0 = (yi - y(j)) / (y(j + 1) - y(j));
    double dy1 = (y(j + 1) - yi) / (y(j + 1) - y(j));
    double result = z(i, j) * dx1 * dy1 + z(i, j + 1) * dx1 * dy0 + z(i + 1, j) * dx0 * dy1
                    + z(i + 1, j + 1) * dx0 * dy0;
    return result;
  }

  LagrangeInterpolator::LagrangeInterpolator(const VecXd& x, double xi, int order)
      : x_(x), xi_(xi), order_(order) {
    if (x.size() <= order) throw std::runtime_error("Invalid size");
    if (xi < x(0) - EPS || xi > x(x.size() - 1) + EPS) throw std::runtime_error("Out of range");
    ComputeFirstIndex();
    ComputeWeights();
  }

  void LagrangeInterpolator::ComputeFirstIndex() {
    // Minimize distance to mid point: |(x(i0) + x(i0+order))/ 2 - xi|
    i0_ = 0;
    double d_best = std::abs(x_(i0_) + x_(i0_ + order_) - 2 * xi_);
    for (int i = i0_ + 1; i < x_.size() - order_; ++i) {
      double d = std::abs(x_(i) + x_(i + order_) - 2 * xi_);
      if (d < d_best) {
        i0_ = i;
        d_best = d;
      } else {
        break;
      }
    }
  }

  void LagrangeInterpolator::ComputeWeights() {
    weights_ = VecXd::Zero(order_);
    for (int i = 0; i < order_; ++i) {
      weights_(i) = 1.0;
      for (int j = 0; j < order_; ++j) {
        if (j == i) continue;
        weights_(i) *= (xi_ - x_(i0_ + j)) / (x_(i0_ + i) - x_(i0_ + j));
      }
    }
  }

  double LagrangeInterpolator::Interpolate(const VecXd& z) {
    if (z.size() != x_.size()) throw std::runtime_error("Invalid input size");
    double result = 0;
    for (int i = 0; i < order_; ++i) {
      result += weights_(i) * z(i0_ + i);
    }
    return result;
  }

}  // namespace lupnt

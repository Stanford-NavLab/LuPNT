#include "lupnt/numerics/interpolation.h"
namespace lupnt {

  /// @brief Linear interpolation in 1D
  /// @param x Vector of x values
  /// @param data Vector of data values
  double LinearInterp1d(const VecXd& x, const VecXd& data, double xi) {
    assert(x.size() == data.size());
    assert(x.size() > 1);
    assert(xi >= x[0] && xi <= x[x.size() - 1]);

    const double* start = x.data();
    const double* end = x.data() + x.size();
    auto it = std::lower_bound(start, end, xi);
    size_t i0 = it - start;
    size_t i1 = i0 + 1;

    double dx0 = xi - x[i0];
    double dx1 = x[i1] - xi;

    dx0 = dx0 / (x[i1] - x[i0]);
    dx1 = dx1 / (x[i1] - x[i0]);

    double result = data[i0] * dx1 + data[i1] * dx0;
    return result;
  }

  /// @brief Linear interpolation in 2D
  /// @param x Vector of x values
  /// @param y Vector of y values
  /// @param data Matrix of data values
  /// @param xi x interpolation point
  double LinearInterp2d(const VecXd& x, const VecXd& y, const MatXd& data, double xi, double yi) {
    assert(x.size() * y.size() == data.size());
    assert(x.size() > 1);
    assert(y.size() > 1);
    assert(xi >= x[0] && xi <= x[x.size() - 1]);
    assert(yi >= y[0] && yi <= y[y.size() - 1]);

    const double *start, *end;

    start = x.data();
    end = x.data() + x.size();
    auto it = std::lower_bound(start, end, xi);
    size_t i0 = it - start - 1;
    size_t i1 = i0 + 1;

    start = y.data();
    end = y.data() + y.size();
    it = std::lower_bound(start, end, yi);
    size_t j0 = it - start - 1;
    size_t j1 = j0 + 1;

    double dx0 = (xi - x[i0]) / (x[i1] - x[i0]);
    double dx1 = (x[i1] - xi) / (x[i1] - x[i0]);
    double dy0 = (yi - y[j0]) / (y[j1] - y[j0]);
    double dy1 = (y[j1] - yi) / (y[j1] - y[j0]);
    double result = data(i0, j0) * dx1 * dy1 + data(i0, j1) * dx1 * dy0 + data(i1, j0) * dx0 * dy1
                    + data(i1, j1) * dx0 * dy0;
    return result;
  }

  LagrangeInterpolator::LagrangeInterpolator(const VecXd& x, double xi, int order)
      : x_(x), xi_(xi), order_(order) {
    assert(x.size() > order);
    assert(xi >= x[0] && xi <= x[x.size() - 1]);
    ComputeFirstIndex();
    ComputeWeights();
  }

  void LagrangeInterpolator::ComputeFirstIndex() {
    const double* start = x_.data();
    const double* end = x_.data() + x_.size();
    auto it = std::lower_bound(start, end, xi_);
    i0_ = (it - start) - order_ / 2;
    if (i0_ < 0) i0_ = 0;

    // minimize abs(x(i0) + x(i0+order) - 2*xi)
    double min_diff = std::abs(x_[i0_] + x_[i0_ + order_] - 2 * xi_);

    double diff;
    // Decrease i0 until the difference is increasing
    for (int i = i0_ - 1; i >= 0; --i) {
      diff = std::abs(x_[i] + x_[i + order_] - 2 * xi_);
      if (diff > min_diff) break;
      min_diff = diff;
      i0_ = i;
    }
    for (int i = i0_ + 1; i < x_.size() - order_; ++i) {
      diff = std::abs(x_[i] + x_[i + order_] - 2 * xi_);
      if (diff > min_diff) break;
      min_diff = diff;
      i0_ = i;
    }

    assert(i0_ >= 0);
    assert(i0_ + order_ < x_.size());
  }

  void LagrangeInterpolator::ComputeWeights() {
    weights_ = VecXd::Zero(order_);
    for (int i = 0; i < order_; ++i) {
      weights_(i) = 1.0;
      for (int j = 0; j < order_; ++j) {
        if (j == i) continue;
        weights_(i) *= (xi_ - x_[i0_ + j]) / (x_[i0_ + i] - x_[i0_ + j]);
      }
    }
  }

  double LagrangeInterpolator::Interpolate(const VecXd& data) {
    assert(data.size() == x_.size());
    double result = 0;
    for (int i = 0; i < order_; ++i) {
      result += weights_(i) * data(i0_ + i);
    }
    return result;
  }

}  // namespace lupnt

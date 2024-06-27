/**
 * @file MathUtils.cpp
 * @author Stanford NAV LAB
 * @brief Math util functions
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lupnt/numerics/math_utils.h"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <autodiff/forward/real.hpp>

namespace lupnt {

Real angleBetweenVecs(const VecX &a, const VecX &b) {
  assert(a.size() == b.size());
  return 2.0 * atan2((a.normalized() - b.normalized()).norm(),
                     (a.normalized() + b.normalized()).norm());
}

Real Wrap2Pi(Real angle) { return atan2(sin(angle), cos(angle)); }

VecX Wrap2Pi(VecX angle) {
  VecX result = angle;
  for (int i = 0; i < angle.size(); i++) {
    result(i) = Wrap2Pi(angle(i));
  }
  return result;
}

/**
 * @brief wrap the angle between 0 and 2pi
 *
 * @param angle   angle in radians
 * @return real  wrapped angle in radians
 */
Real Wrap2TwoPi(Real angle) { return angle - TWO_PI * floor(angle / TWO_PI); }

/**
 * @brief  Wrap the angles between 0 and 2pi
 *
 * @param angle  angle vector in radians
 * @return VecX  wrapped angle vector in radians
 */
VecX WrapToTwoPi(VecX angle) {
  VecX result = angle;
  for (int i = 0; i < angle.size(); i++) {
    result(i) = Wrap2TwoPi(angle(i));
  }
  return result;
}

Real round(Real x, int n) {
  Real y = x;
  y[0] = std::round(x.val() * std::pow(10, n)) / std::pow(10, n);
  return y;
}

Real deg2rad(Real deg) { return (M_PI / 180) * deg; }

Real rad2deg(Real rad) { return (180 / M_PI) * rad; }

double deg2rad(double deg) { return (M_PI / 180) * deg; }

double rad2deg(double rad) { return (180 / M_PI) * rad; }

Real decimal2dB(Real x) { return 10 * log10(x); }

Real dB2decimal(Real x) { return pow(10, x / 10); }

MatX decimal2dB(MatX x) {
  x.array() = x.unaryExpr([](Real x) { return decimal2dB(x); });
  return x;
}

MatX dB2decimal(MatX x) {
  x.array() = x.unaryExpr([](Real x) { return dB2decimal(x); });
  return x;
}

Real floor(Real x) {
  Real y = x;
  y[0] = std::floor(x.val());
  return y;
}

Real frac(Real x) {
  Real y = x;
  y[0] = x.val() - std::floor(x.val());
  return y;
}

Vec3 degrees2dms(Real deg) {
  Real d = floor(deg);
  Real m = floor((deg - d) * 60);
  Real s = (deg - d - m / 60) * 3600;
  return Vec3{d, m, s};
}

Real dms2degrees(Vec3 hms) {
  Real decdeg = hms(0) + hms(1) / 60.0 + hms(2) / 3600.0;
  return decdeg;
}

Real safe_acos(Real x) {
  if (x >= 1.0) {
    return acos(x - 1e-15);
  } else if (x <= -1.0) {
    return acos(x + 1e-15);
  } else {
    return acos(x);
  }
}

Real safe_asin(Real x) {
  if (x >= 1.0) {
    return asin(x - 1e-16);
  } else if (x <= -1.0) {
    return asin(x + 1e-16);
  } else {
    return asin(x);
  }
}

double Rms(VecXd vec) {
  return std::sqrt(vec.array().pow(2).sum() / vec.size());
}

double Percentile(VecXd vec, double p) {
  std::sort(vec.data(), vec.data() + vec.size());
  int index = std::ceil(p * vec.size());
  if (index > (vec.size() - 1)) {
    index = vec.size() - 1;
  }
  return vec(index);
}

double Std(VecXd vec) {
  double mean = vec.sum() / vec.size();
  double sq_sum = 0.0;
  for (int i = 0; i < vec.size(); i++) {
    sq_sum += (vec(i) - mean) * (vec(i) - mean);
  }
  return std::sqrt(sq_sum / vec.size());
}

double LinearInterp1d(VecXd x, VecXd data, double ix) {
  int ix0 = 0;
  int ix1 = 0;
  for (int i = 0; i < x.size() - 1; i++) {
    if (x[i] <= ix && x[i + 1] >= ix) {
      ix0 = i;
      ix1 = i + 1;
      break;
    }
  }

  // Compute the distances from the interpolation point to the nearest sample
  // points
  double dx0 = ix - x[ix0];
  double dx1 = x[ix1] - ix;

  dx0 = dx0 / (x[ix1] - x[ix0]);
  dx1 = dx1 / (x[ix1] - x[ix0]);

  // Compute the interpolated value using the four nearest data points
  double result = data[ix0] * dx1 + data[ix1] * dx0;
  return result;
}

double LinearInterp2d(VecXd x, VecXd y, VecXd data, double ix, double iy) {
  int ix0 = 0;
  int ix1 = 0;
  for (int i = 0; i < x.size() - 1; i++) {
    if (x[i] <= ix && x[i + 1] >= ix) {
      ix0 = i;
      ix1 = i + 1;
      break;
    }
  }
  int iy0 = 0;
  int iy1 = 0;
  for (int i = 0; i < y.size() - 1; i++) {
    if (y[i] <= iy && y[i + 1] >= iy) {
      iy0 = i;
      iy1 = i + 1;
      break;
    }
  }

  // Compute the distances from the interpolation point to the nearest sample
  // points
  double dx0 = ix - x[ix0];
  double dx1 = x[ix1] - ix;
  double dy0 = iy - y[iy0];
  double dy1 = y[iy1] - iy;

  dx0 = dx0 / (x[ix1] - x[ix0]);
  dx1 = dx1 / (x[ix1] - x[ix0]);
  dy0 = dy0 / (y[iy1] - y[iy0]);
  dy1 = dy1 / (y[iy1] - y[iy0]);

  // Compute the interpolated value using the four nearest data points
  double result = data(ix0, iy0) * dx1 * dy1 + data(ix0, iy1) * dx1 * dy0 +
                  data(ix1, iy0) * dx0 * dy1 + data(ix1, iy1) * dx0 * dy0;
  return result;
}

VecXd SampleMVN(const VecXd mean, const VecXd covar, int nn, int seed) {
  // Define random generator with Gaussian distribution
  int xsize = mean.size();
  auto generator = std::mt19937(seed);
  auto dist = std::bind(std::normal_distribution<double>{0.0, 1.0}, generator);

  // Transform Mat
  VecXd normTransform(xsize, xsize);
  Eigen::LLT<VecXd> cholSolver(covar);

  if (cholSolver.info() == Eigen::Success) {
    // Use cholesky solver
    normTransform = cholSolver.matrixL();
  } else {
    std::runtime_error(
        "The covariance matrix must be symmetric and pos-definite.");
  }

  VecXd randN(xsize, nn);
  VecXd mean_samples(xsize, nn);
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < nn; j++) {
      randN(i, j) = dist();
      mean_samples(i, j) = mean(i);
    }
  }

  VecXd samples = normTransform * randN + mean_samples;
  return samples;
};

MatX SampleMVN(const VecX mean, const MatX covar, int nn, int seed) {
  // Define random generator with Gaussian distribution
  int xsize = mean.size();
  auto generator = std::mt19937(seed);
  auto dist = std::bind(std::normal_distribution<double>{0.0, 1.0}, generator);

  // Transform Mat
  MatX normTransform(xsize, xsize);
  Eigen::LLT<MatX> cholSolver(covar);

  if (cholSolver.info() == Eigen::Success) {
    // Use cholesky solver
    normTransform = cholSolver.matrixL();
  } else {
    std::runtime_error(
        "The covariance matrix must be symmetric and pos-definite.");
  }

  MatX randN(xsize, nn);
  MatX mean_samples(xsize, nn);
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < nn; j++) {
      randN(i, j) = dist();
      mean_samples(i, j) = mean(i);
    }
  }

  MatX samples = normTransform * randN + mean_samples;
  return samples;
};

VecXd blkdiag(const VecXd &A, const VecXd &B) {
  VecXd C = VecXd::Zero(A.rows() + B.rows(), A.cols() + B.cols());
  C.topLeftCorner(A.rows(), A.cols()) = A;
  C.bottomRightCorner(B.rows(), B.cols()) = B;
  return C;
}

Mat3 RotX(Real angle) {
  Real c = cos(angle);
  Real s = sin(angle);
  Mat3 R{
      {1.0, 0.0, 0.0},
      {0.0, c, s},
      {0.0, -s, c},
  };
  return R;
}

Mat3 RotY(Real angle) {
  Real c = cos(angle);
  Real s = sin(angle);
  Mat3 R{
      {c, 0.0, -s},
      {0.0, 1.0, 0.0},
      {s, 0.0, c},
  };
  return R;
}

Mat3 RotZ(Real angle) {
  Real c = cos(angle);
  Real s = sin(angle);
  Mat3 R{
      {c, s, 0.0},
      {-s, c, 0.0},
      {0.0, 0.0, 1.0},
  };
  return R;
}

Mat3 Skew(Vec3 x) {
  Mat3 skew{
      {0.0, -x(2), x(1)},
      {x(2), 0.0, -x(0)},
      {-x(1), x(0), 0.0},
  };
  return skew;
}

std::vector<double> Eigen2StdVec(const VecX &vec) {
  std::vector<double> result(vec.size());
  for (int i = 0; i < vec.size(); i++) {
    result[i] = vec(i).val();
  }
  return result;
}
}  // namespace lupnt

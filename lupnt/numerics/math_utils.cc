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

real angleBetweenVectors(const VectorX &a, const VectorX &b) {
  assert(a.size() == b.size());
  return 2.0 * atan2((a.normalized() - b.normalized()).norm(),
                     (a.normalized() + b.normalized()).norm());
}

real wrapToPi(real angle) { return atan2(sin(angle), cos(angle)); }

VectorX wrapToPi(VectorX angle) {
  VectorX result = angle;
  for (int i = 0; i < angle.size(); i++) {
    result(i) = wrapToPi(angle(i));
  }
  return result;
}

/**
 * @brief wrap the angle between 0 and 2pi
 *
 * @param angle   angle in radians
 * @return real  wrapped angle in radians
 */
real wrapTo2Pi(real angle) { return atan2(sin(angle), cos(angle)) + M_PI; }

/**
 * @brief  Wrap the angles between 0 and 2pi
 *
 * @param angle  angle vector in radians
 * @return VectorX  wrapped angle vector in radians
 */
VectorX wrapTo2Pi(VectorX angle) {
  VectorX result = angle;
  for (int i = 0; i < angle.size(); i++) {
    result(i) = wrapTo2Pi(angle(i));
  }
  return result;
}

real deg2rad(real deg) { return (M_PI / 180) * deg; }

real rad2deg(real rad) { return (180 / M_PI) * rad; }

double deg2rad(double deg) { return (M_PI / 180) * deg; }

double rad2deg(double rad) { return (180 / M_PI) * rad; }

real decimal2dB(real x) { return 10 * log10(x); }

real dB2decimal(real x) { return pow(10, x / 10); }

MatrixX decimal2dB(MatrixX x) {
  x.array() = x.unaryExpr([](real x) { return decimal2dB(x); });
  return x;
}

MatrixX dB2decimal(MatrixX x) {
  x.array() = x.unaryExpr([](real x) { return dB2decimal(x); });
  return x;
}

real floor(real x) {
  real y = x;
  y[0] = std::floor(x.val());
  return y;
}

Vector3 degrees2dms(real deg) {
  real d = floor(deg);
  real m = floor((deg - d) * 60);
  real s = (deg - d - m / 60) * 3600;
  return Vector3{d, m, s};
}

real dms2degrees(Vector3 hms) {
  real decdeg = hms(0) + hms(1) / 60.0 + hms(2) / 3600.0;
  return decdeg;
}

real safe_acos(real x) {
  if (x >= 1.0) {
    return acos(x - 1e-15);
  } else if (x <= -1.0) {
    return acos(x + 1e-15);
  } else {
    return acos(x);
  }
}

real safe_asin(real x) {
  if (x >= 1.0) {
    return asin(x - 1e-16);
  } else if (x <= -1.0) {
    return asin(x + 1e-16);
  } else {
    return asin(x);
  }
}

double Rms(VectorXd vec) {
  return std::sqrt(vec.array().pow(2).sum() / vec.size());
}

double Percentile(VectorXd vec, double p) {
  std::sort(vec.data(), vec.data() + vec.size());
  int index = std::ceil(p * vec.size());
  if (index > (vec.size() - 1)) {
    index = vec.size() - 1;
  }
  return vec(index);
}

double Std(VectorXd vec) {
  double mean = vec.sum() / vec.size();
  double sq_sum = 0.0;
  for (int i = 0; i < vec.size(); i++) {
    sq_sum += (vec(i) - mean) * (vec(i) - mean);
  }
  return std::sqrt(sq_sum / vec.size());
}

double LinearInterp1d(VectorXd x, VectorXd data, double ix) {
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

double LinearInterp2d(VectorXd x, VectorXd y, MatrixXd data, double ix,
                      double iy) {
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

MatrixXd SampleMVN(const VectorXd mean, const MatrixXd covar, int nn,
                   int seed) {
  // Define random generator with Gaussian distribution
  int xsize = mean.size();
  auto generator = std::mt19937(seed);
  auto dist = std::bind(std::normal_distribution<double>{0.0, 1.0}, generator);

  // Transform Matrix
  MatrixXd normTransform(xsize, xsize);
  Eigen::LLT<MatrixXd> cholSolver(covar);

  if (cholSolver.info() == Eigen::Success) {
    // Use cholesky solver
    normTransform = cholSolver.matrixL();
  } else {
    std::runtime_error(
        "The covariance matrix must be symmetric and pos-definite.");
  }

  MatrixXd randN(xsize, nn);
  MatrixXd mean_samples(xsize, nn);
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < nn; j++) {
      randN(i, j) = dist();
      mean_samples(i, j) = mean(i);
    }
  }

  MatrixXd samples = normTransform * randN + mean_samples;
  return samples;
};

MatrixX SampleMVN(const VectorX mean, const MatrixX covar, int nn, int seed) {
  // Define random generator with Gaussian distribution
  int xsize = mean.size();
  auto generator = std::mt19937(seed);
  auto dist = std::bind(std::normal_distribution<double>{0.0, 1.0}, generator);

  // Transform Matrix
  MatrixX normTransform(xsize, xsize);
  Eigen::LLT<MatrixX> cholSolver(covar);

  if (cholSolver.info() == Eigen::Success) {
    // Use cholesky solver
    normTransform = cholSolver.matrixL();
  } else {
    std::runtime_error(
        "The covariance matrix must be symmetric and pos-definite.");
  }

  MatrixX randN(xsize, nn);
  MatrixX mean_samples(xsize, nn);
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < nn; j++) {
      randN(i, j) = dist();
      mean_samples(i, j) = mean(i);
    }
  }

  MatrixX samples = normTransform * randN + mean_samples;
  return samples;
};

MatrixXd blkdiag(const MatrixXd &A, const MatrixXd &B) {
  MatrixXd C = MatrixXd::Zero(A.rows() + B.rows(), A.cols() + B.cols());
  C.topLeftCorner(A.rows(), A.cols()) = A;
  C.bottomRightCorner(B.rows(), B.cols()) = B;
  return C;
}

Matrix3 Rot1(real phi) {
  real c = cos(phi);
  real s = sin(phi);
  Matrix3 R1{
      {1.0, 0.0, 0.0},
      {0.0, c, s},
      {0.0, -s, c},
  };
  return R1;
}

Matrix3 Rot2(real phi) {
  real c = cos(phi);
  real s = sin(phi);
  Matrix3 R2{
      {c, 0.0, -s},
      {0.0, 1.0, 0.0},
      {s, 0.0, c},
  };
  return R2;
}

Matrix3 Rot3(real phi) {
  real c = cos(phi);
  real s = sin(phi);
  Matrix3 R3{
      {c, s, 0.0},
      {-s, c, 0.0},
      {0.0, 0.0, 1.0},
  };
  return R3;
}

Matrix3 Skew(Vector3 x) {
  Matrix3 skew{
      {0.0, -x(2), x(1)},
      {x(2), 0.0, -x(0)},
      {-x(1), x(0), 0.0},
  };
  return skew;
}

std::vector<double> EigenToStdVector(const VectorX &vec) {
  std::vector<double> result(vec.size());
  for (int i = 0; i < vec.size(); i++) {
    result[i] = vec(i).val();
  }
  return result;
}
}  // namespace lupnt

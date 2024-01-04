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

#include "lupnt/core/constants.h"

namespace lupnt {

std::tuple<real, real, real> unpack(const Vector3 &vec) {
  return std::make_tuple(vec(0), vec(1), vec(2));
}
std::tuple<real, real, real, real, real, real> unpack(const Vector6 &vec) {
  return std::make_tuple(vec(0), vec(1), vec(2), vec(3), vec(4), vec(5));
}

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

real wrapTo2Pi(real angle) { return atan2(sin(angle), cos(angle)) + M_PI; }
VectorX wrapTo2Pi(VectorX angle) {
  VectorX result = angle;
  for (int i = 0; i < angle.size(); i++) {
    result(i) = wrapTo2Pi(angle(i));
  }
  return result;
}

/**
 * @brief Convert degree to radian
 *
 * @param deg
 * @return real
 */
real deg2rad(real deg) { return (M_PI / 180) * deg; }

/**
 * @brief Convert radian to degree
 *
 * @param rad
 * @return real
 */
real rad2deg(real rad) { return (180 / M_PI) * rad; }

double deg2rad(double deg) { return (M_PI / 180) * deg; }

double rad2deg(double rad) { return (180 / M_PI) * rad; }

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

/**
 * @brief Root Mean Square
 *
 * @param vec  target vector
 * @return double  RMS value
 */
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

/**
 * @brief Linear interpolation of a 2D data set
 *
 * @param x  data x axis
 * @param y  data
 * @param ix  interpolation point x
 * @return double  data at the interpolation point
 */
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

/**
 * @brief Linear interpolation of a 2D data set
 *
 * @param x  data x axis
 * @param y  data y axis
 * @param data data matrix
 * @param ix  interpolation point x
 * @param iy  interpolation point y
 * @return double  data at the interpolation point
 */
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

/**
 * @brief Sample from a multivariate normal distribution
 *
 * @param mean
 * @param cov
 * @param nn
 * @return MatrixXd
 */
MatrixXd SampleMVN(const VectorXd mean, const MatrixXd covar, int nn) {
  // Define random generator with Gaussian distribution
  int xsize = mean.size();
  auto dist = std::bind(std::normal_distribution<double>{0.0, 1.0},
                        std::mt19937(std::random_device{}()));

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

MatrixX SampleMVN(const VectorX mean, const MatrixX covar, int nn) {
  // Define random generator with Gaussian distribution
  int xsize = mean.size();
  auto dist = std::bind(std::normal_distribution<double>{0.0, 1.0},
                        std::mt19937(std::random_device{}()));

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

/**
 * @brief Compute the rotation matrix specified by the input skew vector.
 * @param v
 * @return Matrix3
 */
Matrix3 Skew(Vector3 x) {
  Matrix3 skew{
      {0.0, -x(2), x(1)},
      {x(2), 0.0, -x(0)},
      {-x(1), x(0), 0.0},
  };
  return skew;
}

}  // namespace lupnt

/**
 * @file MathUtils.h
 * @author Keidai Iiyama
 * @brief Utility Mathematical Functions
 * @version 0.1
 * @date 2023-02-07
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lupnt/numerics/MathUtils.h"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

namespace ad = autodiff;

namespace LPT {

Eigen::VectorXd toEigen(ad::VectorXreal x) {
  Eigen::VectorXd v(x.size());
  for (int i = 0; i < x.size(); i++) {
    v[i] = x[i].val();
  }
  return v;
}

ad::real dot(ad::VectorXreal v1, ad::VectorXreal v2) {
  return (v1.transpose() * v2).sum();
}

ad::VectorXreal cross(ad::VectorXreal v1, ad::VectorXreal v2) {
  int n = v1.size();
  ad::VectorXreal vout(3);

  if (n == 3)
    vout << v1(1) * v2(2) - v1(2) * v2(1), v1(2) * v2(0) - v1(0) * v2(2),
        v1(0) * v2(1) - v1(1) * v2(0);
  else
    throw std::invalid_argument(
        "Currently the cross product is only implemetef for n=3");
  return vout;
}

ad::real norm(ad::VectorXreal x) { return sqrt(dot(x, x)); }

/**
 * @brief Wrap the angle between -pi and pi
 *
 * @param angle
 * @return ad::real
 */
ad::real wrapToPi(ad::real angle) { return atan2(sin(angle), cos(angle)); }

/**
 * @brief Wrap the angle between 0 and 2pi
 *
 * @param angle
 * @return ad::real
 */
ad::real wrapTo2Pi(ad::real angle) {
  return atan2(sin(angle), cos(angle)) + M_PI;
}

/**
 * @brief Convert degree to radian
 *
 * @param deg
 * @return ad::real
 */
ad::real degToRad(ad::real deg) { return (M_PI / 180) * deg; }

/**
 * @brief Convert radian to degree
 *
 * @param rad
 * @return ad::real
 */
ad::real radToDeg(ad::real rad) { return (180 / M_PI) * rad; }

double degToRad(double deg) { return (M_PI / 180) * deg; }

double radToDeg(double rad) { return (180 / M_PI) * rad; }

/**
 * @brief Linear interpolation of a 2D data set
 *
 * @param x  data x axis
 * @param y  data
 * @param ix  interpolation point x
 * @return double  data at the interpolation point
 */
double LinearInterp1d(Eigen::VectorXd x, Eigen::VectorXd data, double ix) {
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
double LinearInterp2d(Eigen::VectorXd x, Eigen::VectorXd y,
                      Eigen::MatrixXd data, double ix, double iy) {
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
 * @return Eigen::MatrixXd
 */
Eigen::MatrixXd SampleMVN(const Eigen::VectorXd mean,
                          const Eigen::MatrixXd covar, int nn) {
  // Define random generator with Gaussian distribution
  int xsize = mean.size();
  auto dist = std::bind(std::normal_distribution<double>{0.0, 1.0},
                        std::mt19937(std::random_device{}()));

  // Transform Matrix
  Eigen::MatrixXd normTransform(xsize, xsize);
  Eigen::LLT<Eigen::MatrixXd> cholSolver(covar);

  if (cholSolver.info() == Eigen::Success) {
    // Use cholesky solver
    normTransform = cholSolver.matrixL();
  } else {
    std::runtime_error(
        "The covariance matrix must be symmetric and pos-definite.");
  }

  Eigen::MatrixXd randN(xsize, nn);
  Eigen::MatrixXd mean_samples(xsize, nn);
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < xsize; j++) {
      randN(i, j) = dist();
      mean_samples(i, j) = mean(i);
    }
  }

  Eigen::MatrixXd samples = normTransform * randN + mean_samples;
  return samples;
};

}  // namespace LPT

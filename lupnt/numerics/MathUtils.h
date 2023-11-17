/**
 * @file MathUtils.h
 * @author Stanford NAV LAB
 * @brief
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <random>

namespace ad = autodiff;

namespace LPT {

// Convert autodiff vector to Eigen vector
Eigen::VectorXd toEigen(ad::VectorXreal x);

/**
 * @brief Calculate the dot product of two ad vectors
 *
 * @param v1
 * @param v2
 * @return ad::real
 */
ad::real dot(ad::VectorXreal v1, ad::VectorXreal v2);

/**
 * @brief Calculate the cross product of two ad vectors
 *
 * @param v1
 * @param v2
 * @return ad::VectorXreal
 */
ad::VectorXreal cross(ad::VectorXreal v1, ad::VectorXreal v2);

/**
 * @brief Calculate the 2-norm of ad::vector
 *
 * @param x
 * @return ad::real
 */
ad::real norm(ad::VectorXreal x);

/**
 * @brief Wrap the angle between -pi and pi
 *
 * @param angle
 * @return ad::real
 */
template <typename T>
T wrapToPi(T angle);

/**
 * @brief Convert degree to radian
 *
 * @param deg
 * @return ad::real
 */
ad::real degToRad(ad::real deg);

/**
 * @brief Convert radian to degree
 *
 * @param rad
 * @return ad::real
 */
ad::real radToDeg(ad::real rad);

double LinearInterp1d(Eigen::VectorXd x, Eigen::VectorXd data, double ix);
double LinearInterp2d(Eigen::VectorXd x, Eigen::VectorXd y,
                      Eigen::MatrixXd data, double ix, double iy);

/**
 * @brief Sample from a multivariate normal distribution
 *
 */
Eigen::MatrixXd SampleMVN(const Eigen::VectorXd mean, const Eigen::MatrixXd cov,
                          int nn);

}  // namespace LPT

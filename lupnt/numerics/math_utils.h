/**
 * @file math_utils.h
 * @author Stanford NAV LAB
 * @brief
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <lupnt/core/constants.h>
#include <random>



namespace lupnt {

// Convert autodiff vector to Eigen vector
VectorXd toEigen(VectorXreal x);

/**
 * @brief Wrap the angle between -pi and pi
 *
 * @param angle
 * @return real
 */
real wrapToPi(real angle);

/**
 * @brief Convert degree to radian
 *
 * @param deg
 * @return real
 */
real degToRad(real deg);

/**
 * @brief Convert radian to degree
 *
 * @param rad
 * @return real
 */
real radToDeg(real rad);

double LinearInterp1d(VectorXd x, VectorXd data, double ix);
double LinearInterp2d(VectorXd x, VectorXd y, MatrixXd data, double ix,
                      double iy);

/**
 * @brief Sample from a multivariate normal distribution
 *
 */
MatrixXd SampleMVN(const VectorXd mean, const MatrixXd cov, int nn);

}  // namespace lupnt

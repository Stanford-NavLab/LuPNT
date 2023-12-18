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

std::tuple<real, real, real> unpack(const Vector3 &vec);
std::tuple<real, real, real, real, real, real> unpack(const Vector6 &vec);

/**
 * @brief Wrap the angle between -pi and pi
 *
 * @param angle
 * @return real
 */
real wrapToPi(real angle);
VectorX wrapToPi(VectorX angle);

real wrapTo2Pi(real angle);
VectorX wrapTo2Pi(VectorX angle);

/**
 * @brief Convert degree to radian
 *
 * @param deg
 * @return real
 */
real deg2rad(real deg);

real safe_acos(real x);
real safe_asin(real x);

real angleBetweenVectors(const VectorX &a, const VectorX &b);

/**
 * @brief Convert radian to degree
 *
 * @param rad
 * @return real
 */
real rad2deg(real rad);

double LinearInterp1d(VectorXd x, VectorXd data, double ix);
double LinearInterp2d(VectorXd x, VectorXd y, MatrixXd data, double ix,
                      double iy);

/**
 * @brief Sample from a multivariate normal distribution
 *
 */
MatrixX SampleMVN(const VectorX mean, const MatrixX cov, int nn);

MatrixXd blkdiag(const MatrixXd &A, const MatrixXd &B);

Matrix3 Rot1(real phi);
Matrix3 Rot2(real phi);
Matrix3 Rot3(real phi);
Matrix3 Skew(Vector3 x);

}  // namespace lupnt

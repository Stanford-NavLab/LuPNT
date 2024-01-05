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
 * @brief Compute the angle between two vectors
 *
 * @param a  vector a
 * @param b  vector b
 * @return real   angle between a and b
 */
real angleBetweenVectors(const VectorX &a, const VectorX &b);

/**
 * @brief Wrap the angle between -pi and pi
 *
 * @param angle
 * @return real
 */
real wrapToPi(real angle);
VectorX wrapToPi(VectorX angle);

/**
 * @brief Wrap the angles between 0 and 2pi
 *
 * @param angle
 * @return real
 */
real wrapTo2Pi(real angle);
VectorX wrapTo2Pi(VectorX angle);

/**
 * @brief Convert degree to radian
 *
 * @param deg
 * @return real
 */
real deg2rad(real deg);

/**
 * @brief Arccosine function with input bounds
 *
 * @param x  input
 * @return real  output
 */
real safe_acos(real x);

/**
 * @brief Arcsine function with input bounds
 *
 * @param x  input
 * @return real  output
 */
real safe_asin(real x);

/**
 * @brief Convert radian to degree
 *
 * @param rad
 * @return real
 */
real rad2deg(real rad);

/**
 * @brief Compute the root mean square of a vector
 *
 * @param vec  input vector
 * @return double   rms value
 */
double Rms(VectorXd vec);

/**
 * @brief  Compute the pth percentile of a vector
 *
 * @param vec  input vector
 * @param p  percentile (0-1)
 * @return double  percentile value
 */
double Percentile(VectorXd vec, double p);

/**
 * @brief Standard deviation of a vector
 *
 * @param vec   input vector
 * @return double   standard deviation
 */
double Std(VectorXd vec);

/**
 * @brief Linear interpolation of a 2D data set
 *
 * @param x  data x axis
 * @param y  data
 * @param ix  interpolation point x
 * @return double  data at the interpolation point
 */
double LinearInterp1d(VectorXd x, VectorXd data, double ix);

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
                      double iy);

// /**
//  * @brief Sample from a multivariate normal distribution
//  *
//  * @param mean  mean vector
//  * @param covar  covariance matrix
//  * @param nn   number of samples to generate
//  * @return MatrixXd
//  */
// MatrixXd SampleMVN(const VectorXd mean, const MatrixXd covar, int nn);

/**
 * @brief  Sample from a multivariate normal distribution
 *
 * @param mean  mean vector
 * @param covar  covariance matrix
 * @param nn  number of samples to generate
 * @return MatrixX
 */
MatrixX SampleMVN(const VectorX mean, const MatrixX cov, int nn);

/**
 * @brief block diagonal matrix
 *
 * @param A  Matrix A
 * @param B  Matrix B
 * @return MatrixXd  block diagonal matrix of A and B
 */
MatrixXd blkdiag(const MatrixXd &A, const MatrixXd &B);

/**
 * @brief Rotate a vector by a given angle about the x-axis
 *
 * @param phi  angle in radians
 * @return Matrix3   rotation matrix
 */
Matrix3 Rot1(real phi);

/**
 * @brief  Rotate a vector by a given angle about the y-axis
 *
 * @param phi  angle in radians
 * @return Matrix3  rotation matrix
 */
Matrix3 Rot2(real phi);

/**
 * @brief  Rotate a vector by a given angle about the z-axis
 *
 * @param phi   angle in radians
 * @return Matrix3  rotation matrix
 */
Matrix3 Rot3(real phi);

/**
 * @brief Compute the rotation matrix specified by the input skew vector.
 * @param v
 * @return Matrix3
 */
Matrix3 Skew(Vector3 x);

}  // namespace lupnt

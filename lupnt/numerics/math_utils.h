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
#include <tuple>
#include <utility>

namespace lupnt {

template <typename Vec, std::size_t... Indices>
auto unpackImpl(const Vec &vec, std::index_sequence<Indices...>) {
  return std::make_tuple(vec(Indices)...);
}

template <typename T, int Size>
auto unpack(const Eigen::Matrix<T, Size, 1> &vec) {
  return unpackImpl(vec, std::make_index_sequence<Size>{});
}

/**
 * @brief Compute the angle between two vectors
 *
 * @param a  vector a
 * @param b  vector b
 * @return real   angle between a and b
 */
real angleBetweenVecs(const VecX &a, const VecX &b);

/**
 * @brief Wrap the angle between -pi and pi
 *
 * @param angle
 * @return real
 */
real wrapToPi(real angle);
VecX wrapToPi(VecX angle);

real decimal2dB(real x);
real dB2decimal(real x);

MatX decimal2dB(MatX x);
MatX dB2decimal(MatX x);

/**
 * @brief Wrap the angles between 0 and 2pi
 *
 * @param angle
 * @return real
 */
real wrapTo2Pi(real angle);
VecX wrapTo2Pi(VecX angle);

/**
 * @brief Convert degree to radian
 *
 * @param deg
 * @return real
 */
real deg2rad(real deg);
double rad2deg(double rad);

real floor(real x);
Vec3 degrees2dms(real deg);
real dms2degrees(Vec3 hms);

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
double deg2rad(double deg);

/**
 * @brief Compute the root mean square of a vector
 *
 * @param vec  input vector
 * @return double   rms value
 */
double Rms(VecXd vec);

/**
 * @brief  Compute the pth percentile of a vector
 *
 * @param vec  input vector
 * @param p  percentile (0-1)
 * @return double  percentile value
 */
double Percentile(VecXd vec, double p);

/**
 * @brief Standard deviation of a vector
 *
 * @param vec   input vector
 * @return double   standard deviation
 */
double Std(VecXd vec);

/**
 * @brief Linear interpolation of a 2D data set
 *
 * @param x  data x axis
 * @param y  data
 * @param ix  interpolation point x
 * @return double  data at the interpolation point
 */
double LinearInterp1d(VecXd x, VecXd data, double ix);

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
double LinearInterp2d(VecXd x, VecXd y, VecXd data, double ix, double iy);

// /**
//  * @brief Sample from a multivariate normal distribution
//  *
//  * @param mean  mean vector
//  * @param covar  covariance matrix
//  * @param nn   number of samples to generate
//  * @return VecXd
//  */
// VecXd SampleMVN(const VecXd mean, const VecXd covar, int nn);

/**
 * @brief  Sample from a multivariate normal distribution
 *
 * @param mean  mean vector
 * @param covar  covariance matrix
 * @param nn  number of samples to generate
 * @return MatX
 */
MatX SampleMVN(const VecX mean, const MatX cov, int nn, int seed = 0);

/**
 * @brief block diagonal matrix
 *
 * @param A  Mat A
 * @param B  Mat B
 * @return VecXd  block diagonal matrix of A and B
 */
VecXd blkdiag(const VecXd &A, const VecXd &B);

/**
 * @brief Rotate a vector by a given angle about the x-axis
 *
 * @param phi  angle in radians
 * @return Mat3   rotation matrix
 */
Mat3 Rot1(real phi);

/**
 * @brief  Rotate a vector by a given angle about the y-axis
 *
 * @param phi  angle in radians
 * @return Mat3  rotation matrix
 */
Mat3 Rot2(real phi);

/**
 * @brief  Rotate a vector by a given angle about the z-axis
 *
 * @param phi   angle in radians
 * @return Mat3  rotation matrix
 */
Mat3 Rot3(real phi);

/**
 * @brief Compute the rotation matrix specified by the input skew vector.
 * @param v
 * @return Mat3
 */
Mat3 Skew(Vec3 x);

std::vector<double> EigenToStdVec(const VecX &vec);

}  // namespace lupnt

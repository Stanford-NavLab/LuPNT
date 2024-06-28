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

Real AngleBetweenVecs(const VecX &x, const VecX &y);

Real Wrap2Pi(Real angle);
Real Wrap2TwoPi(Real angle);

Real Decimal2Decibel(Real x);
Real Decibel2Decimal(Real x);

Real round(Real x, int n = 0);
Real frac(Real x);
Real ceil(Real x);
Real floor(Real x);

Vec3 Degrees2DegMinSec(Real deg);
Real DegMinSec2Degrees(Vec3 hms);

Real safe_acos(Real x);
Real safe_asin(Real x);

Real RootMeanSquare(VecX x);
Real Percentile(VecX x, double p);
Real Std(VecX x);

double LinearInterp1d(VecXd x, VecXd data, double ix);
double LinearInterp2d(VecXd x, VecXd y, VecXd data, double ix, double iy);

MatX SampleMVN(const VecX mean, const MatX cov, int nn, int seed = 0);
MatX blkdiag(const MatX &A, const MatX &B);

Mat3 RotX(Real angle);
Mat3 RotY(Real angle);
Mat3 RotZ(Real angle);
Mat3 Skew(Vec3 x);

VecXd ToDouble(const VecX &x);
MatXd ToDouble(const MatX &x);

Real RatioOfSectorToTriangleArea(Vec3 r1, Vec3 r2, Real tau);

}  // namespace lupnt

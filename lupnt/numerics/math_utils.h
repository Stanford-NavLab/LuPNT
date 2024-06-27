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

Real angleBetweenVecs(const VecX &a, const VecX &b);

Real Wrap2Pi(Real angle);
VecX Wrap2Pi(VecX angle);

Real decimal2dB(Real x);
Real dB2decimal(Real x);

MatX decimal2dB(MatX x);
MatX dB2decimal(MatX x);

Real Wrap2TwoPi(Real angle);
VecX WrapToTwoPi(VecX angle);

Real round(Real x, int n = 0);
Real frac(Real x);

Real deg2rad(Real deg);
double rad2deg(double rad);

Real floor(Real x);
Vec3 degrees2dms(Real deg);
Real dms2degrees(Vec3 hms);

Real safe_acos(Real x);

Real safe_asin(Real x);

Real rad2deg(Real rad);
double deg2rad(double deg);

double Rms(VecXd vec);
double Percentile(VecXd vec, double p);
double Std(VecXd vec);

double LinearInterp1d(VecXd x, VecXd data, double ix);

double LinearInterp2d(VecXd x, VecXd y, VecXd data, double ix, double iy);

MatX SampleMVN(const VecX mean, const MatX cov, int nn, int seed = 0);

VecXd blkdiag(const VecXd &A, const VecXd &B);

Mat3 RotX(Real angle);
Mat3 RotY(Real angle);
Mat3 RotZ(Real angle);
Mat3 Skew(Vec3 x);

std::vector<double> Eigen2StdVec(const VecX &vec);

}  // namespace lupnt

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

#include <random>
#include <tuple>
#include <utility>

#include "lupnt/core/constants.h"

namespace lupnt {

  template <typename Vec, std::size_t... Indices>
  auto unpackImpl(const Vec& vec, std::index_sequence<Indices...>) {
    return std::make_tuple(vec(Indices)...);
  }

  template <typename T, int Size> auto unpack(const Eigen::Matrix<T, Size, 1>& vec) {
    return unpackImpl(vec, std::make_index_sequence<Size>{});
  }

  template <typename T> VectorX<T> arange(T start, T stop, T step = 1);

  Real AngleBetweenVecs(const VecX& x, const VecX& y);

  Real Wrap2Pi(Real angle);
  Real Wrap2TwoPi(Real angle);
  VecX Wrap2Pi(VecX angle);
  VecX Wrap2TwoPi(VecX angle);

  Real Decimal2Decibel(Real x);
  VecX Decimal2Decibel(VecX x);
  MatX Decimal2Decibel(MatX x);

  Real Decibel2Decimal(Real x);
  VecX Decibel2Decimal(VecX x);
  MatX Decibel2Decimal(MatX x);

  Real round(Real x, int n = 0);
  Real frac(Real x);
  Real ceil(Real x);
  Real floor(Real x);
  Real mod(Real x, Real y);

  Vec3 Degrees2DegMinSec(Real deg);
  Real DegMinSec2Degrees(Vec3 hms);

  Real sind(Real x);
  Real cosd(Real x);
  Real tand(Real x);

  Real safe_acos(Real x);
  Real safe_asin(Real x);

  double J0Bessel(double x);
  double J1Bessel(double x);

  Real RootMeanSquare(VecX x);
  Real Percentile(VecX x, double p);
  Real Std(VecX x);
  double RootMeanSquareD(VecXd x);
  double PercentileD(VecXd x, double p);
  double StdD(VecXd x);

  double erfc(double x);
  double qfunc(double x);

  MatX SampleMVN(const VecX mean, const MatX cov, int nn, int seed = 0);
  double SampleRandNormal(double mean, double std, int seed = 0);

  MatX BlkDiag(const MatX& A, const MatX& B);
  MatXd BlkDiagD(const MatXd& A, const MatXd& B);

  template <typename T> Matrix<T, 3, 3> RotX(T angle);
  template <typename T> Matrix<T, 3, 3> RotY(T angle);
  template <typename T> Matrix<T, 3, 3> RotZ(T angle);
  template <typename T> Matrix<T, 3, 3> Skew(Vector<T, 3> x);

  VecXd ToDouble(const VecX& x);
  // MatXd ToDouble(const MatX& x);
  std::vector<double> ToDoubleVec(const VecX& x);
  std::vector<double> ToDoubleVec(const VecXd& x);
  std::vector<double> ToDoubleVec(const VecXi& x);

  Real RatioOfSectorToTriangleArea(Vec3 r1, Vec3 r2, Real tau);

  VecX arange(Real start, Real stop, Real step);

}  // namespace lupnt

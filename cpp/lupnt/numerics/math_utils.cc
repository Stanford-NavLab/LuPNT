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
#include <autodiff/forward/real.hpp>
#include <cmath>
#include <random>

#include "lupnt/core/constants.h"
#include "lupnt/core/random_engine.h"

namespace lupnt {
  /// @brief Create a range of values
  /// @param start Start value
  /// @param stop Stop value
  /// @param step Step value
  /// @return Vector of values
  template <typename T> VectorX<T> Arange(T start, T stop, T step) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step) {
      values.push_back(value);
    }
    // Create a proper VectorX that owns its memory (not a Map to temporary data)
    return VectorX<T>(Eigen::Map<VectorX<T>>(values.data(), values.size()));
  }
  template VectorX<int> Arange<int>(int start, int stop, int step);
  template VectorX<double> Arange<double>(double start, double stop, double step);
  template VectorX<Real> Arange<Real>(Real start, Real stop, Real step);

  /// @brief Subsample a vector
  /// @param x Input vector
  /// @param step Step size
  /// @return Subsampled vector
  VecXd Subsample(const VecX& x, int step) {
    int n = x.size();
    VecX out(n / step);
    for (int i = 0; i < n; i += step) out(i / step) = x(i);
    return out;
  }

  /// @brief Subsample a matrix
  /// @param x Input matrix
  /// @param step_row Step size for rows
  /// @param step_col Step size for columns
  /// @return Subsampled matrix
  MatXd Subsample(const MatX& x, int step_row, int step_col) {
    int n_row = x.rows();
    int n_col = x.cols();
    MatX out(n_row / step_row, n_col / step_col);
    for (int i = 0; i < n_row; i += step_row) {
      for (int j = 0; j < n_col; j += step_col) {
        out(i / step_row, j / step_col) = x(i, j);
      }
    }
    return out;
  }

  /// @brief Compute the angle between two vectors
  /// @param x First vector
  /// @param y Second vector
  /// @return Angle between the two vectors in radians
  Real AngleBetweenVecs(const VecX& x, const VecX& y) {
    return 2.0
           * atan2((x.normalized() - y.normalized()).norm(),
                   (x.normalized() + y.normalized()).norm());
  }

  /// @brief Compute the angle between two vectors
  /// @param x First vector
  /// @param y Second vector
  /// @return Angle between the two vectors in radians
  VecX AngleBetweenVecs(const MatX& x, const MatX& y) {
    return (x.rowwise().normalized().array() * y.rowwise().normalized().array())
        .rowwise()
        .sum()
        .cwiseMax(-1.0)
        .cwiseMin(1.0)
        .acos();
  }

  /// @brief Convert degrees, minutes, and seconds to radians
  /// @param degrees Degrees
  /// @param minutes Minutes
  /// @param seconds Seconds
  /// @return Angle in radians
  double DegMinSec2DeciDeg(double degrees, double minutes, double seconds) {
    double decimalDegrees = degrees + minutes / 60.0 + seconds / 3600.0;
    return decimalDegrees;
  }

  /// @brief Wrap angle in radians to [-pi, pi]
  /// @param angle Angle in radians
  /// @return Wrapped angle in radians
  Real WrapToPi(Real angle) { return atan2(sin(angle), cos(angle)); }
  VecX WrapToPi(const VecX& angle) {
    int n = angle.size();
    VecX out(n);
    for (int i = 0; i < n; i++) out(i) = WrapToPi(angle(i));
    return out;
  }

  /// @brief Wrap angle in radians to [0, 2pi]
  /// @param angle Angle in radians
  /// @return Wrapped angle in radians
  Real WrapToTwoPi(Real angle) { return angle - TWO_PI * floor(angle / TWO_PI); }
  VecX WrapToTwoPi(const VecX& angle) {
    int n = angle.size();
    VecX out(n);
    for (int i = 0; i < n; i++) out(i) = WrapToTwoPi(angle(i));
    return out;
  }

  /// @brief Round a number to n decimal places
  /// @param x Number to round
  /// @param n Number of decimal places
  /// @return Rounded number
  Real round(Real x, int n) {
    Real y = x;
    y[0] = std::round(x.val() * std::pow(10, n)) / std::pow(10, n);
    return y;
  }

  /// @brief Get the fractional part of a number
  /// @param x Number
  /// @return Fractional part of the number
  Real frac(Real x) {
    Real y = x;
    y[0] = x.val() - std::floor(x.val());
    return y;
  }

  /// @brief Get the ceiling of a number
  /// @param x Number
  /// @return Ceiling of the number
  Real ceil(Real x) {
    Real y = x;
    y[0] = std::ceil(x.val());
    return y;
  }

  /// @brief Get the floor of a number
  /// @param x Number
  /// @return Floor of the number
  Real floor(Real x) {
    Real y = x;
    y[0] = std::floor(x.val());
    return y;
  }

  /// @brief Get the modulus of two numbers
  /// @param x First number
  /// @param y Second number
  /// @return Modulus of the two numbers
  Real mod(Real x, Real y) {
    Real z = x;
    z[0] = std::fmod(x.val(), y.val());
    return z;
  }

  /// @brief Convert decimal value to decibels
  /// @param x Decimal value
  /// @return Decibel value
  Real DecimalToDecibel(Real x) { return 10 * log10(x); }
  ArrX DecimalToDecibel(const ArrX& x) { return 10 * log10(x); }

  Real Max(Real x, Real y) { return x.val() > y.val() ? x : y; }
  Real Min(Real x, Real y) { return x.val() < y.val() ? x : y; }
  double MaxD(double x, double y) { return x > y ? x : y; }
  double MinD(double x, double y) { return x < y ? x : y; }

  /// @brief Convert decibel value to decimal
  /// @param x Decibel value
  /// @return Decimal value
  Real DecibelToDecimal(Real x) { return pow(10, x / 10); }
  ArrX DecibelToDecimal(const ArrX& x) { return pow(10, x / 10); }

  /// @brief Convert degrees to degrees, minutes, and seconds
  /// @param deg Angle in degrees
  /// @return Angle in degrees, minutes, and seconds
  Vec3 DegToDegMinSec(Real deg) {
    Real d = floor(deg);
    Real m = floor((deg - d) * 60);
    Real s = (deg - d - m / 60) * 3600;
    return Vec3(d, m, s);
  }

  /// @brief Convert degrees, minutes, and seconds to degrees
  /// @param dms Angle in degrees, minutes, and seconds
  /// @return Angle in degrees
  Real DegMinSecToDeg(const Vec3& dms) {
    Real decdeg = dms(0) + dms(1) / 60.0 + dms(2) / 3600.0;
    return decdeg;
  }

  Real sind(Real x) { return sin(x * RAD); }

  Real cosd(Real x) { return cos(x * RAD); }

  Real tand(Real x) { return tan(x * RAD); }

  MatX SampleMvNormal(const VecX& mean, const MatX& cov, int nn, std::mt19937* rng) {
    int n = mean.size();
    Eigen::LLT<MatXd> llt(cov.cast<double>());
    MatXd L = llt.matrixL();
    std::normal_distribution<double> dist(0.0, 1.0);
    std::mt19937& gen = rng ? *rng : lupnt::RandomEngine::Get();
    MatX samples(nn, n);
    for (int i = 0; i < nn; i++) {
      VecXd z(n);
      for (int j = 0; j < n; j++) z(j) = dist(gen);
      samples.row(i) = mean + (L * z).cast<Real>();
    }
    return samples;
  }

  Real SampleNormal(Real mean, Real stddev, std::mt19937* rng) {
    std::normal_distribution<double> dist(0.0, 1.0);
    std::mt19937& gen = rng ? *rng : lupnt::RandomEngine::Get();
    return dist(gen) * stddev + mean;
  }

  /// @brief Compute the safe acos function
  /// @param x Input value
  /// @return acos(x) if x is in [-1, 1], otherwise acos(x - EPS) or acos(x + EPS)
  Real safe_acos(Real x) {
    if (x >= 1.0) {
      return acos(x - EPS);
    } else if (x <= -1.0) {
      return acos(x + EPS);
    } else {
      return acos(x);
    }
  }

  /// @brief Compute the safe asin function
  /// @param x Input value
  /// @return asin(x) if x is in [-1, 1], otherwise asin(x - EPS) or asin(x + EPS)
  Real safe_asin(Real x) {
    if (x >= 1.0) {
      return asin(x - 1e-16);
    } else if (x <= -1.0) {
      return asin(x + 1e-16);
    } else {
      return asin(x);
    }
  }

  /// @brief Passive rotation matrix about the x-axis
  /// @param angle Angle in radians
  /// @return Rotation matrix
  Mat3 RotX(Real angle) {
    Real c = cos(angle);
    Real s = sin(angle);
    Mat3 R{
        {1.0, 0.0, 0.0},
        {0.0, c, s},
        {0.0, -s, c},
    };
    return R;
  }

  /// @brief Passive rotation matrix about the y-axis
  /// @param angle Angle in radians
  /// @return Rotation matrix
  Mat3 RotY(Real angle) {
    Real c = cos(angle);
    Real s = sin(angle);
    Mat3 R{
        {c, 0.0, -s},
        {0.0, 1.0, 0.0},
        {s, 0.0, c},
    };
    return R;
  }

  /// @brief Passive rotation matrix about the z-axis
  /// @param angle Angle in radians
  /// @return Rotation matrix
  Mat3 RotZ(Real angle) {
    Real c = cos(angle);
    Real s = sin(angle);
    Mat3 R{
        {c, s, 0.0},
        {-s, c, 0.0},
        {0.0, 0.0, 1.0},
    };
    return R;
  }

  Mat3 RotXdot(Real angle, Real angle_dot) {
    Real c = cos(angle);
    Real s = sin(angle);
    Mat3 R_dot{
        {0.0, 0.0, 0.0},
        {0.0, -s * angle_dot, c * angle_dot},
        {0.0, -c * angle_dot, -s * angle_dot},
    };
    return R_dot;
  }

  Mat3 RotYdot(Real angle, Real angle_dot) {
    Real c = cos(angle);
    Real s = sin(angle);
    Mat3 R_dot{
        {-s * angle_dot, 0.0, -c * angle_dot},
        {0.0, 0.0, 0.0},
        {c * angle_dot, 0.0, -s * angle_dot},
    };
    return R_dot;
  }

  Mat3 RotZdot(Real angle, Real angle_dot) {
    Real c = cos(angle);
    Real s = sin(angle);
    Mat3 R_dot{
        {-s * angle_dot, c * angle_dot, 0.0},
        {-c * angle_dot, -s * angle_dot, 0.0},
        {0.0, 0.0, 0.0},
    };
    return R_dot;
  }

  /// @brief Skew symmetric matrix from a vector
  /// @param x Input vector
  /// @return Skew symmetric matrix
  Mat3 Skew(const Vec3& x) {
    Mat3 skew{
        {0.0, -x(2), x(1)},
        {x(2), 0.0, -x(0)},
        {-x(1), x(0), 0.0},
    };
    return skew;
  }

  std::vector<double> EigenToVector(const VecX& x) {
    std::vector<double> y(x.size());
    for (int i = 0; i < x.size(); i++) {
      y[i] = x(i).val();
    }
    return y;
  }

  /// @brief Local function for use by RatioOfSectorToTriangleArea
  ///        F = 1 - eta +(m/eta^2)*W(m/eta^2-l)
  /// @param eta
  /// @param m
  /// @param l
  /// @return
  Real F(Real eta, Real m, Real l) {
    const double eps = 100.0 * EPS;
    Real w, W, a, n, g;
    w = m / (eta * eta) - l;
    if (abs(w) < 0.1) {  // Series expansion
      W = a = 4.0 / 3.0;
      n = 0.0;
      do {
        n += 1.0;
        a *= w * (n + 2.0) / (n + 1.5);
        W += a;
      } while (abs(a) >= eps);
    } else {
      if (w > 0.0) {
        g = 2.0 * asin(sqrt(w));
        W = (2.0 * g - sin(2.0 * g)) / pow(sin(g), 3);
      } else {
        g = 2.0 * log(sqrt(-w) + sqrt(1.0 - w));  // =2.0*arsinh(sqrt(-w))
        W = (sinh(2.0 * g) - 2.0 * g) / pow(sinh(g), 3);
      }
    }
    return (1.0 - eta + (w + l) * W);
  }

  /// @brief Compute the ratio of the sector area to the triangle area
  /// @param x First vector
  /// @param y Second vector
  /// @param tau Time of flight
  /// @ref
  /// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
  /// applications. Berlin : New York: Springer, 2000.
  /// doi: 10.1007/978-3-642-58351-3.
  Real RatioOfSectorToTriangleArea(const Vec3& r1, const Vec3& r2, Real tau) {
    const int max_it = 30;
    const double delta = 100 * EPS;

    Real s1 = r1.norm();
    Real s2 = r2.norm();
    Real kappa = sqrt(2.0 * (s1 * s2 + r1.dot(r2)));
    Real m = tau * tau / pow(kappa, 3);
    Real l = (s1 + s2) / (2.0 * kappa) - 0.5;
    Real eta_min = sqrt(m / (l + 1.0));

    // Start with Hansen's approximation
    Real eta2 = (12.0 + 10.0 * sqrt(1.0 + (44.0 / 9.0) * m / (l + 5.0 / 6.0))) / 22.0;
    Real eta1 = eta2 + 0.1;

    // Secant method
    Real F1 = F(eta1, m, l);
    Real F2 = F(eta2, m, l);

    int it = 0;
    while (abs(F2 - F1) > delta && it < max_it) {
      Real d_eta = -F2 * (eta2 - eta1) / (F2 - F1);
      eta1 = eta2;
      F1 = F2;
      while (eta2 + d_eta <= eta_min) d_eta *= 0.5;
      eta2 += d_eta;
      F2 = F(eta2, m, l);
      ++it;
    }
    if (it >= max_it) throw std::runtime_error("Hansen's method did not converge");
    return eta2;
  }

  VecXd SolveLinearEqSVD(const MatXd& A, const VecXd& b) {
    Eigen::JacobiSVD<MatXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    return svd.solve(b);
  }

  // Compute A * X = B
  MatXd SolveLinearEqSVD(const MatXd& A, const MatXd& B) {
    Eigen::JacobiSVD<MatXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    MatXd X(A.cols(), B.cols());
    for (int i = 0; i < B.cols(); i++) {
      X.col(i) = svd.solve(B.col(i));
    }
    return X;
  }

  MatXd PseudoInverse(const MatXd& A) {
    Eigen::JacobiSVD<MatXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    double tol = 1e-6;
    MatXd S = svd.singularValues();
    // `S` is a k x 1 column of singular values; the reciprocal-diagonal matrix must be
    // k x k so `V * S_inv * U^T` is conformable (Zero(S.rows(), S.cols()) would be k x 1).
    MatXd S_inv = MatXd::Zero(S.rows(), S.rows());
    for (int i = 0; i < S.rows(); i++) {
      if (S(i) > tol) {
        S_inv(i, i) = 1.0 / S(i);
      }
    }
    MatXd S_pinv = svd.matrixV() * S_inv * svd.matrixU().transpose();
    return S_pinv;
  }

  VecX JacobianParallel(const std::function<VecX(const VecX&)>& func, const VecX& x0, MatXd& J) {
    int n = x0.size();
    std::vector<MatXd> J_vec(n);
    VecX xf;

#pragma omp parallel for
    for (int i = 0; i < n; i++) {
      VecX xi;
      MatXd Ji;
      VecX x0_tmp = x0.cast<double>();
      jacobian(func, wrt(x0_tmp(i)), at(x0_tmp), xi, Ji);
      if (i == 0) xf = xi;
      J_vec[i] = Ji;
    }

    int m = xf.size();
    J.resize(m, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        J(j, i) = J_vec[i](j);
      }
    }
    return xf;
  }

  /// @brief Compute the vech of a matrix, defined as the lower triangular part of
  /// the matrix
  ///        stacked column-wise
  /// @param x Input matrix
  /// @return Vech of the matrix
  VecX Vech(const MatX& x) {
    int n = x.rows();
    VecX v(n * (n + 1) / 2);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j <= i; j++) {
        v(i * (i + 1) / 2 + j) = x(i, j);
      }
    }
    return v;
  }

  Real RotationAngle(const Mat3& R) { return acos((R.trace() - 1.0) / 2.0); }
}  // namespace lupnt

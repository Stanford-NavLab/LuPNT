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

#include "lupnt/core/constants.h"

namespace lupnt {
  /// @brief Create a range of values
  /// @param start Start value
  /// @param stop Stop value
  /// @param step Step value
  /// @return Vector of values
  template <typename T> VectorX<T> arange(T start, T stop, T step) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step) {
      values.push_back(value);
    }
    return Eigen::Map<VectorX<T>>(values.data(), values.size());
  }
  template VectorX<int> arange<int>(int start, int stop, int step);
  template VectorX<double> arange<double>(double start, double stop, double step);
  template VectorX<Real> arange<Real>(Real start, Real stop, Real step);

  /// @brief Compute the angle between two vectors
  /// @param x First vector
  /// @param y Second vector
  /// @return Angle between the two vectors in radians
  Real AngleBetweenVecs(const VecX& x, const VecX& y) {
    return 2.0
           * atan2((x.normalized() - y.normalized()).norm(),
                   (x.normalized() + y.normalized()).norm());
  }

  /// @brief Wrap angle in radians to [-pi, pi]
  /// @param angle Angle in radians
  /// @return Wrapped angle in radians
  Real Wrap2Pi(Real angle) { return atan2(sin(angle), cos(angle)); }
  VecX Wrap2Pi(VecX angle) {
    int n = angle.size();
    VecX out(n);
    for (int i = 0; i < n; i++) out(i) = Wrap2Pi(angle(i));
    return out;
  }

  /// @brief Wrap angle in radians to [0, 2pi]
  /// @param angle Angle in radians
  /// @return Wrapped angle in radians
  Real Wrap2TwoPi(Real angle) { return angle - TWO_PI * floor(angle / TWO_PI); }
  VecX Wrap2TwoPi(VecX angle) {
    int n = angle.size();
    VecX out(n);
    for (int i = 0; i < n; i++) out(i) = Wrap2TwoPi(angle(i));
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
  Real Decimal2Decibel(Real x) { return 10 * log10(x); }
  VecX Decimal2Decibel(VecX x) {
    int n = x.size();
    VecX out(n);
    for (int i = 0; i < n; i++) out(i) = Decimal2Decibel(x(i));
    return out;
  }
  MatX Decimal2Decibel(MatX x) {
    MatX out(x.rows(), x.cols());
    for (int i = 0; i < x.rows(); i++)
      for (int j = 0; j < x.cols(); j++) out(i, j) = Decimal2Decibel(x(i, j));
    return out;
  }

  /// @brief Convert decibel value to decimal
  /// @param x Decibel value
  /// @return Decimal value
  Real Decibel2Decimal(Real x) { return pow(10, x / 10); }
  VecX Decibel2Decimal(VecX x) {
    VecX out(x.size());
    for (int i = 0; i < x.size(); i++) out(i) = Decibel2Decimal(x(i));
    return out;
  }
  MatX Decibel2Decimal(MatX x) {
    MatX out(x.rows(), x.cols());
    for (int i = 0; i < x.rows(); i++)
      for (int j = 0; j < x.cols(); j++) out(i, j) = Decibel2Decimal(x(i, j));
    return out;
  }

  /// @brief Convert degrees to degrees, minutes, and seconds
  /// @param deg Angle in degrees
  /// @return Angle in degrees, minutes, and seconds
  Vec3 Degrees2DegMinSec(Real deg) {
    Real d = floor(deg);
    Real m = floor((deg - d) * 60);
    Real s = (deg - d - m / 60) * 3600;
    return Vec3(d, m, s);
  }

  /// @brief Convert degrees, minutes, and seconds to degrees
  /// @param dms Angle in degrees, minutes, and seconds
  /// @return Angle in degrees
  Real DegMinSec2Degrees(Vec3 dms) {
    Real decdeg = dms(0) + dms(1) / 60.0 + dms(2) / 3600.0;
    return decdeg;
  }

  Real sind(Real x) { return sin(x * RAD); }

  Real cosd(Real x) { return cos(x * RAD); }

  Real tand(Real x) { return tan(x * RAD); }

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

  double J0Bessel(double x) {
    // Bessel function of the first kind of order 0
    // Reference: https://en.wikipedia.org/wiki/Bessel_function
    double J0 = 0.0;
    double y = 1.0;
    double sum = 1.0;
    for (int i = 1; i < 10; i++) {
      y = y * x * x / (4 * i * i);
      sum += y;
    }

    return sum;
  }

  double J1Bessel(double x) {
    // Bessel function of the first kind of order 1
    // Reference: https://en.wikipedia.org/wiki/Bessel_function
    double J1 = 0.0;
    double y = 1.0;
    double sum = 1.0;
    for (int i = 1; i < 10; i++) {
      y = y * x / (2 * i * (2 * i + 1));
      sum += y;
    }

    return sum;
  }

  /// @brief Compute the root mean square of a vector
  /// @param x Input vector
  /// @return Root mean square of the vector
  Real RootMeanSquare(VecX x) { return sqrt(x.array().pow(2).sum() / x.size()); }
  double RootMeanSquareD(VecXd x) { return sqrt(x.array().pow(2).sum() / x.size()); }

  /// @brief Compute the pth percentile of a vector
  /// @param x Input vector
  /// @param p Percentile value
  Real Percentile(VecX x, double p) {
    Real* start = x.data();
    Real* end = x.data() + x.size();
    std::sort(start, end);
    int index = std::ceil(p * x.size());
    if (index > (x.size() - 1)) {
      index = x.size() - 1;
    }
    return x(index);
  }

  double PercentileD(VecXd x, double p) {
    double* start = x.data();
    double* end = x.data() + x.size();
    std::sort(start, end);
    int index = std::ceil(p * x.size());
    if (index > (x.size() - 1)) {
      index = x.size() - 1;
    }
    return x(index);
  }

  /// @brief Compute the standard deviation of a vector
  /// @param x Input vector
  /// @return Standard deviation of the vector
  Real Std(VecX x) {
    Real mean = x.sum() / x.size();
    Real sq_sum = 0.0;
    for (int i = 0; i < x.size(); i++) {
      sq_sum += (x(i) - mean) * (x(i) - mean);
    }
    return sqrt(sq_sum / x.size());
  }

  double StdD(VecXd x) {
    double mean = x.sum() / x.size();
    double sq_sum = 0.0;
    for (int i = 0; i < x.size(); i++) {
      sq_sum += (x(i) - mean) * (x(i) - mean);
    }
    return sqrt(sq_sum / x.size());
  }

  double erfc(double x) { return 1 - erf(x); }

  double qfunc(double x) { return 0.5 * erfc(x / sqrt(2)); }

  /// @brief Sample from a multivariate normal distribution
  /// @param mean Mean vector
  /// @param covar Covariance matrix
  /// @param nn Number of samples
  /// @param seed Random seed
  MatX SampleMVN(const VecX mean, const MatX covar, int nn, int seed) {
    // Define random generator with Gaussian distribution
    int xsize = mean.size();
    auto generator = std::mt19937(seed);
    auto dist = std::bind(std::normal_distribution<double>{0.0, 1.0}, generator);

    // Transform Mat
    MatX normTransform(xsize, xsize);
    Eigen::LLT<MatX> cholSolver(covar);

    if (cholSolver.info() == Eigen::Success) {
      // Use cholesky solver
      normTransform = cholSolver.matrixL();
    } else {
      std::runtime_error("The covariance matrix must be symmetric and pos-definite.");
    }

    MatX randN(xsize, nn);
    MatX mean_samples(xsize, nn);
    for (int i = 0; i < xsize; i++) {
      for (int j = 0; j < nn; j++) {
        randN(i, j) = dist();
        mean_samples(i, j) = mean(i);
      }
    }

    MatX samples = normTransform * randN + mean_samples;
    return samples;
  };

  double SampleRandNormal(double mean, double std, int seed) {
    auto generator = std::mt19937(seed);
    auto dist = std::bind(std::normal_distribution<double>{mean, std}, generator);
    return dist();
  }

  /// @brief Create a block diagonal matrix from two matrices
  /// @param A First matrix
  /// @param B Second matrix
  /// @return Block diagonal matrix
  MatX BlkDiag(const MatX& A, const MatX& B) {
    MatX C(A.rows() + B.rows(), A.cols() + B.cols());
    C << A, MatX::Zero(A.rows(), B.cols()), MatX::Zero(B.rows(), A.cols()), B;
    return C;
  }

  MatXd BlkDiagD(const MatXd& A, const MatXd& B) {
    int rows = A.rows() + B.rows();

    MatXd C(rows, rows);
    C = MatXd::Zero(rows, rows);
    C.block(0, 0, A.rows(), A.cols()) = A;
    C.block(A.rows(), A.cols(), B.rows(), B.cols()) = B;
    return C;
  }

  /// @brief Active rotation matrix about the x-axis
  /// @param angle Angle in radians
  /// @return Rotation matrix
  template <typename T> Matrix<T, 3, 3> RotX(T angle) {
    T c = cos(angle);
    T s = sin(angle);
    Matrix<T, 3, 3> R{
        {1.0, 0.0, 0.0},
        {0.0, c, s},
        {0.0, -s, c},
    };
    return R;
  }
  template Mat3d RotX(double angle);
  template Mat3 RotX(Real angle);

  /// @brief Active rotation matrix about the y-axis
  /// @param angle Angle in radians
  /// @return Rotation matrix
  template <typename T> Matrix<T, 3, 3> RotY(T angle) {
    T c = cos(angle);
    T s = sin(angle);
    Matrix<T, 3, 3> R{
        {c, 0.0, -s},
        {0.0, 1.0, 0.0},
        {s, 0.0, c},
    };
    return R;
  }
  template Mat3d RotY(double angle);
  template Mat3 RotY(Real angle);

  /// @brief Active rotation matrix about the z-axis
  /// @param angle Angle in radians
  /// @return Rotation matrix
  template <typename T> Matrix<T, 3, 3> RotZ(T angle) {
    T c = cos(angle);
    T s = sin(angle);
    Matrix<T, 3, 3> R{
        {c, s, 0.0},
        {-s, c, 0.0},
        {0.0, 0.0, 1.0},
    };
    return R;
  }
  template Mat3d RotZ(double angle);
  template Mat3 RotZ(Real angle);

  /// @brief Skew symmetric matrix from a vector
  /// @param x Input vector
  /// @return Skew symmetric matrix
  template <typename T> Matrix<T, 3, 3> Skew(Vector<T, 3> x) {
    Matrix<T, 3, 3> skew{
        {0.0, -x(2), x(1)},
        {x(2), 0.0, -x(0)},
        {-x(1), x(0), 0.0},
    };
    return skew;
  }
  template Mat3d Skew(Vec3d x);
  template Mat3 Skew(Vec3 x);

  /// @brief Convert a vector of floats to a vector of doubles
  /// @param x Input vector
  /// @return Vector of doubles
  VecXd ToDouble(const VecX& x) {
    VecXd y = x.cast<double>();
    return y;
  }

  /// @brief Convert a matrix of floats to a matrix of doubles
  /// @param x Input matrix
  /// @return Matrix of doubles
  MatXd ToDouble(const MatX& x) {
    MatXd y = x.cast<double>();
    return y;
  }

  std::vector<double> ToDoubleVec(const VecX& x) {
    std::vector<double> y(x.size());
    for (size_t i = 0; i < x.size(); i++) {
      y[i] = x(i).val();
    }
    return y;
  }

  std::vector<double> ToDoubleVec(const VecXd& x) {
    std::vector<double> y(x.size());
    for (size_t i = 0; i < x.size(); i++) {
      y[i] = x(i);
    }
    return y;
  }

  std::vector<double> ToDoubleVec(const VecXi& x) {
    std::vector<double> y(x.size());
    for (size_t i = 0; i < x.size(); i++) {
      y[i] = x(i);
    }
    return y;
  }

  std::vector<double> ToDoubleVec(const std::vector<Real>& x) {
    std::vector<double> y(x.size());
    for (size_t i = 0; i < x.size(); i++) {
      y[i] = x[i].val();
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
  Real RatioOfSectorToTriangleArea(Vec3 r1, Vec3 r2, Real tau) {
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
    assert(it < max_it && "Hansen's method did not converge");
    return eta2;
  }

  /// @brief Create a vector of evenly spaced values
  VecX arange(Real start, Real stop, Real step) {
    int n = static_cast<int>((stop - start) / step);
    VecX v(n);
    for (int i = 0; i < n; i++) {
      v(i) = start + i * step;
    }
    return v;
  }

}  // namespace lupnt

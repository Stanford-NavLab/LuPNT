#pragma once

#include <random>
#include <tuple>
#include <utility>

#include "lupnt/core/definitions.h"

namespace lupnt {

  /// @brief Implementation helper for `Unpack`: expands a fixed-size Eigen vector into a
  /// `std::tuple` of its elements via an index sequence.
  template <typename Vec, std::size_t... Indices>
  auto UnpackImpl(const Vec& vec, std::index_sequence<Indices...>) {
    return std::make_tuple(vec(Indices)...);
  }

  /// @brief Unpack a fixed-size Eigen vector into a `std::tuple` of its scalar elements,
  /// for convenient structured-binding access (e.g. `auto [a, e, i, Om, om, M] =
  /// Unpack(Vec6(coe))`).
  ///
  /// Used across `conversions/state_conversions.cc` and example scripts to destructure
  /// orbital-element / state vectors into named scalars.
  ///
  /// @param vec  Fixed-size Eigen column vector (`Size` > 0)
  /// @return     `std::tuple` of `Size` scalars, one per element of `vec`
  template <typename T, int Size> auto Unpack(const Eigen::Matrix<T, Size, 1>& vec) {
    static_assert(Size > 0, "Cannot unpack a zero-sized vector.");
    return UnpackImpl(vec, std::make_index_sequence<Size>{});
  }

  /// @brief Generate a vector of evenly-spaced values `[start, start+step, ..., < stop)`,
  /// analogous to NumPy's `arange`.
  ///
  /// Used throughout example/scenario setup code (e.g. `ex_integration.cc`,
  /// `ex_n_body_dynamics.cc`) to build the time-span vector `tspan`/`times` driving a
  /// propagation loop.
  ///
  /// @param start  First value
  /// @param stop   Exclusive upper bound
  /// @param step   Increment (default 1)
  /// @return       Vector of values from `start` up to (excluding) `stop`
  template <typename T> VectorX<T> Arange(T start, T stop, T step = 1);

  /// @brief Take every `step`-th element of vector `x`.
  ///
  /// Used by plotting/post-processing code (e.g. `ex_imu_dynamics.cc`, `ex_ekf_2d.cc`) to
  /// thin a dense time-history vector before plotting.
  ///
  /// @param x     Input vector
  /// @param step  Stride (every `step`-th element is kept)
  /// @return      Subsampled vector of size `x.size() / step`
  VecXd Subsample(const VecX& x, int step);

  /// @brief Take every `step_row`-th row and `step_col`-th column of matrix `x`.
  ///
  /// Used alongside the vector overload to thin a dense state-history matrix (rows =
  /// time samples, columns = state components) before plotting.
  ///
  /// @param x         Input matrix
  /// @param step_row  Row stride
  /// @param step_col  Column stride
  /// @return          Subsampled matrix of size `(x.rows()/step_row) x (x.cols()/step_col)`
  MatXd Subsample(const MatX& x, int step_row, int step_col);

  /// @brief Angle between two vectors [rad], computed via the numerically robust
  /// `2*atan2(|n1-n2|, |n1+n2|)` formula (avoids the precision loss of `acos(dot)` near 0
  /// and pi), where `n1`, `n2` are the normalized inputs.
  ///
  /// @param x  First vector
  /// @param y  Second vector
  /// @return   Angle between `x` and `y` [rad], in `[0, pi]`
  Real AngleBetweenVecs(const VecX& x, const VecX& y);

  /// @brief Row-wise angle between corresponding rows of `x` and `y` [rad].
  ///
  /// Matrix overload of `AngleBetweenVecs`; computed via `acos` of the row-wise
  /// dot-product of normalized rows (clamped to `[-1, 1]`).
  ///
  /// @param x  Matrix whose rows are vectors to compare
  /// @param y  Matrix whose rows are vectors to compare (same shape as `x`)
  /// @return   Vector of per-row angles [rad], in `[0, pi]`
  VecX AngleBetweenVecs(const MatX& x, const MatX& y);

  /// @brief Convert an angle given in degrees/arcminutes/arcseconds to decimal degrees.
  ///
  /// @param degrees  Whole-degree component
  /// @param minutes  Arcminute component
  /// @param seconds  Arcsecond component
  /// @return         Angle in decimal degrees
  double DegMinSec2DeciDeg(double degrees, double minutes, double seconds);

  /// @brief Wrap an angle to the range `(-pi, pi]` via `atan2(sin(angle), cos(angle))`.
  ///
  /// Used pervasively in `conversions/anomaly_conversions.cc` (e.g.
  /// `EccToMeanAnomaly`, `MeanToEccAnomaly`) and time conversions (e.g.
  /// `EarthRotationAngle`) to keep angular quantities in a canonical range.
  ///
  /// @param angle  Angle [rad]
  /// @return       Wrapped angle [rad] in `(-pi, pi]`
  Real WrapToPi(Real angle);

  /// @brief Wrap an angle to the range `[0, 2*pi)`.
  ///
  /// @param angle  Angle [rad]
  /// @return       Wrapped angle [rad] in `[0, 2*pi)`
  Real WrapToTwoPi(Real angle);

  /// @brief Element-wise `VecX` overload of `WrapToPi`.
  VecX WrapToPi(const VecX& angle);

  /// @brief Element-wise `VecX` overload of `WrapToTwoPi`; used e.g. by
  /// `Antenna::ComputeGain` to normalize a query azimuth angle before pattern lookup.
  VecX WrapToTwoPi(const VecX& angle);

  /// @brief Convert a decimal (linear) ratio to decibels: `10*log10(x)`.
  ///
  /// @param x  Decimal (linear-scale) ratio, e.g. a power ratio
  /// @return   Value in decibels [dB]
  Real DecimalToDecibel(Real x);

  /// @brief Element-wise `ArrX` overload of `DecimalToDecibel`.
  ArrX DecimalToDecibel(const ArrX& x);

  /// @brief Convert a value in decibels to a decimal (linear) ratio: `10^(x/10)`.
  ///
  /// Used by the link-budget / C/N0 computations (e.g. `GnssConstellation`'s C/N0 models,
  /// `comms_utils.cc`'s FLL/PLL tracking-loop equations) to convert a C/N0 value from dB-Hz
  /// to linear (Hz) units before further computation.
  ///
  /// @param x  Value in decibels [dB]
  /// @return   Decimal (linear-scale) ratio
  Real DecibelToDecimal(Real x);

  /// @brief Element-wise `ArrX` overload of `DecibelToDecimal`.
  ArrX DecibelToDecimal(const ArrX& x);

  /// @brief Return the larger of two autodiff scalars, comparing by underlying value
  /// (`.val()`) and returning the corresponding `Real` (preserving its derivative
  /// information).
  ///
  /// Used e.g. by `environment/forces.cc` to clamp an intermediate quantity to be
  /// non-negative before taking its square root (`sqrt(Max(y2, 0.0))`).
  ///
  /// @param x  First value
  /// @param y  Second value
  /// @return   `x` if `x.val() > y.val()`, else `y`
  Real Max(Real x, Real y);

  /// @brief Return the smaller of two autodiff scalars; see `Max` for comparison
  /// semantics.
  Real Min(Real x, Real y);

  /// @brief Plain-`double` overload of `Max`.
  double MaxD(double x, double y);

  /// @brief Plain-`double` overload of `Min`.
  double MinD(double x, double y);

  /// @brief Round `x` to `n` decimal places, preserving `x`'s autodiff derivative
  /// information (only the value is rounded).
  ///
  /// @param x  Value to round
  /// @param n  Number of decimal places (default 0)
  /// @return   Rounded value
  Real round(Real x, int n = 0);

  /// @brief Fractional part of `x` (`x - floor(x)`), preserving autodiff derivative
  /// information.
  Real frac(Real x);

  /// @brief Ceiling of `x`, preserving autodiff derivative information.
  Real ceil(Real x);

  /// @brief Floor of `x`, preserving autodiff derivative information.
  Real floor(Real x);

  /// @brief Floating-point modulus `x mod y` (`std::fmod`), preserving `x`'s autodiff
  /// derivative information.
  Real mod(Real x, Real y);

  /// @brief Convert an angle in decimal degrees to degrees/arcminutes/arcseconds.
  ///
  /// @param deg  Angle in decimal degrees
  /// @return     `(degrees, arcminutes, arcseconds)` packed as `Vec3`
  Vec3 DegToDegMinSec(Real deg);

  /// @brief Convert an angle given as degrees/arcminutes/arcseconds to decimal degrees.
  ///
  /// @param hms  `(degrees, arcminutes, arcseconds)` packed as `Vec3`
  /// @return     Angle in decimal degrees
  Real DegMinSecToDeg(const Vec3& hms);

  /// @brief Sine of an angle given in degrees.
  Real sind(Real x);

  /// @brief Cosine of an angle given in degrees.
  Real cosd(Real x);

  /// @brief Tangent of an angle given in degrees.
  Real tand(Real x);

  /// @brief Numerically robust `acos`, clamping the input to `[-1, 1]` (offset by `EPS`)
  /// before evaluating, to avoid NaNs from values that are out of range only due to
  /// floating-point round-off.
  ///
  /// Used throughout geometry/elevation-angle computations -- e.g.
  /// `GnssMeasurement`'s elevation/off-boresight angles, `GnssAttitude`'s yaw/attitude
  /// angles, and orbital-element conversions (inclination from angular-momentum
  /// components) -- wherever the `acos` argument is a dot product of unit vectors that may
  /// slightly exceed +/-1 due to round-off.
  ///
  /// @param x  Cosine value, nominally in `[-1, 1]`
  /// @return   `acos(x)`, with `x` clamped to `(-1, 1)` if it lies (slightly) outside
  Real safe_acos(Real x);

  /// @brief Numerically robust `asin`, analogous to `safe_acos`.
  ///
  /// @param x  Sine value, nominally in `[-1, 1]`
  /// @return   `asin(x)`, with `x` clamped to `(-1, 1)` if it lies (slightly) outside
  Real safe_asin(Real x);

  /// @brief Bessel function of the first kind, order 0, via its truncated (10-term) power
  /// series.
  ///
  /// @param x  Argument
  /// @return   J0(x)
  /// @note https://en.wikipedia.org/wiki/Bessel_function
  template <typename T> T J0Bessel(T x) {
    // J0(x) = sum_{k>=0} (-1)^k (x/2)^{2k} / (k!)^2; term ratio t_k/t_{k-1} = -x^2/(4 k^2).
    T y = 1.0;
    T sum = 1.0;
    for (int i = 1; i < 10; i++) {
      y = -y * x * x / (4 * i * i);
      sum += y;
    }
    return sum;
  }

  /// @brief Bessel function of the first kind, order 1, via its truncated (10-term) power
  /// series.
  ///
  /// @param x  Argument
  /// @return   J1(x)
  /// @note https://en.wikipedia.org/wiki/Bessel_function
  template <typename T> T J1Bessel(T x) {
    // J1(x) = sum_{k>=0} (-1)^k (x/2)^{2k+1} / (k! (k+1)!); first term x/2,
    // term ratio t_k/t_{k-1} = -x^2 / (4 k (k+1)).
    T y = x / 2;
    T sum = y;
    for (int i = 1; i < 10; i++) {
      y = -y * x * x / (4 * i * (i + 1));
      sum += y;
    }
    return sum;
  }
  /// @brief Root-mean-square value of a vector/array: `sqrt(sum(x.^2) / numel(x))`.
  ///
  /// Used by `filters/filter_print.cc` to summarize per-component estimation-error
  /// statistics (alongside `Std` and `Percentile`) when printing filter results.
  ///
  /// @param x  Input vector/array
  /// @return   Root-mean-square value of `x`
  template <typename T> typename T::Scalar RootMeanSquare(const DenseBase<T>& x) {
    return sqrt(x.derived().squaredNorm() / x.size());
  }

  /// @brief `p`-th percentile of a vector/array, computed by sorting the values and
  /// indexing at `ceil(p * (n-1))`.
  ///
  /// Used by `filters/filter_print.cc` to report e.g. the 68th/95th/99th percentile of
  /// filter estimation errors.
  ///
  /// @param x  Input vector/array
  /// @param p  Percentile as a fraction in `[0, 1]` (e.g. 0.95 for the 95th percentile)
  /// @return   The `p`-th percentile value of `x`
  template <typename T> typename T::Scalar Percentile(const DenseBase<T>& x, double p) {
    auto evaluated = x.eval();
    std::vector<typename T::Scalar> data(evaluated.data(), evaluated.data() + evaluated.size());
    std::sort(data.begin(), data.end());
    size_t index = std::ceil(p * (data.size() - 1));
    return data[index];
  }

  /// @brief Sample standard deviation of a vector/array (Bessel-corrected, dividing by
  /// `n-1`).
  ///
  /// Used by `filters/filter_print.cc` to report the dispersion of filter estimation
  /// errors.
  ///
  /// @param x  Input vector/array
  /// @return   Sample standard deviation of `x`
  template <typename T> typename T::Scalar Std(const DenseBase<T>& x) {
    auto mean_val = x.mean();
    return sqrt((x.derived().array() - mean_val).square().sum() / (x.size() - 1));
  }

  /// @brief Complementary error function, `erfc(x) = 1 - erf(x)`.
  template <typename T> T erfc(T x) { return 1 - erf(x.val()); }

  /// @brief Gaussian Q-function (upper-tail probability of a standard normal),
  /// `qfunc(x) = 0.5*erfc(x/sqrt(2))`.
  template <typename T> T qfunc(T x) { return 0.5 * erfc(x / sqrt(2)); }

  /// @brief Draw `nn` samples from a multivariate normal distribution `N(mean, cov)`.
  ///
  /// Used by measurement-simulation and filter-test code (e.g. `example_adaptive.cc`) to
  /// generate synthetic state/measurement-error realizations from a given covariance
  /// matrix, via a Cholesky factorization of `cov`.
  ///
  /// @param mean  Mean vector [size n]
  /// @param cov   Covariance matrix [size n x n], must be positive-definite
  /// @param nn    Number of samples to draw (default 1)
  /// @param rng   Optional random engine; defaults to `lupnt::RandomEngine::Get()`
  /// @return      `nn x n` matrix, each row a sample from `N(mean, cov)`
  MatX SampleMvNormal(const VecX& mean, const MatX& cov, int nn = 1, std::mt19937* rng = nullptr);

  /// @brief Draw a single sample from a scalar normal distribution `N(mean, std^2)`.
  ///
  /// Used throughout example/scenario code to generate scalar noise realizations (e.g. for
  /// clock or measurement noise models).
  ///
  /// @param mean  Mean (default 0.0)
  /// @param std   Standard deviation (default 1.0)
  /// @param rng   Optional random engine; defaults to `lupnt::RandomEngine::Get()`
  /// @return      Sampled value
  Real SampleNormal(Real mean = 0.0, Real std = 1.0, std::mt19937* rng = nullptr);

  /// @brief Form the block-diagonal matrix `[[A, 0], [0, B]]` from two matrices.
  ///
  /// Used by `filters/filter_utils.cc::InitialCovariancePosVelClock` to combine an
  /// independent position/velocity covariance block with a clock bias/drift covariance
  /// block into the joint filter initial-covariance matrix.
  ///
  /// @param A  Top-left block (size N1 x M1)
  /// @param B  Bottom-right block (size N2 x M2)
  /// @return   Block-diagonal matrix of size (N1+N2) x (M1+M2), with off-diagonal blocks
  ///           zero
  template <typename T, int N1, int M1, int N2, int M2>
  Matrix<T, N1 + N2, M1 + N2> BlockDiagonal(const Matrix<T, N1, M1>& A,
                                            const Matrix<T, N2, M2>& B) {
    Matrix<T, N1 + N2, M1 + M2> C;
    C << A, Matrix<T, N1, M2>::Zero(), Matrix<T, N2, M1>::Zero(), B;
    return C;
  }

  /// @brief Passive (frame-rotation) rotation matrix about the x-axis by `angle`.
  ///
  /// Rotates the coordinate frame (not the vector) by `angle` about x, i.e. transforms a
  /// vector's components from the unrotated frame to the frame rotated by `angle` about x.
  /// Used as a building block for composite Euler-angle rotations throughout
  /// `conversions/frame_conversions.cc` (e.g. polar-motion matrix `RotPolarMotion`, body-to-
  /// ENU rotations, lunar PA/CI Euler-angle reconstruction).
  ///
  /// @param angle  Rotation angle [rad]
  /// @return       3x3 passive rotation matrix about x
  Mat3 RotX(Real angle);

  /// @brief Passive (frame-rotation) rotation matrix about the y-axis by `angle`. See
  /// `RotX` for sign/passive convention.
  Mat3 RotY(Real angle);

  /// @brief Passive (frame-rotation) rotation matrix about the z-axis by `angle`. See
  /// `RotX` for sign/passive convention. Used e.g. by `EarthSiderealRotation`/
  /// `RotSideralMotion` to build the Earth-rotation-angle rotation `Rz(theta_ERA)`.
  Mat3 RotZ(Real angle);

  /// @brief Time derivative of `RotX(angle)` with respect to time, given the angle's time
  /// derivative `angle_dot`: `d/dt RotX(angle(t))`.
  ///
  /// @param angle      Rotation angle [rad]
  /// @param angle_dot  Time derivative of `angle` [rad/s]
  /// @return           3x3 matrix `d(RotX)/dt`
  Mat3 RotXdot(Real angle, Real angle_dot);

  /// @brief Time derivative of `RotY(angle)`; see `RotXdot`.
  Mat3 RotYdot(Real angle, Real angle_dot);

  /// @brief Time derivative of `RotZ(angle)`; see `RotXdot`. Used by
  /// `RotSideralMotionDot` to compute `d/dt Rz(theta_ERA)` for the Earth's sidereal
  /// rotation rate, and by the lunar PA/CI Euler-angle rotation-rate reconstruction in
  /// `frame_conversions.cc`.
  Mat3 RotZdot(Real angle, Real angle_dot);

  /// @brief Build the 3x3 skew-symmetric ("cross-product") matrix `[x]_x` of a vector `x`,
  /// such that `[x]_x * v == x.cross(v)` for any vector `v`.
  ///
  /// Exposed to Python via `lupnt.skew` for use in attitude-dynamics / angular-velocity
  /// computations (e.g. forming `omega x r` as a matrix-vector product).
  ///
  /// @param x  Input 3-vector
  /// @return   3x3 skew-symmetric matrix of `x`
  Mat3 Skew(const Vec3& x);

  /// @brief Copy the elements of an Eigen dense vector/array into a `std::vector<T>`,
  /// converting element type via `static_cast`.
  ///
  /// @param x  Input Eigen dense vector/array
  /// @return   `std::vector<T>` with the same elements, cast to `T`
  template <typename T, typename Derived>
  std::vector<T> EigenToVector(const Eigen::DenseBase<Derived>& x) {
    std::vector<T> y(x.size());
    for (int i = 0; i < x.size(); i++) y[i] = static_cast<T>(x(i));
    return y;
  }

  /// @brief Copy the elements of a fixed-size Eigen column vector into a `std::array<T,
  /// N>`, converting element type via `static_cast`.
  ///
  /// @param x  Input fixed-size Eigen column vector (size N)
  /// @return   `std::array<T, N>` with the same elements, cast to `T`
  template <typename T, int N, typename Derived>
  std::array<T, N> EigenToArray(const Eigen::Matrix<T, N, 1>& x) {
    std::array<T, N> y;
    for (int i = 0; i < N; i++) y[i] = static_cast<T>(x(i));
    return y;
  }

  /// @brief Ratio `eta` of the area of the orbital sector swept between two position
  /// vectors to the area of the triangle they form with the focus, for a given
  /// time-of-flight `tau` (Gauss's sector-triangle ratio, solved via Hansen's method).
  ///
  /// Used by `conversions/state_conversions.cc` (Gauss/Lambert-style position-velocity
  /// conversion, e.g. converting two position vectors plus a transfer time into a state
  /// vector) as part of the Herrick-Gibbs/Gauss orbit-determination calculation.
  ///
  /// @param r1   First position vector [any consistent length unit]
  /// @param r2   Second position vector [same unit as `r1`]
  /// @param tau  Time of flight between `r1` and `r2`, scaled per Montenbruck & Eberhard's
  ///             convention [consistent time unit]
  /// @return     Sector-to-triangle area ratio `eta`; throws if Hansen's iteration does
  ///             not converge
  /// @ref
  /// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
  /// applications. Berlin : New York: Springer, 2000.
  /// doi: 10.1007/978-3-642-58351-3.
  Real RatioOfSectorToTriangleArea(const Vec3& r1, const Vec3& r2, Real tau);

  /// @brief Solve the linear least-squares system `A*x = b` via the SVD-based
  /// pseudo-inverse (`Eigen::JacobiSVD`), robust to rank-deficient/ill-conditioned `A`.
  ///
  /// @param A  Coefficient matrix
  /// @param b  Right-hand-side vector
  /// @return   Least-squares (minimum-norm) solution `x`
  VecXd SolveLinearEqSVD(const MatXd& A, const VecXd& b);  // Ax = b

  /// @brief Matrix-right-hand-side overload of `SolveLinearEqSVD`: solves `A*X = B`
  /// column-by-column via the SVD of `A`.
  ///
  /// @param A  Coefficient matrix
  /// @param B  Right-hand-side matrix
  /// @return   Least-squares (minimum-norm) solution `X`
  MatXd SolveLinearEqSVD(const MatXd& A, const MatXd& B);  // AX = B

  /// @brief Moore-Penrose pseudo-inverse of `A` via SVD, with singular values below
  /// `1e-6` treated as zero.
  ///
  /// @param A  Input matrix
  /// @return   Pseudo-inverse of `A`
  MatXd PseudoInverse(const MatXd& A);

  /// @brief Compute `func(x0)` together with its Jacobian `J = d(func)/d(x0)`, evaluating
  /// each column of `J` in parallel via finite differences/autodiff (`jacobian`/`wrt`/`at`
  /// from autodiff).
  ///
  /// Called by `Integrator::Propagate`/`PropagateEx` (the `MatXd* J`/state-transition-
  /// matrix overloads) to compute the sensitivity of a propagated state with respect to
  /// its initial condition, by treating the whole propagation as `func`.
  ///
  /// @param func  Function whose Jacobian is to be evaluated
  /// @param x0    Point at which to evaluate `func` and its Jacobian [size n]
  /// @param J     Output Jacobian matrix, resized to `m x n` where `m = func(x0).size()`
  /// @return      `func(x0)` [size m]
  VecX JacobianParallel(const std::function<VecX(const VecX&)>& func, const VecX& x0, MatXd& J);

  /// @brief Half-vectorization of a square matrix: stack the lower-triangular elements
  /// (including the diagonal) column-wise into a vector.
  ///
  /// @param x  Square input matrix (size n x n)
  /// @return   Vector of length `n*(n+1)/2` containing the lower-triangular elements of
  ///           `x`
  VecX Vech(const MatX& x);

  /// @brief Rotation angle [rad] of a 3x3 rotation matrix `R`, computed from its trace via
  /// `acos((trace(R) - 1) / 2)`.
  ///
  /// @param R  3x3 rotation matrix
  /// @return   Rotation angle [rad] in `[0, pi]`
  Real RotationAngle(const Mat3& R);

}  // namespace lupnt

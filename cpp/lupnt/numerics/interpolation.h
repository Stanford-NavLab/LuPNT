#include <iostream>
#include <vector>

#include "lupnt/core/constants.h"

namespace lupnt {

  /// @brief Piecewise-linear interpolation of a 1D tabulated function `z(x)` at point `ix`.
  ///
  /// Used throughout the simulator for lookups in tabulated data: e.g.
  /// `Sp3Loader::GetPositionEcef` interpolates GNSS satellite ECEF position/clock-bias
  /// columns vs. epoch, `GnssConstellation`/`GnssMeasurement` interpolate fitted ECI state
  /// history vs. epoch, and `Antenna::ComputeGain` interpolates a 1D antenna gain pattern
  /// vs. angle.
  ///
  /// @param x  Strictly monotonic table of independent-variable samples
  /// @param z  Table of dependent-variable values, same size as `x`
  /// @param ix Query point; must lie within `[x(0), x(end)]` (within `EPS`), else throws
  /// @return   Linearly interpolated value of `z` at `ix`
  double LinearInterp1d(const VecXd& x, const VecXd& z, double ix);

  /// @brief Bilinear interpolation of a 2D tabulated function `z(x, y)` at point `(ix, iy)`.
  ///
  /// Used by `Antenna::ComputeGain` to look up a 2D antenna gain pattern as a function of
  /// boresight angle (`phi`) and azimuth (`theta`).
  ///
  /// @param x  Strictly monotonic table of samples along the first axis (rows of `z`)
  /// @param y  Strictly monotonic table of samples along the second axis (columns of `z`)
  /// @param z  Table of values, size `x.size() x y.size()`
  /// @param ix Query point along `x`; must lie within `[x(0), x(end)]` (within `EPS`)
  /// @param iy Query point along `y`; must lie within `[y(0), y(end)]` (within `EPS`)
  /// @return   Bilinearly interpolated value of `z` at `(ix, iy)`
  double LinearInterp2d(const VecXd& x, const VecXd& y, const MatXd& z, double ix, double iy);

  /// @brief Lagrange-polynomial interpolator for evaluating one or more tabulated functions
  /// at a fixed query point `xi`, using a fixed-order local polynomial centered near `xi`.
  ///
  /// Used by the EOP (`GetEopData`, lupnt/data/eop.cc) and IAU SOFA constants
  /// (`GetIauSofaData`, lupnt/data/iau_sofa.cc) loaders to interpolate tabulated
  /// Earth-orientation/precession-nutation series to an arbitrary epoch: constructing the
  /// interpolator once selects the `order`+1 table points nearest `xi` and precomputes the
  /// Lagrange weights, and `Interpolate` is then called once per data column (e.g.
  /// x_pole, y_pole, UT1-UTC, ...) sharing the same abscissa `x`.
  class LagrangeInterpolator {
  public:
    /// @brief Construct an interpolator for query point `xi` over abscissa table `x`,
    /// using a local Lagrange polynomial of degree `order` (i.e. `order` table points).
    ///
    /// Selects the `order` consecutive points of `x` whose span is best centered on `xi`
    /// (via `ComputeFirstIndex`) and precomputes the corresponding Lagrange basis weights
    /// (via `ComputeWeights`).
    ///
    /// @param x      Strictly monotonic table of abscissa samples (e.g. MJD/JD epochs)
    /// @param xi     Query point; must lie within `[x(0), x(end)]` (within `EPS`), else
    ///                throws
    /// @param order  Number of table points used by the local interpolant (interpolant
    ///                degree = `order` - 1); `x.size()` must be > `order`
    LagrangeInterpolator(const VecXd& x, double xi, int order);

    /// @brief Evaluate the Lagrange interpolant of tabulated values `z` (sharing the
    /// abscissa `x` passed to the constructor) at the query point `xi`.
    ///
    /// @param z  Table of dependent-variable values, same size as `x`
    /// @return   Interpolated value of `z` at `xi`
    double Interpolate(const VecXd& z);

  private:
    VecXd x_;
    double xi_;
    int order_;
    VecXd weights_;
    int i0_;

    /// @brief Find the starting index `i0_` of the `order_` consecutive points of `x_`
    /// whose midpoint is closest to `xi_`.
    ///
    /// Called once by the constructor before `ComputeWeights`.
    void ComputeFirstIndex();

    /// @brief Precompute the Lagrange basis weights for the `order_` points
    /// `x_(i0_) ... x_(i0_ + order_ - 1)` evaluated at `xi_`.
    ///
    /// Called once by the constructor; the resulting `weights_` are reused by every
    /// `Interpolate` call.
    void ComputeWeights();
  };

}  // namespace lupnt

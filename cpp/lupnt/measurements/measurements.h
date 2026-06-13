#pragma once

#include "lupnt/core/definitions.h"

namespace lupnt {

  /// @brief Compute the geometric (true) range between two position vectors.
  ///
  /// A minimal generic two-way/one-way range measurement model, used e.g. as
  /// a `FilterMeasurementFunction` building block in EKF/UKF examples
  /// (`example_adaptive.cc`) when a full GnssMeasurement is not needed.
  ///
  /// @param r1 First position [m], any common frame
  /// @param r2 Second position [m], same frame as `r1`
  /// @return   Range `||r2 - r1||` [m] as a length-1 vector
  Vec1 Range(const VecX& r1, const VecX& r2);

  /// @brief Compute the range-rate (radial relative velocity) between two
  /// point masses given their positions and velocities.
  ///
  /// Companion to Range() for generic range-rate/Doppler measurement models
  /// used in filter examples and tests.
  ///
  /// @param r1 First position [m], any common frame
  /// @param r2 Second position [m], same frame as `r1`
  /// @param v1 First velocity [m/s], same frame as `r1`
  /// @param v2 Second velocity [m/s], same frame as `r1`
  /// @return   Range rate `(v2-v1).dot(r2-r1) / ||r2-r1||` [m/s] as a length-1 vector
  Vec1 RangeRate(const VecX& r1, const VecX& r2, const VecX& v1, const VecX& v2);

  /// @brief Overload of RangeRate() taking stacked position/velocity state
  /// vectors `rv = [r; v]` instead of separate position and velocity vectors.
  Vec1 RangeRate(const VecX& rv1, const VecX& rv2);

  /// @brief Compute both the geometric range and range-rate between two point
  /// masses in one call.
  ///
  /// Convenience combination of Range() and RangeRate(), e.g. for a
  /// range+range-rate `FilterMeasurementFunction` measurement model.
  ///
  /// @param r1 First position [m], any common frame
  /// @param r2 Second position [m], same frame as `r1`
  /// @param v1 First velocity [m/s], same frame as `r1`
  /// @param v2 Second velocity [m/s], same frame as `r1`
  /// @return   [range [m], range rate [m/s]]
  Vec2 RangeAndRangeRate(const VecX& r1, const VecX& r2, const VecX& v1, const VecX& v2);

  /// @brief Overload of RangeAndRangeRate() taking stacked position/velocity
  /// state vectors `rv = [r; v]` instead of separate position and velocity vectors.
  Vec2 RangeAndRangeRate(const VecX& rv1, const VecX& rv2);

  /// @brief Compute a pseudorange: the geometric range biased by the
  /// transmitter/receiver clock offsets, in light-time-equivalent distance.
  ///
  /// Generic pseudorange measurement model `Range(r1,r2) + C*(b1-b2)`, used as
  /// a lightweight alternative to GnssMeasurement::ComputeValue's pseudorange
  /// term in simplified filter examples/tests.
  ///
  /// @param r1 First (e.g. receiver) position [m], any common frame
  /// @param r2 Second (e.g. transmitter) position [m], same frame as `r1`
  /// @param b1 Clock bias of point 1 [s]
  /// @param b2 Clock bias of point 2 [s]
  /// @return   Pseudorange [m] as a length-1 vector
  Vec1 Pseudorange(const VecX& r1, const VecX& r2, Real b1, Real b2);

  /// @brief Compute a pseudorange-rate: the range-rate biased by the
  /// transmitter/receiver clock drifts, in light-time-equivalent velocity.
  ///
  /// Generic Doppler/pseudorange-rate measurement model
  /// `RangeRate(...) + C*(d1-d2)`, the rate analog of Pseudorange().
  ///
  /// @param r1 First position [m], any common frame
  /// @param r2 Second position [m], same frame as `r1`
  /// @param v1 First velocity [m/s], same frame as `r1`
  /// @param v2 Second velocity [m/s], same frame as `r1`
  /// @param d1 Clock drift of point 1 [s/s]
  /// @param d2 Clock drift of point 2 [s/s]
  /// @return   Pseudorange rate [m/s] as a length-1 vector
  Vec1 PseudorangeRate(const VecX& r1, const VecX& r2, const VecX& v1, const VecX& v2, Real d1,
                       Real d2);

  /// @brief Overload of PseudorangeRate() taking stacked position/velocity
  /// state vectors `rv = [r; v]` instead of separate position and velocity vectors.
  Vec1 PseudorangeRate(const VecX& rv1, const VecX& rv2, Real d1, Real d2);
}  // namespace lupnt

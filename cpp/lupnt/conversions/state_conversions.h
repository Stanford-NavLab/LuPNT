#pragma once

#include "lupnt/core/definitions.h"
#include "lupnt/numerics/vector_macros.h"
#include "lupnt/states/state.h"
#include "lupnt/states/tle.h"

namespace lupnt {

  /**
   * @brief Convert a Cartesian position-velocity state to classical orbital elements.
   * @param rv Cartesian state in the same distance/time units as `GM`.
   * @param GM Gravitational parameter consistent with the state units.
   * @return Classical orbital element state `[a, e, i, RAAN, arg_periapsis, anomaly]`.
   */
  State CartToClassical(const State& rv, Real GM);

  /**
   * @brief Estimate classical orbital elements from two Cartesian position states.
   * @param dt Time elapsed between `r1` and `r2`.
   * @param r1 First position state.
   * @param r2 Second position state.
   * @param GM Gravitational parameter consistent with the state units.
   * @return Classical orbital element state.
   */
  State CartToClassical(Real dt, const State& r1, const State& r2, Real GM);

  /**
   * @brief Convert an inertial deputy state to a synodic state relative to a chief state.
   * @param rv_c Chief inertial Cartesian state.
   * @param rv_d Deputy inertial Cartesian state.
   * @return Deputy state in the chief-centered synodic frame.
   */
  State InertialToSynodic(const State& rv_c, const State& rv_d);

  /**
   * @brief Convert a synodic deputy state back to inertial coordinates.
   * @param rv_c Chief inertial Cartesian state.
   * @param rv_syn_d Deputy state in the chief-centered synodic frame.
   * @return Deputy inertial Cartesian state.
   */
  State SynodicToInertial(const State& rv_c, const State& rv_syn_d);

  /**
   * @brief Convert classical orbital elements to Cartesian position and velocity.
   * @param coe Classical orbital element state `[a, e, i, RAAN, arg_periapsis, anomaly]`.
   * @param GM Gravitational parameter consistent with the semi-major-axis units.
   * @return Cartesian state in the same frame as `coe`.
   */
  State ClassicalToCart(const State& coe, Real GM);

  /// @brief Convert classical orbital elements to quasi-nonsingular elements.
  State ClassicalToQuasiNonsing(const State& coe, Real GM = 0);

  /// @brief Convert classical orbital elements to equinoctial elements.
  State ClassicalToEquinoctial(const State& coe, Real GM = 0);

  /// @brief Convert classical orbital elements to Delaunay elements.
  State ClassicalToDelaunay(const State& coe, Real GM);

  /// @brief Convert quasi-nonsingular elements to classical orbital elements.
  State QuasiNonsingToClassical(const State& qnsoeVec, Real GM = 0);

  /// @brief Convert equinoctial elements to classical orbital elements.
  State EquinoctialToClassical(const State& eqoe, Real GM = 0);

  /// @brief Convert Delaunay elements to classical orbital elements.
  State DelaunayToClassical(const State& deloe, Real GM);

  /// @brief Convert relative quasi-nonsingular elements to classical orbital elements.
  State RelQuasiNonsingToClassical(const State& coe, const State& rqnsoe);

  /// @brief Convert a two-line element set to classical orbital elements.
  State TleToClassical(const TLE& tle, Real GM);

  // Vector definitions
  VEC_DEF_VECTOR_REAL(ClassicalToCart, 6);
  VEC_DEF_VECTOR_VECTOR(InertialToSynodic, 6);
  VEC_DEF_VECTOR_VECTOR(SynodicToInertial, 6);
  VEC_DEF_VECTOR_REAL(CartToClassical, 6);
  VEC_DEF_VECTOR_REAL(ClassicalToQuasiNonsing, 6);
  VEC_DEF_VECTOR_REAL(ClassicalToEquinoctial, 6);
  VEC_DEF_VECTOR_REAL(ClassicalToDelaunay, 6);
  VEC_DEF_VECTOR_REAL(QuasiNonsingToClassical, 6);
  VEC_DEF_VECTOR_REAL(EquinoctialToClassical, 6);
  VEC_DEF_VECTOR_REAL(DelaunayToClassical, 6);
  VEC_DEF_VECTOR_VECTOR(RelQuasiNonsingToClassical, 6);

}  // namespace lupnt

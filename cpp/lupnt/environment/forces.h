/**
 * @file gravity_field.h
 * @author Stanford NAVLAB
 * @brief Gravity field model and acceleration calculation
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

// C++ includes
#include <string>
#include <vector>

// lupnt includes
#include "lupnt/core/constants.h"
#include "lupnt/core/file.h"
#include "lupnt/environment/body.h"
#include "lupnt/numerics/math_utils.h"

template <typename T> using Vector3 = Eigen::Vector3<T>;
template <typename T> using Matrix3 = Eigen::Matrix<T, 3, 3>;
template <typename T> using MatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

namespace lupnt {

  /// @brief Compute the gravitational acceleration of a spherical-harmonic
  /// gravity field at a body-fixed position, via Montenbruck & Gill's harmonic
  /// (Pines/Cunningham) recursion.
  ///
  /// Called once per integration step by `NumericalOrbitDynamics::CalcContrib`
  /// (templated on `Real` for autodiff or `double` for fast numerical
  /// propagation) to add the central body's high-fidelity gravity contribution
  /// -- e.g. an EGM96 Earth field or GRGM900C lunar field loaded via
  /// `ReadHarmonicGravityField` -- once `r` has been rotated into the body-fixed
  /// frame (e.g. ITRF or MOON_PA).
  ///
  /// @tparam T     Numeric type (`double` for fast propagation, `Real` for autodiff)
  /// @param r      Position vector in the body-fixed frame [m]
  /// @param GM     Gravitational parameter of the central body [m^3/s^2]
  /// @param R_ref  Reference radius of the gravity field model [m]
  /// @param CS     Unnormalized spherical-harmonic coefficients (C on/above diagonal, S below),
  /// from `GravityField::CS`
  /// @param n_max  Maximum degree of the spherical-harmonic expansion to evaluate
  /// @param m_max  Maximum order of the spherical-harmonic expansion to evaluate (m_max <= n_max;
  /// m_max = 0 for zonals only)
  /// @return       Gravitational acceleration in the body-fixed frame [m/s^2]
  template <typename T> Vector3<T> AccelarationGravityField(const Vector3<T>& r, T GM, T R_ref,
                                                            const MatrixX<T>& CS, int n_max,
                                                            int m_max);

  /// @brief Compute the third-body (point-mass) gravitational acceleration on a
  /// spacecraft, including the indirect (central-body recoil) term.
  ///
  /// Called by `NumericalOrbitDynamics::CalcContrib` for each perturbing body
  /// that does NOT use a spherical-harmonic gravity field (`body.use_gravity_field
  /// == false`), e.g. the Sun, Moon, or outer planets perturbing an Earth- or
  /// Moon-centered orbit. Implements `a = -GM * (d/|d|^3 + s/|s|^3)` where `d = r
  /// - s`, i.e. the direct attraction on the spacecraft minus the central body's
  /// own acceleration toward the same point mass.
  ///
  /// @param r  Spacecraft position relative to the central body [m]
  /// @param s  Perturbing point-mass position relative to the central body [m]
  /// @param GM Gravitational parameter of the perturbing point mass [m^3/s^2]
  /// @return   Perturbing acceleration on the spacecraft [m/s^2]
  Vec3 AccelerationPointMass(const Vec3& r, const Vec3& s, Real GM);

  /// @brief Compute the cannonball solar radiation pressure (SRP) acceleration on
  /// a spacecraft.
  ///
  /// Called by `NumericalOrbitDynamics::CalcContrib` once per step (when SRP is
  /// enabled), scaled by the eclipse `Illumination` factor for the relevant
  /// occulting body, to add the SRP contribution to the total acceleration.
  ///
  /// @param r          Spacecraft position relative to the Sun-illuminated central body [m]
  /// @param r_sun      Sun position relative to the same central body [m]
  /// @param bcoeff_srp Spacecraft SRP ballistic coefficient, `CR * area / mass` [m^2/kg]
  /// @param P0         Solar radiation pressure at 1 AU [N/m^2]
  /// @param AU         Length of one Astronomical Unit [m]
  /// @return           SRP acceleration on the spacecraft [m/s^2]
  Vec3 AccelerationSolarRadiation(const Vec3& r, const Vec3& r_sun, Real bcoeff_srp, Real P0,
                                  Real AU);

  /// @brief Compute the first-order post-Newtonian (Schwarzschild) relativistic
  /// acceleration correction relative to a central body.
  ///
  /// Called by `NumericalOrbitDynamics::CalcContrib` (when relativistic
  /// corrections are enabled) once for the Sun-relative term and again for an
  /// optional "relativity center" body (e.g. Earth or Moon), summing both into
  /// the total acceleration. Implements Montenbruck & Gill Sec. 3.7.3 Eq. 3.146.
  ///
  /// @param r       Spacecraft position relative to the central body [m]
  /// @param v       Spacecraft velocity relative to the central body [m/s]
  /// @param GM      Gravitational parameter of the central body [m^3/s^2]
  /// @param c_light Speed of light [m/s]
  /// @return        Post-Newtonian acceleration correction [m/s^2]
  Vec3 AccelerationRelativisticCorrection(const Vec3& r, const Vec3& v, Real GM, Real c_light = C);

  /**
   * @brief Computes the acceleration due to atmospheric drag.
   *
   * @param mjd_tt Modified Julian Date (Terrestrial Time).
   * @param rv State vector (position and velocity) of the spacecraft.
   * @param T Transformation matrix from inertial to body frame.
   * @param bcoeff_drag Ballistic coefficient for atmospheric drag.
   * @return Vec3 Acceleration vector due to atmospheric drag.
   */
  Vec3 AccelerationDrag(Real mjd_tt, const Vec6& rv, const Mat3& T, Real bcoeff_drag);

  /**
   * @brief Computes the solar-radiation shadow function using apparent disk overlap.
   *
   * @param r Position vector of the spacecraft relative to the occulting body.
   * @param r_sun Position vector of the Sun relative to the occulting body.
   * @param R_body Radius of the occulting body.
   * @param R_sun Radius of the Sun.
   * @return Real Shadow function nu: 0 in umbra, 1 in sunlight, and between 0 and 1 in penumbra.
   */
  Real ShadowFunction(const Vec3& r, const Vec3& r_sun, Real R_body, Real R_sun = R_SUN);

  /**
   * @brief Computes the illumination factor for solar radiation pressure.
   *
   * @param r Position vector of the spacecraft relative to the occulting body.
   * @param r_sun Position vector of the Sun relative to the occulting body.
   * @param R_body Radius of the occulting body.
   * @return Real Illumination factor (0 to 1).
   */
  Real Illumination(const Vec3& r, const Vec3& r_sun, Real R_body);

  /**
   * @brief Computes the atmospheric density using the Harris-Priester model.
   *
   * @param mjd_tt Modified Julian Date (Terrestrial Time).
   * @param r_tod Position vector in the true of date frame.
   * @return Real Atmospheric density.
   */
  Real DensityHarrisPriester(Real mjd_tt, const Vec3& r_tod);

  /**
   * @brief Computes the total acceleration on a spacecraft due to Earth's gravity, solar radiation
   * pressure, and atmospheric drag.
   *
   * @param mjd_tt Modified Julian Date (Terrestrial Time).
   * @param rv State vector (position and velocity) of the spacecraft.
   * @param bcoeff_srp Ballistic coefficient for solar radiation pressure.
   * @param bcoeff_drag Ballistic coefficient for atmospheric drag.
   * @param grav Gravity field model.
   * @return Vec3 Total acceleration vector on the spacecraft.
   */
  Vec3 AccelerationEarthSpacecraft(Real mjd_tt, const Vec6& rv, Real bcoeff_srp, Real bcoeff_drag,
                                   GravityField<Real> grav);

}  // namespace lupnt

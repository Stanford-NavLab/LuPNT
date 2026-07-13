#pragma once

#include <utility>
#include <vector>

#include "lupnt/core/constants.h"
#include "lupnt/core/definitions.h"

namespace lupnt {

  /// @file ground_station_corrections.h
  /// @brief Signal-path and station-location corrections for Earth ground-station tracking
  /// of a distant (e.g. lunar) spacecraft.
  ///
  /// These model the physical delays and displacements that a real Deep Space Network
  /// station experiences but that a purely geometric range model omits. They are meant to be
  /// applied to the *truth* observable a station sensor generates, so that an estimator which
  /// does not model them sees them as realistic measurement errors.

  /// @brief Tropospheric slant delay via the Saastamoinen zenith model with a `1/sin(el)`
  /// obliquity mapping.
  ///
  /// The zenith delay is the sum of a hydrostatic ("dry") and a wet component; both scale to
  /// the slant path by the elevation-dependent obliquity factor. Water-vapour partial pressure
  /// is derived from relative humidity via a Tetens saturation formula.
  ///
  /// @param elevation_rad Topocentric elevation of the spacecraft [rad] (must be > 0).
  /// @param latitude_rad  Geodetic latitude of the station [rad].
  /// @param height_m      Station height above the ellipsoid [m].
  /// @param pressure_hPa  Total surface pressure [hPa].
  /// @param temperature_K Surface temperature [K].
  /// @param humidity_pct  Relative humidity [%] in [0, 100].
  /// @return One-way excess range delay [m] (>= 0).
  double TroposphereDelaySaastamoinen(double elevation_rad, double latitude_rad, double height_m,
                                      double pressure_hPa, double temperature_K,
                                      double humidity_pct);

  /// @brief Surface pressure of the ISO standard atmosphere at a given height [hPa].
  double StandardAtmospherePressureHPa(double height_m);

  /// @brief Surface temperature of the ISO standard atmosphere at a given height [K].
  double StandardAtmosphereTemperatureK(double height_m);

  /// @brief Ionospheric slant delay from a single-layer (thin-shell) vertical-TEC model.
  ///
  /// The slant TEC is obtained by mapping the vertical TEC through the ionospheric pierce
  /// point of a thin shell at `shell_height_m`, and the group delay follows the dispersive
  /// `40.3 * STEC / f^2` relation. The delay is positive (a group-path lengthening) and
  /// scales as `1/f^2`, so a dual-frequency link can cancel it.
  ///
  /// @param elevation_rad   Topocentric elevation of the spacecraft [rad].
  /// @param vtec_tecu        Vertical total electron content [TECU] (1 TECU = 1e16 el/m^2).
  /// @param frequency_hz     Carrier frequency [Hz].
  /// @param station_height_m Station height above the ellipsoid [m].
  /// @param shell_height_m   Ionospheric thin-shell height [m].
  /// @param earth_radius_m   Earth radius used for the pierce-point geometry [m].
  /// @return One-way group delay [m] (>= 0).
  double IonosphereDelayThinShell(double elevation_rad, double vtec_tecu, double frequency_hz,
                                  double station_height_m, double shell_height_m = 350.0e3,
                                  double earth_radius_m = R_EARTH);

  /// @brief Relativistic (Shapiro) one-way range delay for a signal travelling between
  /// `r_tx` and `r_rx`, summed over the given gravitating bodies.
  ///
  /// Each body contributes `2 GM / c^2 * ln[(r1 + r2 + r12)/(r1 + r2 - r12)]`, with `r1`, `r2`
  /// the transmitter/receiver distances to the body and `r12` the transmitter-receiver range.
  /// For an Earth-station to lunar-spacecraft link the Earth term dominates and the Sun term
  /// is a secondary contribution.
  ///
  /// @param r_tx    Transmitter position [m], in a common inertial frame.
  /// @param r_rx    Receiver position [m], same frame.
  /// @param bodies  List of `(body position [m], GM [m^3/s^2])` in the same frame.
  /// @return One-way range delay [m] (>= 0 for the usual geometry).
  double ShapiroRangeDelay(const Vec3d& r_tx, const Vec3d& r_rx,
                           const std::vector<std::pair<Vec3d, double>>& bodies);

  /// @brief Degree-2 in-phase solid Earth tide station displacement (IERS Conventions form).
  ///
  /// The crust displacement is `sum_j (GM_j/GM_E)(R_E^4/R_j^3){ h2 r_hat[(3(R_hat_j.r_hat)^2 -
  /// 1)/2] + 3 l2 (R_hat_j.r_hat)[R_hat_j - (R_hat_j.r_hat) r_hat] }` summed over the
  /// tide-raising bodies (Moon and Sun). Because it is built from frame-invariant dot products
  /// of unit vectors, evaluating it with all vectors expressed in an inertial frame yields the
  /// displacement directly in that frame.
  ///
  /// @param r_station   Geocentric station position [m].
  /// @param tide_bodies List of `(geocentric body position [m], GM [m^3/s^2])` (Moon, Sun).
  /// @param GM_earth    Earth gravitational parameter [m^3/s^2].
  /// @param earth_radius_m Earth equatorial radius [m].
  /// @param h2          Degree-2 vertical (radial) Love number.
  /// @param l2          Degree-2 horizontal Shida number.
  /// @return Station displacement vector [m], same frame as the inputs.
  Vec3d SolidEarthTideDisplacement(const Vec3d& r_station,
                                   const std::vector<std::pair<Vec3d, double>>& tide_bodies,
                                   double GM_earth, double earth_radius_m, double h2 = 0.6078,
                                   double l2 = 0.0847);

}  // namespace lupnt

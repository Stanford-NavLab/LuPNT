/**
 * @file atmosphere.h
 * @author Stanford NAV Lab
 * @brief  Atmosphere and ionosphere models
 * @version 0.1
 * @date 2024-12-04
 *
 * @copyright Copyright (c) 2024
 *
 */

#pragma once

#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

#include "lupnt/core/definitions.h"

namespace lupnt {

  /// @brief Compute the GPS broadcast (Klobuchar) ionospheric slant delay for a
  /// given receiver-satellite geometry and frequency.
  ///
  /// Used by pseudorange/carrier-phase measurement models to apply (or remove)
  /// the single-frequency ionospheric delay correction broadcast in the GPS
  /// navigation message, scaling the standard L1 model to an arbitrary
  /// frequency via the (L1_FREQ / freq_Hz)^2 factor.
  ///
  /// @param t_gps      GPS time of day [s]
  /// @param elevation  Satellite elevation angle as seen from the receiver [rad]
  /// @param azimuth    Satellite azimuth angle as seen from the receiver [rad]
  /// @param latitude_u  Receiver geodetic latitude [rad]
  /// @param longitude_u Receiver geodetic longitude [rad]
  /// @param freq_Hz    Signal frequency to scale the delay to [Hz]
  /// @param alpha      Broadcast ionospheric correction amplitude coefficients (4x1)
  /// @param beta       Broadcast ionospheric correction period coefficients (4x1)
  /// @return           Ionospheric slant delay at `freq_Hz` [s]
  double Klobucher(double t_gps, double elevation, double azimuth, double latitude_u,
                   double longitude_u, double freq_Hz, const Vec4d& alpha, const Vec4d& beta);

}  // namespace lupnt

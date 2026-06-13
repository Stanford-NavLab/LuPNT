/**
 * @file tec/klobuchar.h
 * @author Keidai Iiyama
 * @brief This file contains the interface for the Klobuchar model.
 * @version 0.1
 * @date 2025-02-17
 */

#include "lupnt/environment/plasma/tec/klobuchar.h"

namespace pecsim {

  double Klobucher(double t_gps, double elevation, double azimuth, double latitude_u,
                   double longitude_u, double freq_Hz, const Vec4d& alpha, const Vec4d& beta) {
    // Compute the earth-centered angle
    double psi = 0.0137 / (elevation + 0.11) - 0.022;

    // compute the latitude of the ionospheric pierce point
    double phi_i = latitude_u + psi * std::cos(azimuth);
    if (phi_i > 0.416) {
      phi_i = 0.416;
    } else if (phi_i < -0.416) {
      phi_i = -0.416;
    }

    // compute the longitude of the ionospheric pierce point
    double lambda_i = longitude_u + psi * std::sin(azimuth) / std::cos(phi_i);

    // compute the geometric latitude of the IPP
    double phi_m = phi_i + 0.064 * std::cos(lambda_i - 1.617);

    // compute the local time at the IPP
    double t = 43200 * lambda_i + t_gps;  // sec

    // wrap to be within 0 - 86400
    t = fmod(t, 86400);

    // compute the amplitude of delay
    double A_I = 0.0;
    for (int i = 0; i < 4; i++) {
      A_I += alpha[i] * std::pow(phi_m, i);
    }
    if (A_I < 0.0) A_I = 0.0;

    // compute the period of delay
    double P_I = 0.0;
    for (int i = 0; i < 4; i++) {
      P_I += beta[i] * std::pow(phi_m, i);
    }
    if (P_I < 72000.0) P_I = 72000.0;

    // compute the phase of the delay
    double X_I = 2 * M_PI * (t - 50400) / P_I;

    // comute the slant factor
    double F = 1.0 + 16.0 * std::pow(0.53 - elevation, 3);

    // compute the ionospheric delay at L1
    double iono_delay = 0.0;
    if (std::abs(X_I) < 1.57) {
      iono_delay = F * (5.0e-9 + A_I * (1 - std::pow(X_I, 2) / 2 + std::pow(X_I, 4) / 24));
    } else {
      iono_delay = F * 5.0e-9;
    }

    // convert to the frequency of interest
    double L1_FREQ = 1575.42e6;  // L1 frequency in Hz
    iono_delay = iono_delay * std::pow(L1_FREQ / freq_Hz, 2);
    return iono_delay;
  }

}  // namespace pecsim

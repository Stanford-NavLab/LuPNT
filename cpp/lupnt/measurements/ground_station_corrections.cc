#include "lupnt/measurements/ground_station_corrections.h"

#include <algorithm>
#include <cmath>

namespace lupnt {

  double StandardAtmospherePressureHPa(double height_m) {
    // ISO standard atmosphere troposphere (0-11 km): P = P0 (1 - L h / T0)^(g M / R L).
    constexpr double P0 = 1013.25;      // [hPa]
    constexpr double kExp = 5.2558797;  // g0 M / (R* L)
    constexpr double kLapse = 2.25577e-5;
    double base = 1.0 - kLapse * height_m;
    if (base <= 0.0) return 0.0;
    return P0 * std::pow(base, kExp);
  }

  double StandardAtmosphereTemperatureK(double height_m) {
    constexpr double T0 = 288.15;  // [K]
    constexpr double L = 6.5e-3;   // [K/m]
    return T0 - L * height_m;
  }

  double TroposphereDelaySaastamoinen(double elevation_rad, double latitude_rad, double height_m,
                                      double pressure_hPa, double temperature_K,
                                      double humidity_pct) {
    double sin_el = std::sin(elevation_rad);
    if (sin_el <= 1e-3) sin_el = 1e-3;  // guard against the horizon singularity

    // Latitude/height gravity correction factor shared by both components.
    double f_lat_h = 1.0 - 0.00266 * std::cos(2.0 * latitude_rad) - 0.00000028 * height_m;
    if (f_lat_h <= 0.0) f_lat_h = 1.0;

    // Water-vapour partial pressure from relative humidity (Tetens saturation formula).
    double humidity = std::clamp(humidity_pct, 0.0, 100.0);
    double t_c = temperature_K - 273.15;
    double e_sat = 6.11 * std::exp(17.27 * t_c / (t_c + 237.3));  // [hPa]
    double e = (humidity / 100.0) * e_sat;

    double zhd = 0.0022768 * pressure_hPa / f_lat_h;                         // hydrostatic [m]
    double zwd = 0.0022768 * (1255.0 / temperature_K + 0.05) * e / f_lat_h;  // wet [m]
    return (zhd + zwd) / sin_el;
  }

  double IonosphereDelayThinShell(double elevation_rad, double vtec_tecu, double frequency_hz,
                                  double station_height_m, double shell_height_m,
                                  double earth_radius_m) {
    if (frequency_hz <= 0.0 || vtec_tecu <= 0.0) return 0.0;
    // Thin-shell obliquity: sin(z') = (Re + hs)/(Re + hshell) * sin(z), z = zenith angle.
    double z = M_PI / 2.0 - elevation_rad;
    double ratio = (earth_radius_m + station_height_m) / (earth_radius_m + shell_height_m);
    double sin_zp = ratio * std::sin(z);
    sin_zp = std::clamp(sin_zp, -1.0, 1.0);
    double cos_zp = std::sqrt(std::max(0.0, 1.0 - sin_zp * sin_zp));
    if (cos_zp <= 1e-6) cos_zp = 1e-6;
    double mapping = 1.0 / cos_zp;
    double stec_el_m2 = vtec_tecu * 1.0e16 * mapping;  // slant TEC [electrons/m^2]
    return 40.3 * stec_el_m2 / (frequency_hz * frequency_hz);
  }

  double ShapiroRangeDelay(const Vec3d& r_tx, const Vec3d& r_rx,
                           const std::vector<std::pair<Vec3d, double>>& bodies) {
    double r12 = (r_rx - r_tx).norm();
    double delay = 0.0;
    for (const auto& [r_body, gm] : bodies) {
      double r1 = (r_tx - r_body).norm();
      double r2 = (r_rx - r_body).norm();
      double num = r1 + r2 + r12;
      double den = r1 + r2 - r12;
      if (den <= 0.0 || num <= 0.0) continue;
      delay += 2.0 * gm / (C * C) * std::log(num / den);
    }
    return delay;
  }

  Vec3d SolidEarthTideDisplacement(const Vec3d& r_station,
                                   const std::vector<std::pair<Vec3d, double>>& tide_bodies,
                                   double GM_earth, double earth_radius_m, double h2, double l2) {
    double r_sta = r_station.norm();
    if (r_sta <= 0.0 || GM_earth <= 0.0) return Vec3d::Zero();
    Vec3d r_hat = r_station / r_sta;
    double re4 = std::pow(earth_radius_m, 4);

    Vec3d disp = Vec3d::Zero();
    for (const auto& [r_body, gm] : tide_bodies) {
      double rj = r_body.norm();
      if (rj <= 0.0) continue;
      Vec3d rj_hat = r_body / rj;
      double c = rj_hat.dot(r_hat);  // cos(angle between body and station)
      double scale = (gm / GM_earth) * re4 / (rj * rj * rj);
      Vec3d radial = h2 * r_hat * (1.5 * c * c - 0.5);
      Vec3d transverse = 3.0 * l2 * c * (rj_hat - c * r_hat);
      disp += scale * (radial + transverse);
    }
    return disp;
  }

}  // namespace lupnt

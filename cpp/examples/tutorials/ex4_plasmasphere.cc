// Example 4: Plasmasphere and Ionosphere
// ------------------------------------------------------------------------------
// C++ counterpart of python/examples/ex4_plasmasphere.ipynb.
//
// Demonstrates LuPNT's plasma/ionosphere module (the `pecsim` backend behind
// `pylupnt.plasma`):
//
//   Part 1  Sample the GCPM v2.4 electron-density model along the equatorial
//           radial direction (a 1-D slice through the plasmasphere).
//   Part 2  Ray-trace a GNSS signal through the ionosphere/plasmasphere and
//           report the total electron content (TEC) and dispersive delay.
//
// References: GCPM v2.4 (Gallagher et al. 2000), IRI-2007 (Bilitza & Reinisch
// 2008), ray-tracing (Haselgrove 1955 / pecsim implementation).

#include <iomanip>
#include <iostream>
#include <vector>

#include "lupnt/environment/plasma/plasma.h"

using namespace pecsim;

int main() {
  // --- Epoch and space-weather index ----------------------------------------
  const int year = 2026, month = 1, day = 1, hour = 0;
  const double mjd = gregorian_to_mjd(year, month, day, hour, 0, 0.0);
  const double epoch_j2000 = mjd_to_tj2000(mjd);
  const DateTime dt = mjd_to_datetime(mjd);

  double kp = get_kp_index(dt);
  if (kp < 0) kp = 2.0;  // fall back if the Kp table is unavailable

  set_iri_model(IRIModel::IRI_2007);
  std::cout << "Epoch : " << year << "-" << month << "-" << day << " " << hour << ":00 UTC\n";
  std::cout << "J2000 : " << std::fixed << std::setprecision(0) << epoch_j2000 << " s\n";
  std::cout << "Kp    : " << std::setprecision(2) << kp << "\n";
  std::cout << "RE    : " << RE << " km\n\n";

  // --- Part 1: equatorial electron-density profile ---------------------------
  // gcpm_v24(datetime, r_RE, amlt, alatr, kp) -> [ne_cm3, H+, He+, O+].
  //   r_RE  : geocentric distance in Earth radii
  //   amlt  : magnetic local time [h]; 12 = noon
  //   alatr : magnetic latitude [rad]; 0 = magnetic equator
  const double amlt = 12.0;  // noon meridian
  const double alatr = 0.0;  // magnetic equator
  std::cout << "Equatorial electron density (noon meridian), GCPM v2.4:\n";
  std::cout << "  r [RE]   ne [cm^-3]\n";
  for (double r_RE = 1.5; r_RE <= 6.0 + 1e-9; r_RE += 0.5) {
    std::vector<double> n = gcpm_v24(dt, r_RE, amlt, alatr, kp);
    std::cout << "  " << std::setw(5) << std::setprecision(1) << r_RE << "   " << std::scientific
              << std::setprecision(3) << n[0] << std::fixed << "\n";
  }

  // --- Part 2: ray-trace a GNSS link through the plasmasphere -----------------
  // Illustrative geometry: a horizontal ray skimming the plasmasphere at ~1000 km
  // altitude, from a GNSS-altitude transmitter toward a distant (lunar-range)
  // receiver, so the whole path samples the modelled electron content.
  RayTraceConfig config;
  config.freq_Hz = freq_L1;    // GPS L1
  config.step_size = 50.0;     // ray-trace step [km]
  config.cutoff_r = 6.0 * RE;  // stop integrating beyond 6 RE [km]
  config.gradn_dx = 1.0;       // refractive-index gradient step [km]
  config.integ_method = "Euler";
  config.correction = false;  // geometric straight-ray TEC (no bending correction)
  config.fine_correction = false;
  config.use_fortran_gcpm = false;
  config.kp = kp;

  const double h_km = RE + 1000.0;          // tangent altitude ~1000 km
  const Vec3d pos_tx{-20000.0, h_km, 0.0};  // transmitter [km, ECI]
  const Vec3d pos_rx{380000.0, h_km, 0.0};  // receiver (lunar range) [km, ECI]

  std::cout << "\nRay trace (GPS L1) from " << pos_tx.transpose() << " km\n"
            << "                 to   " << pos_rx.transpose() << " km\n";
  PathProfile pp = trace_ray(epoch_j2000, pos_tx, pos_rx, config, /*debug_prop=*/false,
                             /*debug_corr=*/false);
  std::cout << std::setprecision(3);
  std::cout << "  TEC        : " << pp.tecu << " TECU\n";
  std::cout << "  TEC delay  : " << pp.tec_delay_m << " m\n";
  std::cout << "  Total delay: " << pp.total_delay_m << " m\n";

  return 0;
}

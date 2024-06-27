#include <lupnt/lupnt.h>

#include <iomanip>
#include <iostream>

using namespace lupnt;
using namespace std;

int main() {
  // Time
  const Real mjd0_utc = Calendar2ModJulianDate(1997, 01, 01);

  // Ground station
  const Real lon_gs = 11 * RAD;  // [rad]
  const Real lat_gs = 48 * RAD;  // [rad]
  const Real alt_gs = 0;         // [km]

  const Vec3 r_gs =
      LatLonAlt2Cart(Vec3(lat_gs, lon_gs, alt_gs), R_EARTH, WGS84_F);

  // Spacecraft
  const Real a = 960 + R_EARTH;    // Semimajor axis [m]
  const Real e = 0;                // Eccentricity
  const Real i = 97 * RAD;         // Inclination [rad]
  const Real Omega = 130.7 * RAD;  // RA ascend. node [rad]
  const Real omega = 0 * RAD;      // Argument of latitude [rad]
  const Real M0 = 0 * RAD;         // Mean anomaly at epoch [rad]

  const Vec6 coe0(a, e, i, Omega, omega, M0);  // Keplerian elements

  cout << "Exercise 2-4: Topocentric satellite motion" << endl
       << endl
       << "   Date         UTC           Az         El      Dist" << endl
       << "yyyy/mm/dd  hh:mm:ss.sss     [deg]     [deg]     [km]" << endl;

  for (int minute = 6; minute <= 24; minute++) {
    Real mjd_utc = mjd0_utc + minute / MINUTES_DAY;  // [days]
    Real dt = (mjd_utc - mjd0_utc) * SECS_DAY;       // [s]

    Vec6 coe = KeplerianDynamics::PropagateClassicalOE(coe0, dt, GM_EARTH);
    Vec6 rv_eci = Classical2Cart(coe, GM_EARTH);

    Real theta = GreenwichMeanSiderealTime(mjd_utc);  // Note: it should be UT1
    Mat3 R_eci2ecef = RotZ(theta);
    Vec3 r_sat = R_eci2ecef * rv_eci.segment(0, 3);

    auto [az, el, range] =
        unpack(Cart2AzElRange(r_sat, r_gs, R_EARTH, WGS84_F));

    Real mjd_utc_ms =
        round(mjd_utc * SECS_DAY * 1e3, 3) / (SECS_DAY * 1e3) + EPS;
    auto [year, month, day, hour, min, sec] =
        ModJulianDate2Calendar(mjd_utc_ms);

    cout << FormatDate(year, month, day, hour, min, sec, 3);
    cout << setw(10) << fixed << setprecision(1);
    cout << setw(10) << az * DEG;
    cout << setw(10) << el * DEG;
    cout << setw(10) << range << endl;
  };

  return 0;
}
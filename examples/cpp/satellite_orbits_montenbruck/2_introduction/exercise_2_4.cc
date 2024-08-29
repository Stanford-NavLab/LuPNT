#include <lupnt/lupnt.h>

#include <iomanip>
#include <iostream>

using namespace lupnt;
using namespace std;

// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
// applications. Berlin : New York: Springer, 2000.
// doi: 10.1007/978-3-642-58351-3.
int main() {
  // Time
  const Real mjd0_utc = Gregorian2MJD(1997, 01, 01);

  // Ground station
  const Real lon_gs = 11 * RAD;  // [rad]
  const Real lat_gs = 48 * RAD;  // [rad]
  const Real alt_gs = 0;         // [km]

  const Vec3 r_gs = LatLonAlt2Cart(Vec3(lat_gs, lon_gs, alt_gs), R_EARTH, WGS84_F);

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
       << "   Date        UTC           Az        El       Dist" << endl
       << "yyyy/mm/dd hh:mm:ss.sss     [deg]     [deg]     [km]" << endl;

  for (int minute = 6; minute <= 24; minute++) {
    Real mjd_utc = mjd0_utc + minute / MINS_DAY;  // [days]
    Real dt = (mjd_utc - mjd0_utc) * SECS_DAY;    // [s]

    KeplerianDynamics dyn(GM_EARTH);
    Vec6 coe = dyn.PropagateClassicalOE(coe0, 0, dt);
    Vec6 rv_eci = Classical2Cart(coe, GM_EARTH);

    // Note: it should be UT1
    Mat3 R = RotZ(GreenwichMeanSiderealTime(mjd_utc));
    Vec3 r_ecef = R * rv_eci.segment(0, 3);

    Vec3 aer = Cart2AzElRange(r_ecef, r_gs, R_EARTH, WGS84_F);
    auto [az, el, range] = unpack(aer);

    cout << MJD2GregorianString(mjd_utc, 3);
    cout << setw(10) << fixed << setprecision(1);
    cout << setw(10) << az * DEG;
    cout << setw(10) << el * DEG;
    cout << setw(10) << range << endl;
  };

  return 0;
}

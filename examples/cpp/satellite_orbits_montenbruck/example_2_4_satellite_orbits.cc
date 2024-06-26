#include <lupnt/lupnt.h>

#include <iomanip>
#include <iostream>

using namespace lupnt;
using std::cout, std::endl, std::fixed, std::setprecision, std::setw,
    std::setfill;

int main() {
  // Load EOP data
  std::shared_ptr<EOPData> eop_data =
      LoadEOPData(GetFilePath("eopc04_08.62-now"));

  // Time
  real mjd0_utc = DateToModifiedJulianDate(1997, 1, 1, 0, 0, 0);  // Epoch

  // Ground station
  real lon_gs = 11 * RAD_PER_DEG;  // [rad]
  real lat_gs = 48 * RAD_PER_DEG;  // [rad]
  real alt_gs = 0;                 // [m]
  Vector3 r_gs = GeodeticToCartesian(Vector3(lat_gs, lon_gs, alt_gs), R_EARTH,
                                     FLATTENING_EARTH_WGS84);

  // Spacecraft
  real a = 960 + R_EARTH;                   // Semimajor axis [Km]
  real e = 0;                               // Eccentricity
  real i = 97 * RAD_PER_DEG;                // Inclination [rad]
  real Omega = 130.7 * RAD_PER_DEG;         // RA ascend. node [rad]
  real omega = 0 * RAD_PER_DEG;             // Argument of latitude [rad]
  real M0 = 0 * RAD_PER_DEG;                // Mean anomaly at epoch [rad]
  Vector6 coe0(a, e, i, Omega, omega, M0);  // Classical orbital elements

  // Header
  int w = 10;
  cout << "Exercise 2-4 (Topocentric satellite motion)\n\n";
  cout << setw(w) << "UTC" << setw(w) << "az" << setw(w) << "el" << setw(w)
       << "rho" << endl;
  cout << setw(w) << "hh:mm:ss.s" << setw(w) << "[deg]" << setw(w) << "[deg]"
       << setw(w) << "[km]" << endl;

  for (int minute = 6; minute <= 24; ++minute) {
    // Time
    real mjd_utc = mjd0_utc + minute * SECS_PER_MINUTE / SECS_PER_DAY;
    real mjd_ut1 = UTCtoUT1(mjd_utc);
    real t = (mjd_utc - mjd0_utc) * SECS_PER_DAY;

    Vector6 coe = KeplerianDynamics::PropagateClassicalOE(coe0, t, GM_EARTH);
    Vector6 rv_eci = ClassicalToCartesian(coe, GM_EARTH);
    real theta_era = GreenwichMeanSiderealTime(mjd_ut1);
    Matrix3 eci_to_ecef = Rot3(theta_era);
    Vector3 r_ecef = eci_to_ecef * rv_eci.head(3);
    auto [az, el, rho] = unpack(
        CartesianToAzimuthElevationRange(r_gs, r_ecef, FLATTENING_EARTH_WGS84));

    auto [year, month, day, hour, min, sec] = ModifiedJulianDateToDate(mjd_utc);

    cout << setw(2) << setfill('0') << fixed << setprecision(3) << hour << ":";
    cout << setw(2) << setfill('0') << fixed << setprecision(3) << min << ":";
    cout << setw(4) << setfill('0') << fixed << setprecision(1) << sec;
    cout << setfill(' ');
    cout << setw(w) << fixed << setprecision(1) << az * DEG_PER_RAD;
    cout << setw(w) << fixed << setprecision(1) << el * DEG_PER_RAD;
    cout << setw(w) << fixed << setprecision(1) << rho;
    cout << "\n";
  }

  return 0;
}
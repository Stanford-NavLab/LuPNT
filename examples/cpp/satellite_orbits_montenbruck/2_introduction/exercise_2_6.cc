#include <lupnt/lupnt.h>

using namespace lupnt;
using namespace std;

// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
// applications. Berlin : New York: Springer, 2000.
// doi: 10.1007/978-3-642-58351-3.
int main() {
  // Ground station
  const Vec3 r_gs(1344.143, 6068.601, 1429.311);            // (x, y, z) [km]
  const Vec3 lla = Cart2LatLonAlt(r_gs, R_EARTH, WGS84_F);  // (lat, lon, alt) [rad, rad, km]

  // Measurements
  struct Measurement {
    Real mjd_utc;
    Real az, el, range;
  };
  const Measurement meas[2]
      = {{Gregorian2MJD(1999, 04, 02, 00, 30, 00.0), 132.67 * RAD, 32.44 * RAD, 16.945450e3},
         {Gregorian2MJD(1999, 04, 02, 03, 00, 00.0), 123.08 * RAD, 50.06 * RAD, 37.350340e3}};

  // Convert observations
  Vec3 r_sat[2];
  for (int i = 0; i < 2; i++) {
    // Earth rotation
    Mat3 R = RotZ(GreenwichMeanSiderealTime(meas[i].mjd_utc));
    Vec3 aer{meas[i].az, meas[i].el, meas[i].range};
    Vec3 r_ecef = AzElRange2Cart(aer, r_gs, R_EARTH, WGS84_F);
    r_sat[i] = R.transpose() * r_ecef;
  }

  // Orbital elements
  Real dt = (meas[1].mjd_utc - meas[0].mjd_utc) * SECS_DAY;
  Vec6 coe = Cart2Classical(dt, r_sat[0], r_sat[1], GM_EARTH);

  int date_prec = 3;
  cout << "Exercise 2-6: Initial orbit determination" << endl << endl;
  cout << "Inertial positions:" << endl << endl;
  cout << setw(36) << "[km]" << setw(12) << "[km]" << setw(12) << "[km]" << endl;
  for (int i = 0; i < 2; i++) {
    cout << "  " << MJD2GregorianString(meas[i].mjd_utc, date_prec) << "  ";
    for (int j = 0; j < 3; j++) {
      cout << setw(12) << r_sat[i](j);
    };
    cout << endl;
  };
  cout << endl;

  cout << "Orbital elements:" << endl << endl;
  cout << "  Epoch (1st obs.)  " << MJD2GregorianString(meas[0].mjd_utc, date_prec) << endl;
  cout << fixed << setprecision(3);
  cout << "  Semimajor axis   " << setw(10) << coe(0) / 1000.0 << " km" << endl;
  cout << setprecision(7);
  cout << "  Eccentricity     " << setw(10) << coe(1) << endl;
  cout << setprecision(3);
  cout << "  Inclination      " << setw(10) << coe(2) * DEG << " deg" << endl;
  cout << "  RA ascend. node  " << setw(10) << coe(3) * DEG << " deg" << endl;
  cout << "  Arg. of perigee  " << setw(10) << coe(4) * DEG << " deg" << endl;
  cout << "  Mean anomaly     " << setw(10) << coe(5) * DEG << " deg" << endl;
  cout << endl;

  return 0;
}

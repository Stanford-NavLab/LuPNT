#include <lupnt/lupnt.h>

#include <iostream>
#include <string>

using namespace lupnt;
using namespace std;

// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
// applications. Berlin : New York: Springer, 2000.
// doi: 10.1007/978-3-642-58351-3.
int main() {
  // Constants
  const int N_step = 8;
  const double dt = 0.5;  // [days]
  const auto fmt = Eigen::IOFormat(9, 0, " ", " ", "", "", "", "");

  cout << "Exercise 3-2: Lunar Ephemerides " << endl << endl;
  cout << "  Moon position from low precision analytical theory" << endl;
  cout << "  Date [TT]                 " << " Position [km] " << endl;

  Real mjd0 = Gregorian2MJD(2006, 03, 14, 00, 00, 0.0);
  Real mjd_tt;
  Vec3 r;
  for (int i = 0; i <= N_step; i++) {
    mjd_tt = mjd0 + i * dt;
    r = MoonPositionLowPrecision(mjd_tt) / 1000.0;
    cout << "  " << MJD2GregorianString(mjd_tt, 1) << "      " << r.transpose().format(fmt) << endl;
  };

  cout << endl << " Moon position from DE405" << endl;
  cout << " Date [TT]                 " << " Position [km] " << endl;

  for (int i = 0; i <= N_step; i++) {
    mjd_tt = mjd0 + i * dt;
    Real tt = (mjd_tt - MJD_J2000) * SECS_DAY;
    Real t_tai = ConvertTime(tt, Time::TT, Time::TAI);
    r = GetBodyPosVel(t_tai, NaifId::EARTH, NaifId::MOON, Frame::GCRF).head(3);
    cout << " " << MJD2GregorianString(mjd_tt, 1) << "      " << r.transpose().format(fmt) << endl;
  };

  return 0;
}

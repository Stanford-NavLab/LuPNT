#include <lupnt/lupnt.h>

using namespace lupnt;
using namespace std;

int main() {
  Real t_tt = 0.0;  // [s] J2000.0 TT time
  Real t_tai = TT2TAI(t_tt);
  Real t_utc = TAI2UTC(t_tai);
  Real t_tdb = TT2TDB(t_tt);
  Real t_tai_spice = ConvertTime(t_tt, Time::TT, Time::TAI);
  Real t_tdb_spice = ConvertTime(t_tt, Time::TT, Time::TDB);
  Real mjd_tt = t_tt / SECS_DAY + MJD_J2000;
  Real mjd_tai = t_tai / SECS_DAY + MJD_J2000;
  Real mjd_utc = t_utc / SECS_DAY + MJD_J2000;
  Real mjd_tdb = t_tdb / SECS_DAY + MJD_J2000;
  int prec1 = 12;
  int prec2 = 13;
  int prec3 = 9;
  int w1 = 18;
  int w2 = 21;
  cout << fixed;
  cout << "TT  = ";
  cout << setw(w1) << setprecision(prec1) << t_tt << " s  ";
  cout << setw(w2) << setprecision(prec2) << mjd_tt << " days  ";
  cout << MJD2GregorianString(mjd_tt, prec3) << endl;
  cout << "TAI = ";
  cout << setw(w1) << setprecision(prec1) << t_tai << " s  ";
  cout << setw(w2) << setprecision(prec2) << mjd_tai << " days  ";
  cout << MJD2GregorianString(mjd_tai, prec3) << endl;
  cout << "UTC = ";
  cout << setw(w1) << setprecision(prec1) << t_utc << " s  ";
  cout << setw(w2) << setprecision(prec2) << mjd_utc << " days  ";
  cout << MJD2GregorianString(mjd_utc, prec3) << endl;
  cout << "TDB = ";
  cout << setw(w1) << setprecision(prec1) << t_tdb << " s  ";
  cout << setw(w2) << setprecision(prec2) << mjd_tdb << " days  ";
  cout << MJD2GregorianString(mjd_tdb, prec3) << endl;

  cout << "SPICE:" << endl;
  cout << "TAI = ";
  cout << setw(w1) << setprecision(prec1) << t_tai_spice << " s  " << endl;
  cout << "TDB = ";
  cout << setw(w1) << setprecision(prec1) << t_tdb_spice << " s  " << endl;

  return 0;
}

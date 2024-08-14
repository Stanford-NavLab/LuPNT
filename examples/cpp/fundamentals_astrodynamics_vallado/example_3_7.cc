#include <lupnt/lupnt.h>

using namespace lupnt;
using namespace std;

int main() {
  Real mjd_utc = GregorianToMJD(2004, 5, 14, 16, 43, 0);
  Real t_utc = MJDtoTime(mjd_utc);
  Real t_ut1 = UTCtoUT1(t_utc);
  Real t_tai = UTCtoTAI(t_utc);
  Real t_gps = TAItoGPS(t_tai);
  Real t_tt = TAItoTT(t_tai);
  Real t_tdb = TTtoTDB(t_tt);

  int p1 = 4;
  int w1 = 12;
  int p2 = 12;
  int w2 = 15;
  cout << fixed << setprecision(p1);
  cout << "UT1 = " << setw(w1) << TimeToGregorianString(t_ut1, p1);
  cout << fixed << setprecision(p2);
  cout << "  T_UT1 = " << t_ut1 / SECS_DAY / DAYS_CENTURY << endl;
  cout << fixed << setprecision(p1);
  cout << "UTC = " << setw(w1) << TimeToGregorianString(t_utc, p1) << endl;
  cout << "TAI = " << setw(w1) << TimeToGregorianString(t_tai, p1) << endl;
  cout << "GPS = " << setw(w1) << TimeToGregorianString(t_gps, p1) << endl;
  cout << fixed << setprecision(p1);
  cout << "TT  = " << setw(w1) << TimeToGregorianString(t_tt, p1);
  cout << fixed << setprecision(p2);
  cout << "  T_TT  = " << t_tt / SECS_DAY / DAYS_CENTURY << endl;
  cout << fixed << setprecision(p1);
  cout << "TDB = " << setw(w1) << TimeToGregorianString(t_tdb, p1);
  cout << fixed << setprecision(p2);
  cout << "  T_TDB = " << t_tdb / SECS_DAY / DAYS_CENTURY << endl;
  return 0;
}

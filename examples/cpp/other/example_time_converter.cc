#include <lupnt/lupnt.h>

#include <iostream>
#include <string>

using namespace lupnt;
using namespace std;

int main() {
  Real t_tai, t_tdb, t_tt, t_gps;              // [s]
  Real t_tai_sp, t_tdb_sp, t_tt_sp, t_gps_sp;  // [s]
  Real jd_tdb, jd_tt;                          // [days]
  Real jd_tdb_sp, jd_tt_sp;                    // [days]

  // J2000
  t_tt_sp = 1.23456789e12;
  t_tdb_sp = spice::ConvertTime(t_tt_sp, Time::TT, Time::TDB);
  t_tai_sp = spice::ConvertTime(t_tt_sp, Time::TT, Time::TAI);
  t_gps_sp = spice::ConvertTime(t_tt_sp, Time::TT, Time::GPS);
  jd_tt_sp = spice::ConvertTime(t_tt_sp, Time::TT, Time::JD_TT);
  jd_tdb_sp = spice::ConvertTime(t_tt_sp, Time::TT, Time::JD_TDB);

  t_tt = t_tt_sp;
  t_tdb = ConvertTime(t_tt, Time::TT, Time::TDB);
  t_tai = ConvertTime(t_tt, Time::TT, Time::TAI);
  t_gps = ConvertTime(t_tt, Time::TT, Time::GPS);
  jd_tt = Time2JD(t_tt);
  jd_tdb = Time2JD(t_tdb);

  int w = 21;
  cout << fixed << setprecision(14);
  cout << endl << "Spice/LuPNT time conversion:" << endl << endl;
  cout << "t_tt   = " << t_tt_sp << " s" << endl;
  cout << "       = " << t_tt << " s" << endl;
  cout << "    Δt = " << t_tt_sp - t_tt << " s" << endl << endl;

  cout << "t_tdb  = " << t_tdb_sp << " s" << endl;
  cout << "       = " << t_tdb << " s" << endl;
  cout << "    Δt = " << t_tdb_sp - t_tdb << " s" << endl << endl;

  cout << "t_tai  = " << t_tai_sp << " s" << endl;
  cout << "       = " << t_tai << " s" << endl;
  cout << "    Δt = " << t_tai_sp - t_tai << " s" << endl << endl;

  cout << "t_gps  = " << t_gps_sp << " s" << endl;
  cout << "       = " << t_gps << " s" << endl;
  cout << "    Δt = " << t_gps_sp - t_gps << " s" << endl << endl;

  cout << "jd_tt  = " << jd_tt_sp << " days" << endl;
  cout << "       = " << jd_tt << " days" << endl;
  cout << "   Δjd = " << (jd_tt_sp - jd_tt) << " days" << endl << endl;

  cout << "jd_tdb = " << jd_tdb << " days" << endl;
  cout << "       = " << jd_tdb_sp << " days" << endl;
  cout << "   Δjd = " << (jd_tdb_sp - jd_tdb) << " days" << endl << endl;

  return 0;
}

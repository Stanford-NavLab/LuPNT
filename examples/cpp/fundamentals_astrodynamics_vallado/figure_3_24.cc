#include <lupnt/lupnt.h>

using namespace lupnt;
using namespace std;
using namespace matplot;

int main() {
  int year0 = 1961;
  int yearf = 2017;
  int N_points = 100;
  int dyears = 4;

  Real t0_tai = GregorianToTime(year0, 1, 1, 0, 0, 0);
  Real tf_tai = GregorianToTime(yearf, 1, 1, 0, 0, 0);
  VecX t_tai = VecX::LinSpaced(N_points, t0_tai, tf_tai);
  VecX t_utc = TAItoUTC(t_tai);
  VecX t_ut1 = UTCtoUT1(t_utc);
  VecX t_gps = TAItoGPS(t_tai);
  VecX t_tt = TAItoTT(t_tai);
  VecX t_tdb = TTtoTDB(t_tt);
  VecX t_tcg = TTtoTCG(t_tt);
  VecX t_tcb = TTtoTCB(t_tt);

  VecX years_plot = VecX::LinSpaced(N_points, year0, yearf);

  vector<double> years_xticks;
  for (int y = year0; y <= yearf; y += dyears) years_xticks.push_back(y);
  vector<std::string> years_xticklabels;  // Jan-XX
  for (int y = year0; y <= yearf; y += dyears)
    years_xticklabels.push_back(string("Jan-") + to_string(y).substr(2, 2));

  // Figure with size
  auto f = figure(true);
  f->size(800, 600);

  hold(on);
  line_handle p;
  double lw = 2;
  lupnt::plot(years_plot, t_tai - t_tai)->line_width(lw).display_name("TAI");
  lupnt::plot(years_plot, t_utc - t_tai)->line_width(lw).display_name("UTC");
  lupnt::plot(years_plot, t_ut1 - t_tai)->line_width(lw).display_name("UT1");
  lupnt::plot(years_plot, t_gps - t_tai)->line_width(lw).display_name("GPS");
  lupnt::plot(years_plot, t_tt - t_tai)->line_width(lw).display_name("TT");
  lupnt::plot(years_plot, t_tdb - t_tai)->line_width(lw).display_name("TDB");
  lupnt::plot(years_plot, t_tcg - t_tai)->line_width(lw).display_name("TCG");
  lupnt::plot(years_plot, t_tcb - t_tai)->line_width(lw).display_name("TCB");

  ylabel("Difference in Time to TAI [s]");
  xticks(years_xticks);
  xticklabels(years_xticklabels);
  grid(on);
  ::matplot::legend("northwest")->location(legend::general_alignment::topleft);

  show();
}

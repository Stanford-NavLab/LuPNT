#include <lupnt/lupnt.h>

using namespace lupnt;
using namespace std;
using namespace matplot;

int main() {
  int year0 = 1961;
  int yearf = 2017;
  int N_points = 100;
  int dyears = 4;

  Real t0_tai = Gregorian2Time(year0, 1, 1, 0, 0, 0);
  Real tf_tai = Gregorian2Time(yearf, 1, 1, 0, 0, 0);
  VecX t_tai = VecX::LinSpaced(N_points, t0_tai, tf_tai);
  VecX t_utc = TAI2UTC(t_tai);
  VecX t_ut1 = UTC2UT1(t_utc);
  VecX t_gps = TAI2GPS(t_tai);
  VecX t_tt = TAI2TT(t_tai);
  VecX t_tdb = TT2TDB(t_tt);
  VecX t_tcg = TT2TCG(t_tt);
  VecX t_tcb = TT2TCB(t_tt);

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
  lupnt::Plot(years_plot, t_tai - t_tai)->line_width(lw).display_name("TAI");
  lupnt::Plot(years_plot, t_utc - t_tai)->line_width(lw).display_name("UTC");
  lupnt::Plot(years_plot, t_ut1 - t_tai)->line_width(lw).display_name("UT1");
  lupnt::Plot(years_plot, t_gps - t_tai)->line_width(lw).display_name("GPS");
  lupnt::Plot(years_plot, t_tt - t_tai)->line_width(lw).display_name("TT");
  lupnt::Plot(years_plot, t_tdb - t_tai)->line_width(lw).display_name("TDB");
  lupnt::Plot(years_plot, t_tcg - t_tai)->line_width(lw).display_name("TCG");
  lupnt::Plot(years_plot, t_tcb - t_tai)->line_width(lw).display_name("TCB");

  ylabel("Difference in Time to TAI [s]");
  xticks(years_xticks);
  xticklabels(years_xticklabels);
  grid(on);
  ::matplot::legend("northwest")->location(legend::general_alignment::topleft);

  show();
}

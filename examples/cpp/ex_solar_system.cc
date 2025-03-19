#include <lupnt/lupnt.h>
#include <matplot/matplot.h>

using namespace lupnt;
using namespace std;
using namespace matplot;

int main() {
  // Solar System
  Real t0 = Gregorian2Time(2000, 1, 1, 0, 0, 0);
  Real dt_total = 10 * DAYS_YEAR * SECS_DAY;
  Real dt_step = 1 * SECS_DAY;
  VecX tspan = arange(0, dt_total + dt_step, dt_step);
  VecX tfs = t0 + tspan.array();

  vector<NaifId> bodies
      = {NaifId::MERCURY_BARYCENTER, NaifId::VENUS_BARYCENTER,   NaifId::EARTH,
         NaifId::MARS_BARYCENTER,    NaifId::JUPITER_BARYCENTER, NaifId::SATURN_BARYCENTER,
         NaifId::URANUS_BARYCENTER,  NaifId::NEPTUNE_BARYCENTER};

  figure_handle fig;

  fig = figure(true);
  title("Solar System");
  hold(true);
  grid(true);
  for (auto body : bodies) {
    MatX6 rv = GetBodyPosVel(tfs, body, Frame::GCRF);
    Plot3(rv.col(0), rv.col(1), rv.col(2), "-", 9);
  }
  matplot::legend({"Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"});
  SetLim(1.5e9, 9);
  xlabel("X [10^9 km]");
  ylabel("Y [10^9 km]");
  zlabel("Z [10^9 km]");
  fig->draw();

  // Moon Angles
  t0 = Gregorian2Time(2000, 1, 1, 0, 0, 0);
  dt_total = 1 * DAYS_YEAR * SECS_DAY;
  dt_step = 1 * SECS_DAY;
  tspan = arange(0, dt_total + dt_step, dt_step);
  tfs = t0 + tspan.array();

  MatX6 rv = GetLunarMantleData(tfs);

  fig = figure(true);
  title("Lunar Mantle Angles");
  hold(true);
  grid(true);
  Plot(tspan / SECS_DAY, Wrap2Pi(rv.col(0)) * DEG, "-");
  Plot(tspan / SECS_DAY, Wrap2Pi(rv.col(1)) * DEG, "-");
  Plot(tspan / SECS_DAY, Wrap2Pi(rv.col(2)), "-");
  // legend({"phi", "theta", "psi"});
  matplot::legend({"phi", "theta", "psi"});
  xlabel("Time [days]");
  ylabel("Angle [deg]");
  fig->draw();

  fig = figure(true);
  title("Lunar Mantle Angular Velocity");
  hold(true);
  grid(true);
  Plot(tspan / SECS_DAY, rv.col(3) * DEG, "-");
  Plot(tspan / SECS_DAY, rv.col(4) * DEG, "-");
  Plot(tspan / SECS_DAY, rv.col(5) * DEG, "-");
  // Plot(tspan / SECS_DAY, rv.col(2), "-");
  // legend({"phi", "theta", "psi"});
  matplot::legend({"phi dot", "theta dot", "psi dot"});
  xlabel("Time [days]");
  ylabel("Angle [deg]]");
  fig->draw();

  show();
  return 0;
}

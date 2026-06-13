#include <lupnt/lupnt.h>
#include <matplot/matplot.h>

using namespace lupnt;
using namespace std;
using namespace matplot;

int main() {
  // Plot barycentric/inertial planet trajectories over ten years. The example
  // is primarily a visual smoke test for SPICE body queries and matplot output.
  Real t0 = GregorianToTime(2000, 1, 1, 0, 0, 0);
  Real dt_total = 10 * DAYS_YEAR * SECS_DAY;
  Real dt_step = 1 * SECS_DAY;
  VecX tspan = Arange(Real(0.0), dt_total + dt_step, dt_step);
  VecX tfs = t0 + tspan.array();

  vector<BodyId> bodies
      = {BodyId::MERCURY_BARYCENTER, BodyId::VENUS_BARYCENTER,   BodyId::EARTH,
         BodyId::MARS_BARYCENTER,    BodyId::JUPITER_BARYCENTER, BodyId::SATURN_BARYCENTER,
         BodyId::URANUS_BARYCENTER,  BodyId::NEPTUNE_BARYCENTER};

  figure_handle fig;

  fig = figure(true);
  title("Solar System");
  hold(true);
  grid(true);
  for (auto body : bodies) {
    MatX6 rv = GetBodyPosVel(tfs, body, Frame::ICRF);
    rv /= AU;
    Plot3(rv.col(0).eval(), rv.col(1).eval(), rv.col(2).eval(), "-");
  }
  matplot::legend({"Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"});
  // SetLim(1.5);
  xlabel("X [AU]");
  ylabel("Y [AU]");
  zlabel("Z [AU]");
  fig->draw();

  // Query the lunar mantle Euler angles over one year.
  t0 = GregorianToTime(2000, 1, 1, 0, 0, 0);
  dt_total = 1 * DAYS_YEAR * SECS_DAY;
  dt_step = 1 * SECS_DAY;
  tspan = Arange(Real(0.0), dt_total + dt_step, dt_step);
  tfs = t0 + tspan.array();

  MatX6 ang = GetLunarMantleData(tfs);

  fig = figure(true);
  title("Lunar Mantle Angles");
  hold(true);
  grid(true);
  Plot(tspan / SECS_DAY, WrapToPi(ang.col(0)) * DEG, "-");
  Plot(tspan / SECS_DAY, WrapToPi(ang.col(1)) * DEG, "-");
  Plot(tspan / SECS_DAY, WrapToPi(ang.col(2)), "-");
  // legend({"phi", "theta", "psi"});
  matplot::legend({"phi", "theta", "psi"});
  xlabel("Time [days]");
  ylabel("Angle [deg]");
  fig->draw();

  fig = figure(true);
  title("Lunar Mantle Angular Velocity");
  hold(true);
  grid(true);
  Plot(tspan / SECS_DAY, ang.col(3) * DEG, "-");
  Plot(tspan / SECS_DAY, ang.col(4) * DEG, "-");
  Plot(tspan / SECS_DAY, ang.col(5) * DEG, "-");
  // Plot(tspan / SECS_DAY, rv.col(2), "-");
  // legend({"phi", "theta", "psi"});
  matplot::legend({"phi dot", "theta dot", "psi dot"});
  xlabel("Time [days]");
  ylabel("Angle [deg]");
  fig->draw();

  show();
  return 0;
}

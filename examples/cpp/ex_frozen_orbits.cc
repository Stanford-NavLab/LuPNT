#include <lupnt/lupnt.h>
#include <matplot/matplot.h>

using namespace lupnt;
using namespace std;
namespace plt = matplot;

class Case1Dynamics : public CartesianTwoBodyDynamics {
public:
  Case1Dynamics(Real t0, Vec6 coe0_moon_ci, IntegratorType integ = IntegratorType::RK4)
      : CartesianTwoBodyDynamics(GM_MOON, integ),
        t0_(t0),
        coe0_moon_ci_(coe0_moon_ci),
        dyn_moon_(GM_EARTH) {
    SetODEFunction([this](Real t, const Vec6 &rv) { return ComputeRates(t, rv); });
  }

  Vec6 ComputeRates(Real t, const Vec6 &rv) {
    Vec6 coe_moon = dyn_moon_.Propagate(coe0_moon_ci_, t0_, t);
    Vec6 rv_earth = -Classical2Cart(coe_moon, GM_EARTH);

    Vec6 rv_dot = CartesianTwoBodyDynamics::ComputeRates(t, rv);
    rv_dot.tail(3) += AccelerationPointMass(rv.head(3), rv_earth.head(3), GM_EARTH);
    return rv_dot;
  }

private:
  Real t0_;
  Vec6 coe0_moon_ci_;
  KeplerianDynamics dyn_moon_;
};

int main() {
  Real t0 = Gregorian2Time(2009, 7, 1, 1, 0, 0);  // [s] Start time (TAI)

  Real a = 6541.4;      // [km] Semi-major axis
  Real e = 0.60;        // [-] Eccentricity
  Real i = 56.2 * RAD;  // [deg] Inclination
  Real O = 0.00 * RAD;  // [deg] Right ascension of the ascending node
  Real w = 90.0 * RAD;  // [deg] Argument of perigee
  Real M = 0.00 * RAD;  // [deg] Mean anomaly

  Vec6 coe0_op(a, e, i, O, w, M);
  Vec6 rv0_op = Classical2Cart(coe0_op, GM_MOON);
  Vec6 rv0_ci = ConvertFrame(t0, rv0_op, MOON_OP, MOON_CI);
  Vec6 coe0_ci = Cart2Classical(rv0_ci, GM_MOON);

  // **************************************************************************
  // Case 0
  // **************************************************************************

  Real dt_total = 2 * SECS_DAY;     // [s] Total propagation time
  Real dt_step = 20 * SECS_MINUTE;  // [s] Time step
  Real dt_prop = 60;                // [s] Propagation time step

  VecX tspan = arange(0, dt_total, dt_step);  // [s] Time span
  VecX tfs = t0 + tspan.array();              // [s] Final times

  CartesianTwoBodyDynamics dyn0(GM_MOON);
  dyn0.SetTimeStep(dt_prop);
  MatX6 rv_case0_ci = dyn0.Propagate(rv0_ci, t0, tfs);

  Vec6 rv0_earth2moon_ci = GetBodyPosVel(t0, EARTH, MOON, MOON_CI);
  Vec3 e_earth2moon_ci = rv0_earth2moon_ci.head(3).normalized();

  // Plot
  auto fig0 = plt::figure(true);
  plt::hold(true);
  plot_body(MOON);
  plot3(rv_case0_ci.col(0), rv_case0_ci.col(1), rv_case0_ci.col(2));
  plot3(e_earth2moon_ci * 2 * R_MOON, "r");
  set_lim(12e3);
  fig0->draw();

  // **************************************************************************
  // Case 1
  // **************************************************************************

  Real moon_period = 2 * PI * sqrt(pow(D_EARTH_MOON, 3) / GM_EARTH);
  cout << "Moon period: " << moon_period / SECS_DAY << " days" << endl;

  // Time
  dt_total = 1.25 * DAYS_YEAR * SECS_DAY;  // [s] Total propagation time
  dt_step = 60 * SECS_MINUTE;              // [s] Time step
  dt_prop = 180;                           // [s] Propagation time step

  tspan = arange(0, dt_total, dt_step);  // [s] Time span
  tfs = t0 + tspan.array();              // [s] Final times

  // Moon
  KeplerianDynamics dyn_moon(GM_EARTH);
  Vec6 coe0_moon_ci = Cart2Classical(rv0_earth2moon_ci, GM_EARTH);
  MatX6 coe_moon_ci = dyn_moon.Propagate(coe0_moon_ci, t0, tfs);
  MatX6 rv_moon_ci = Classical2Cart(coe_moon_ci, GM_EARTH);

  // Plot
  auto fig11 = plt::figure(true);
  plt::hold(true);
  plot_body(EARTH);
  plot_body(MOON, rv0_earth2moon_ci.head(3));
  plot3(rv_moon_ci.col(0), rv_moon_ci.col(1), rv_moon_ci.col(2));
  set_lim(D_EARTH_MOON * 1.1);
  fig11->draw();

  Case1Dynamics dyn_case1(t0, coe0_moon_ci, IntegratorType::RK4);
  dyn_case1.SetTimeStep(dt_prop);
  cout << "Propagating Case 1..." << endl;
  MatX6 rv_case1_ci = dyn_case1.Propagate(rv0_ci, t0, tfs, true);

  // Plot
  auto fig12 = plt::figure(true);
  plt::hold(true);
  plot_body(MOON);
  plot3(rv_case1_ci.col(0), rv_case1_ci.col(1), rv_case1_ci.col(2));
  set_lim(12e3);
  fig12->draw();

  MatX6 rv_case1_op = ConvertFrame(t0, rv_case1_ci, MOON_CI, MOON_OP);
  MatX6 coe_case1_op = Cart2Classical(rv_case1_op, GM_MOON);
  VecX e_vec = coe_case1_op.col(1).array();
  VecX w_vec = coe_case1_op.col(4).array();

  // Plot
  auto fig13 = plt::figure(true);
  plot(w_vec, e_vec);
  plt::xlabel("Argument of perigee, w [rad]");
  plt::ylabel("Eccentricity, e [-]");

  plt::show();
  return 0;
}

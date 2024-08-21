#include <lupnt/lupnt.h>
#include <matplot/matplot.h>

#include <highfive/highfive.hpp>

using namespace lupnt;
using namespace std;
using namespace HighFive;
namespace plt = matplot;

class Case1Dynamics : public CartesianTwoBodyDynamics {
public:
  Case1Dynamics(Real t0, Vec6 coe0_moon_op, IntegratorType integ = IntegratorType::RK4)
      : CartesianTwoBodyDynamics(GM_MOON, integ),
        t0_(t0),
        coe0_moon_op_(coe0_moon_op),
        dyn_moon_(GM_EARTH) {
    SetODEFunction([this](Real t, const Vec6 &rv) { return ComputeRates(t, rv); });
  }

  Vec6 ComputeRates(Real t, const Vec6 &rv) {
    Vec6 coe_moon_op = dyn_moon_.Propagate(coe0_moon_op_, t0_, t);
    Vec6 rv_earth_op = -Classical2Cart(coe_moon_op, GM_EARTH);

    Vec6 rv_dot = CartesianTwoBodyDynamics::ComputeRates(t, rv);
    rv_dot.tail(3) += AccelerationPointMass(rv.head(3), rv_earth_op.head(3), GM_EARTH);
    return rv_dot;
  }

private:
  Real t0_;
  Vec6 coe0_moon_op_;
  KeplerianDynamics dyn_moon_;
};

int main() {
  // Options
  File file("ex_frozen_orbits.h5", File::OpenOrCreate);
  // Get the root group
  Group root = file.getGroup("/");
  // Iterate over the groups in the root group
  for (const auto &group : root.listObjectNames()) {
    std::cout << "Group: " << group << std::endl;
  }

  bool recompute_case1 = false;
  bool recompute_case2 = false;

  // Time
  Real t0 = Gregorian2Time(2009, 7, 1, 1, 0, 0);  // [s] Start time (TAI)

  // Orbital elements
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

  Real sat_period = 2 * PI * sqrt(pow(a, 3) / GM_MOON);
  cout << "Satellite period: " << sat_period / SECS_MINUTE << " minutes" << endl;

  // **************************************************************************
  // Case 0
  // **************************************************************************
  cout << endl << endl << "*********** Case 0 ***********" << endl << endl;

  Real dt_total = sat_period;                           // [s] Total propagation time
  Real dt_step = 5 * SECS_MINUTE;                       // [s] Time step
  Real dt_prop = 20;                                    // [s] Propagation time step
  VecX tspan = arange(0, dt_total + dt_step, dt_step);  // [s] Time span
  VecX tfs = t0 + tspan.array();                        // [s] Final times
  int n_steps = tspan.size();
  cout << "Total duration   " << dt_total / SECS_HOUR << " hours" << endl;
  cout << "Time step        " << dt_step / SECS_MINUTE << " minutes" << endl;
  cout << "Propagation step " << dt_prop << " seconds" << endl;
  cout << "Start epoch      " << Time2GregorianString(t0) << endl;
  cout << "End epoch        " << Time2GregorianString(t0 + dt_total) << endl;
  cout << "Number of steps  " << n_steps << endl;

  CartesianTwoBodyDynamics dyn0(GM_MOON);
  dyn0.SetTimeStep(dt_prop);
  MatX6 rv_case0_ci = dyn0.Propagate(rv0_ci, t0, tfs);

  Vec6 rv0_moon2earth_ci = GetBodyPosVel(t0, MOON, EARTH, MOON_CI);
  Vec3 e_moon2earth_ci = rv0_moon2earth_ci.head(3).normalized();

  // Plot
  plt::line_handle p;
  auto fig0 = plt::figure(true);
  plt::hold(true);
  PlotBody(MOON);
  p = Plot3(rv_case0_ci.col(0), rv_case0_ci.col(1), rv_case0_ci.col(2), "b-.");
  p->line_width(2);
  p->marker_indices({0});
  p = PlotArrow3(e_moon2earth_ci * 2 * R_MOON, "g-");
  p->line_width(2);
  SetLim(12e3);
  fig0->draw();

  // **************************************************************************
  // Case 1
  // **************************************************************************
  cout << endl << endl << "*********** Case 1 ***********" << endl << endl;

  Real moon_period = 2 * PI * sqrt(pow(D_EARTH_MOON, 3) / GM_EARTH);
  cout << "Moon period: " << moon_period / SECS_DAY << " days" << endl << endl;

  // Time
  dt_total = 2 * DAYS_YEAR * SECS_DAY;   // [s] Total propagation time
  dt_step = 30 * SECS_MINUTE;            // [s] Time step
  dt_prop = 60;                          // [s] Propagation time step
  tspan = arange(0, dt_total, dt_step);  // [s] Time span
  tfs = t0 + tspan.array();              // [s] Final times
  n_steps = tspan.size();
  cout << "Total duration   " << dt_total / SECS_DAY << " days" << endl;
  cout << "Time step        " << dt_step / SECS_MINUTE << " minutes" << endl;
  cout << "Propagation step " << dt_prop << " seconds" << endl;
  cout << "Start epoch      " << Time2GregorianString(t0) << endl;
  cout << "End epoch        " << Time2GregorianString(t0 + dt_total) << endl;
  cout << "Number of steps  " << n_steps << endl;

  // Moon
  Vec6 rv0_moon_op = GetBodyPosVel(t0, EARTH, MOON, MOON_OP);
  Vec6 coe0_moon_op = Cart2Classical(rv0_moon_op, GM_EARTH);

  KeplerianDynamics dyn_moon(GM_EARTH);
  MatX6 coe_moon_op = dyn_moon.Propagate(coe0_moon_op, t0, tfs);
  MatX6 rv_moon_op = Classical2Cart(coe_moon_op, GM_EARTH);

  // Plot
  auto fig11 = plt::figure(true);
  plt::hold(true);
  PlotBody(EARTH);
  PlotBody(MOON, rv0_moon_op.head(3));
  p = Plot3(rv_moon_op.col(0), rv_moon_op.col(1), rv_moon_op.col(2), "b-o");
  p->marker_indices({0});
  p->line_width(2);
  SetLim(D_EARTH_MOON * 1.1);
  fig11->draw();

  // Propagate
  Case1Dynamics dyn_case1(t0, coe0_moon_op, IntegratorType::RK4);
  dyn_case1.SetTimeStep(dt_prop);
  MatX6 rv_case1_op;
  if (recompute_case1 || !file.exist("/case1/rv_case1_op")) {
    cout << endl << "Propagating" << endl;
    rv_case1_op = dyn_case1.Propagate(rv0_op, t0, tfs, true);
    dump(file, "/case1/rv_case1_op", rv_case1_op.cast<double>());
  } else {
    rv_case1_op = load<MatX6d>(file, "/case1/rv_case1_op");
    cout << "Loaded from file" << endl;
  }
  MatX6 coe_case1_op = Cart2Classical(rv_case1_op, GM_MOON);

  // Plo

  // Plot orbital elements
  auto fig12 = plt::figure(true);
  plt::hold(true);
  vector<string> labels = {"a [km]", "e [-]", "i [deg]", "O [deg]", "w [deg]", "M [deg]"};
  for (int i = 0; i < 6; i++) {
    plt::subplot(3, 2, i);
    VecX x = tspan / SECS_DAY;
    VecX y;
    if (i < 2)
      y = coe_case1_op.col(i);
    else
      y = coe_case1_op.col(i) * DEG;

    Plot(x, y);
    plt::xlabel("Days past " + Time2GregorianString(t0) + " TAI");
    plt::ylabel(labels[i]);
    if (i == 5)
      plt::xlim({0, 10});
    else
      plt::xlim({0, dt_total.val() / SECS_DAY});
    plt::grid(true);
  }
  fig12->draw();

  // Plot orbits
  int n_plot = 20;
  tspan = arange(0, sat_period + dt_prop, dt_prop);
  vector<MatX6> rv_plot_op;
  for (int i = 0; i < n_plot; i++) {
    int j = i * n_steps / n_plot;
    VecX tfs_plot = tfs[j] + tspan.array();
    Vec6 rv0_plot = rv_case1_op.row(j);
    rv_plot_op.push_back(dyn_case1.Propagate(rv0_plot, tfs[j], tfs_plot));
  }

  auto fig13 = plt::figure(true);
  plt::hold(true);
  PlotBody(MOON);
  for (int i = 0; i < n_plot; i++) {
    p = Plot3(rv_plot_op[i].col(0), rv_plot_op[i].col(1), rv_plot_op[i].col(2), "b-");
    p->line_width(2);
  }
  SetLim(12e3);
  fig13->draw();

  // **************************************************************************
  // Case 2
  // **************************************************************************
  cout << endl << endl << "*********** Case 2 ***********" << endl << endl;

  NBodyDynamics dyn_case2(IntegratorType::RK4);
  dyn_case2.AddBody(Body::Moon());
  dyn_case2.AddBody(Body::Earth());
  dyn_case2.SetTimeStep(dt_prop);
  dyn_case2.SetFrame(MOON_CI);

  // Propagate
  MatX6 rv_case2_ci;
  if (recompute_case2 || !file.exist("/case2/rv_case2_ci")) {
    cout << endl << "Propagating" << endl;
    rv_case2_ci = dyn_case2.Propagate(rv0_ci, t0, tfs, true);
    dump(file, "/case2/rv_case2_ci", rv_case2_ci.cast<double>());
  } else {
    rv_case2_ci = load<MatX6d>(file, "/case2/rv_case2_ci");
    cout << "Loaded from file" << endl;
  }
  MatX6 rv_case2_op = ConvertFrame(tfs, rv_case2_ci, Frame::MOON_CI, Frame::MOON_OP);
  MatX6 coe_case2_op = Cart2Classical(rv_case2_op, GM_MOON);

  // Plot orbital elements
  auto fig21 = plt::figure(true);
  plt::hold(true);
  for (int i = 0; i < 6; i++) {
    plt::subplot(3, 2, i);
    VecX x = tspan / SECS_DAY;
    VecX y;
    if (i < 2)
      y = coe_case2_op.col(i);
    else
      y = coe_case2_op.col(i) * DEG;

    Plot(x, y);
    plt::xlabel("Days past " + Time2GregorianString(t0) + " TAI");
    plt::ylabel(labels[i]);
    if (i == 5)
      plt::xlim({0, 10});
    else
      plt::xlim({0, dt_total.val() / SECS_DAY});
    plt::grid(true);
  }
  fig21->draw();

  plt::show();
  return 0;
  // Plot
  // auto fig12 = plt::figure(true);
  // plt::hold(true);
  // PlotBody(MOON);
  // Plot3(rv_case1_ci.col(0), rv_case1_ci.col(1), rv_case1_ci.col(2));
  // SetLim(12e3);
  // fig12->draw();

  VecX e_vec = coe_case1_op.col(1).array();
  VecX w_vec = coe_case1_op.col(4).array();

  // Plot
  auto fig14 = plt::figure(true);
  Plot(w_vec, e_vec);
  plt::xlabel("Argument of perigee, w [rad]");
  plt::ylabel("Eccentricity, e [-]");

  plt::show();
  return 0;
}

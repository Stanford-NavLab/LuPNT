#include <lupnt/lupnt.h>
#include <matplot/matplot.h>
#include <omp.h>

#include <highfive/highfive.hpp>

using namespace lupnt;
using namespace std;
using namespace HighFive;
using namespace matplot;

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

struct {
  bool recompute_part1 = true;
  bool recompute_part2 = true;
  bool plot_case0 = false;
  bool plot_case1 = false;
  bool plot_case2 = false;
  bool plot_case3 = false;
  bool plot_case4 = false;
} config;

int main() {
  auto begin = GetSystemTime();

  auto output_path = GetOutputPath("ex_frozen_orbits");
  cout << "Output path: " << output_path << endl;

  auto open_mode_part1 = (config.recompute_part1) ? File::Truncate : File::OpenOrCreate;
  File file_part1(output_path / "data_part1.h5", open_mode_part1);

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

  Real sat_period = GetOrbitalPeriod(a, GM_MOON);
  cout << "Satellite period: " << sat_period / SECS_MINUTE << " minutes" << endl;

  // **************************************************************************
  // Case 0
  // **************************************************************************
  cout << endl << endl << "*********** Case 0 ***********" << endl;

  // Time
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

  // Initial state
  Vec6 rv0_op = Classical2Cart(coe0_op, GM_MOON);
  Vec6 rv0_ci = ConvertFrame(t0, rv0_op, MOON_OP, MOON_CI);

  // Dynamics
  CartesianTwoBodyDynamics dyn0(GM_MOON);
  dyn0.SetTimeStep(dt_prop);

  // Propagate
  MatX6 rv_case0_ci = dyn0.Propagate(rv0_ci, t0, tfs);

  // Plot
  line_handle p;
  figure_handle fig;
  if (config.plot_case0) {
    fig = figure(true);
    title("Case 0: Moon orbit in CI frame");
    hold(true);
    PlotBody(MOON);
    p = Plot3(rv_case0_ci.col(0), rv_case0_ci.col(1), rv_case0_ci.col(2), "b-");
    p->line_width(2);
    p->marker_indices({0});
    SetLim(12e3);
    fig->draw();
  }

  // **************************************************************************
  // Case 1
  // **************************************************************************
  cout << endl << endl << "*********** Case 1 ***********" << endl;

  Real moon_period = GetOrbitalPeriod(D_EARTH_MOON, GM_EARTH);  // [s] Moon period
  cout << "Moon period: " << moon_period / SECS_DAY << " days" << endl;

  // Time
  dt_total = 2 * DAYS_YEAR * SECS_DAY;             // [s] Total propagation time
  dt_step = 60 * SECS_MINUTE;                      // [s] Time step
  dt_prop = 2 * SECS_MINUTE;                       // [s] Propagation time step
  tspan = arange(0, dt_total + dt_step, dt_step);  // [s] Time span
  tfs = t0 + tspan.array();                        // [s] Final times
  n_steps = tspan.size();
  cout << "Total duration   " << dt_total / SECS_DAY << " days" << endl;
  cout << "Time step        " << dt_step / SECS_MINUTE << " minutes" << endl;
  cout << "Propagation step " << dt_prop << " seconds" << endl;
  cout << "Start epoch      " << Time2GregorianString(t0) << endl;
  cout << "End epoch        " << Time2GregorianString(t0 + dt_total) << endl;
  cout << "Number of steps  " << n_steps << endl;

  // Initial state
  Vec6 rv0_moon_op = GetBodyPosVel(t0, EARTH, MOON, MOON_OP);
  Vec6 coe0_moon_op = Cart2Classical(rv0_moon_op, GM_EARTH);

  // Dynamics
  KeplerianDynamics dyn_moon(GM_EARTH);

  // Propagate
  MatX6 coe_moon_op = dyn_moon.Propagate(coe0_moon_op, t0, tfs);
  MatX6 rv_moon_op = Classical2Cart(coe_moon_op, GM_EARTH);

  if (config.plot_case1) {
    // Plot
    fig = figure(true);
    title("Case 1: Moon orbit in OP frame");
    hold(true);
    PlotBody(EARTH);
    PlotBody(MOON, rv0_moon_op.head(3));
    p = Plot3(rv_moon_op.col(0), rv_moon_op.col(1), rv_moon_op.col(2), "b-o");
    p->marker_indices({0});
    p->line_width(2);
    SetLim(D_EARTH_MOON * 1.1);
    fig->draw();
  }

  // Dynamics
  Case1Dynamics dyn_3body_circ(t0, coe0_moon_op, IntegratorType::RK4);
  dyn_3body_circ.SetTimeStep(dt_prop);

  // Propagate
  MatX6 rv_case1_op;
  if (config.recompute_part1 || !file_part1.exist("/rv_case1_op")) {
    cout << endl << "Propagating" << endl;
    rv_case1_op = dyn_3body_circ.Propagate(rv0_op, t0, tfs, true);
    dump(file_part1, "/rv_case1_op", rv_case1_op.cast<double>(), DumpMode::Overwrite);
  } else {
    rv_case1_op = load<MatX6d>(file_part1, "/rv_case1_op");
    cout << endl << "Loaded from file" << endl;
  }
  MatX6 coe_case1_op = Cart2Classical(rv_case1_op, GM_MOON);

  int n_plot = 20;
  VecX tspan_plot = arange(0, sat_period + dt_prop, dt_prop);
  if (config.plot_case1) {
    // Plot orbital elements
    fig = figure(true);
    title("Case 1: Orbital elements in OP frame");
    hold(true);
    VecX x = tspan / SECS_DAY;
    VecX y;
    // Eccentricity
    subplot(1, 2, 0);
    y = coe_case1_op.col(1);
    Plot(x, y);
    xlabel("Time [days]");
    ylabel("e [-]");
    xlim({0, dt_total.val() / SECS_DAY});
    grid(true);
    // Inclination
    subplot(1, 2, 1);
    y = coe_case1_op.col(2) * DEG;
    Plot(x, y);
    xlabel("Time [days]");
    ylabel("i [deg]");
    xlim({0, dt_total.val() / SECS_DAY});
    grid(true);
    fig->draw();

    // Plot orbits
    vector<MatX6> rv_plot_op;
    for (int i = 0; i < n_plot; i++) {
      int j = i * n_steps / n_plot;
      VecX tfs_plot = tfs[j] + tspan_plot.array();
      Vec6 rv0_plot = rv_case1_op.row(j);
      rv_plot_op.push_back(dyn_3body_circ.Propagate(rv0_plot, tfs[j], tfs_plot));
    }
    fig = figure(true);
    title("Case 1: Satellite orbit in OP frame");
    hold(true);
    PlotBody(MOON);
    for (int i = 0; i < n_plot; i++) {
      p = Plot3(rv_plot_op[i].col(0), rv_plot_op[i].col(1), rv_plot_op[i].col(2), "b-");
      p->line_width(2);
    }
    SetLim(12e3);
    fig->draw();
  }

  // **************************************************************************
  // Case 2
  // **************************************************************************
  cout << endl << endl << "*********** Case 2 ***********" << endl;

  // Dynamics
  NBodyDynamics dyn_3body(IntegratorType::RK4);
  dyn_3body.AddBody(Body::Moon());
  dyn_3body.AddBody(Body::Earth());
  dyn_3body.SetTimeStep(dt_prop);
  dyn_3body.SetFrame(MOON_CI);

  // Propagate
  MatX6 rv_case2_ci;
  if (config.recompute_part1 || !file_part1.exist("/rv_case2_ci")) {
    cout << endl << "Propagating" << endl;
    rv_case2_ci = dyn_3body.Propagate(rv0_ci, t0, tfs, true);
    dump(file_part1, "/rv_case2_ci", rv_case2_ci.cast<double>(), DumpMode::Overwrite);
  } else {
    rv_case2_ci = load<MatX6d>(file_part1, "/rv_case2_ci");
    cout << endl << "Loaded from file" << endl;
  }
  MatX6 rv_case2_op = ConvertFrame(tfs, rv_case2_ci, Frame::MOON_CI, Frame::MOON_OP);
  MatX6 coe_case2_op = Cart2Classical(rv_case2_op, GM_MOON);

  // **************************************************************************
  // Case 3
  // **************************************************************************
  cout << endl << endl << "*********** Case 3 ***********" << endl;

  NBodyDynamics dyn_nbody(IntegratorType::RK4);
  dyn_nbody.AddBody(Body::Moon(7, 1));
  dyn_nbody.AddBody(Body::Earth());
  dyn_nbody.AddBody(Body::Sun());
  dyn_nbody.SetTimeStep(dt_prop);
  dyn_nbody.SetFrame(MOON_CI);

  // Propagate
  MatX6 rv_case3_ci;
  if (config.recompute_part1 || !file_part1.exist("/rv_case3_ci")) {
    cout << endl << "Propagating" << endl;
    rv_case3_ci = dyn_nbody.Propagate(rv0_ci, t0, tfs, true);
    dump(file_part1, "/rv_case3_ci", rv_case3_ci.cast<double>(), DumpMode::Overwrite);
  } else {
    rv_case3_ci = load<MatX6d>(file_part1, "/rv_case3_ci");
    cout << endl << "Loaded from file" << endl;
  }
  MatX6 rv_case3_op = ConvertFrame(tfs, rv_case3_ci, Frame::MOON_CI, Frame::MOON_OP);
  MatX6 coe_case3_op = Cart2Classical(rv_case3_op, GM_MOON);

  // **************************************************************************
  // Case 4
  // **************************************************************************
  cout << endl << endl << "*********** Case 4 ***********" << endl;
  MatX6 rv_case4_me = ConvertFrame(tfs, rv_case3_ci, Frame::MOON_CI, Frame::MOON_ME, true);
  MatX6 coe_case4_me = Cart2Classical(rv_case4_me, GM_MOON);

  // **************************************************************************
  // Case 5
  // **************************************************************************

  NBodyDynamics dyn_nbody50(IntegratorType::RK4);
  dyn_nbody50.AddBody(Body::Moon(10, 10));
  dyn_nbody50.AddBody(Body::Earth());
  dyn_nbody50.AddBody(Body::Sun());
  dyn_nbody50.SetTimeStep(dt_prop);
  dyn_nbody50.SetFrame(MOON_CI);

  // Propagate
  MatX6 rv_case5_ci;
  VecX tfs_ = tfs.head(1000);
  if (config.recompute_part1 || !file_part1.exist("/rv_case5_ci")) {
    cout << endl << "Propagating" << endl;
    rv_case5_ci = dyn_nbody50.Propagate(rv0_ci, t0, tfs_, true);
    dump(file_part1, "/rv_case5_ci", rv_case5_ci.cast<double>(), DumpMode::Overwrite);
  } else {
    rv_case5_ci = load<MatX6d>(file_part1, "/rv_case5_ci");
    cout << endl << "Loaded from file" << endl;
  }
  MatX6 rv_case5_op = ConvertFrame(tfs_, rv_case5_ci, Frame::MOON_CI, Frame::MOON_OP, true);
  MatX6 coe_case5_op = Cart2Classical(rv_case5_op, GM_MOON);

  // **************************************************************************
  // Case 6
  // **************************************************************************
  cout << endl << endl << "*********** Case 6 ***********" << endl;

  // Time
  dt_total = 10 * DAYS_YEAR * SECS_DAY;            // [s] Total propagation time
  tspan = arange(0, dt_total + dt_step, dt_step);  // [s] Time span
  tfs = t0 + tspan.array();                        // [s] Final times
  n_steps = tspan.size();
  cout << "Total duration   " << dt_total / SECS_DAY << " days" << endl;
  cout << "Time step        " << dt_step / SECS_MINUTE << " minutes" << endl;
  cout << "Propagation step " << dt_prop << " seconds" << endl;
  cout << "Start epoch      " << Time2GregorianString(t0) << endl;
  cout << "End epoch        " << Time2GregorianString(t0 + dt_total) << endl;
  cout << "Number of steps  " << n_steps << endl;

  // Propagate
  MatX6 rv_case6_ci;
  if (config.recompute_part1 || !file_part1.exist("/rv_case6_ci")) {
    cout << endl << "Propagating" << endl;
    rv_case6_ci = dyn_nbody.Propagate(rv0_ci, t0, tfs, true);
    dump(file_part1, "/rv_case6_ci", rv_case6_ci.cast<double>(), DumpMode::Overwrite);
  } else {
    rv_case6_ci = load<MatX6d>(file_part1, "/rv_case6_ci");
    cout << endl << "Loaded from file" << endl;
  }
  MatX6 rv_case6_op = ConvertFrame(tfs, rv_case6_ci, Frame::MOON_CI, Frame::MOON_OP, true);
  MatX6 coe_case6_op = Cart2Classical(rv_case6_op, GM_MOON);

  // **************************************************************************
  // e-w plots
  // **************************************************************************

  std::vector<MatX6> coe_cases
      = {coe_case1_op, coe_case2_op, coe_case3_op, coe_case4_me, coe_case5_op, coe_case6_op};

  // Plot
  fig = figure(true);
  title("e-w plots");
  for (size_t i = 0; i < coe_cases.size(); i++) {
    fig->add_subplot(2, 2, i);
    hold(true);
    VecX e_vec = coe_cases[i].col(1);
    VecX w_vec = coe_cases[i].col(4) * DEG;
    Plot(w_vec, e_vec, "b");
    xlabel("w [deg]");
    ylabel("e [-]");
    grid(true);
    xlim({70, 110});
    ylim({0.5, 0.75});
    title("Case " + std::to_string(i + 1));
  }
  fig->draw();

  // **************************************************************************
  // Constellation stability
  // **************************************************************************
  cout << endl << endl << "*********** Constellation stability ***********" << endl;

  auto open_mode_part2 = (config.recompute_part2) ? File::Truncate : File::OpenOrCreate;
  File file_part2(output_path / "data_part2.h5", open_mode_part1);

  array<Vec6, 3> coes0_op = {
      Vec6(6541.4, 0.6, 56.2 * DEG, 0, 90 * DEG, 0),
      Vec6(6541.4, 0.6, 56.2 * DEG, 0, 90 * DEG, 120 * DEG),
      Vec6(6541.4, 0.6, 56.2 * DEG, 0, 90 * DEG, 240 * DEG),
  };
  array<MatX6, 3> rvs_ci;
  array<MatX6, 3> coes_ci;
  for (int i = 0; i < 3; i++) {
    if (config.recompute_part1 || !file_part2.exist("/rvs_ci" + std::to_string(i))
        || !file_part2.exist("/coes_ci" + std::to_string(i))) {
      cout << endl << "Propagating" << endl;

      // Initial state
      Vec6 coe0_op_ = coes0_op[i];
      Vec6 rv0_op_ = Classical2Cart(coe0_op_, GM_MOON);
      Vec6 rv0_ci_ = ConvertFrame(t0, rv0_op_, MOON_OP, MOON_CI);

      NBodyDynamics dyn_nbody_(IntegratorType::RK4);
      dyn_nbody_.AddBody(Body::Moon(7, 1));
      dyn_nbody_.AddBody(Body::Earth());
      dyn_nbody_.AddBody(Body::Sun());
      dyn_nbody_.SetTimeStep(dt_prop);
      dyn_nbody_.SetFrame(MOON_CI);

      // Propagate
      VecX tfs_ = tfs.head(1000);
      MatX6 rv_ci_ = dyn_nbody_.Propagate(rv0_ci_, t0, tfs, true);
      MatX6 coe_ci_ = Cart2Classical(rv_ci_, GM_MOON);

      rvs_ci[i] = rv_ci_;
      coes_ci[i] = coe_ci_;

      dump(file_part2, "/rvs_ci" + std::to_string(i), rvs_ci[i].cast<double>(),
           DumpMode::Overwrite);
      dump(file_part2, "/coes_ci" + std::to_string(i), coes_ci[i].cast<double>(),
           DumpMode::Overwrite);
    } else {
      rvs_ci[i] = load<MatX6d>(file_part2, "/rvs_ci" + std::to_string(i));
      coes_ci[i] = load<MatX6d>(file_part2, "/coes_ci" + std::to_string(i));
      cout << endl << "Loaded from file" << endl;
    }
  }

  auto end = GetSystemTime();
  cout << "Total elapsed time: " << PrintDuration(end - begin) << endl;

  show();
  return 0;
}

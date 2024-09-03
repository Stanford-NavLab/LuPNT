#include <lupnt/lupnt.h>
#include <matplot/matplot.h>
#include <omp.h>

#include <filesystem>
#include <highfive/H5Easy.hpp>

using namespace lupnt;
using namespace std;
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
  bool recompute_part1 = false;
  bool recompute_part2 = false;
  bool plot_case0 = false;
  bool plot_case1 = false;
  bool plot_case2 = false;
  bool plot_case3 = false;
  bool plot_case4 = false;
  bool plot_ew = false;
  bool plot_delta_M = false;
} config;

int main() {
  auto begin = GetSystemTime();

  auto output_path = GetOutputPath("ex_frozen_orbits");
  cout << "Output path: " << output_path << endl;

  auto path1 = output_path / "data_part1.h5";
  if (config.recompute_part1 && filesystem::exists(path1)) filesystem::remove(path1);
  auto open_mode1 = (config.recompute_part1 || !filesystem::exists(path1))
                        ? H5Easy::File::OpenOrCreate
                        : H5Easy::File::ReadOnly;
  H5Easy::File file_part1(path1, open_mode1);

  // Time
  Real t0 = Gregorian2Time(2009, 7, 15, 1, 0, 0);  // [s] Start time (TAI)

  // Orbital elements
  Real a = 6541.4;      // [km] Semi-major axis
  Real e = 0.60;        // [-] Eccentricity
  Real i = 56.2 * RAD;  // [deg] Inclination
  Real O = 0.00 * RAD;  // [deg] Right ascension of the ascending node
  Real w = 90.0 * RAD;  // [deg] Argument of perigee
  Real M = 0.00 * RAD;  // [deg] Mean anomaly

  Vec6 coe0_op(a, e, i, O, w, M);

  Real sat_period = GetOrbitalPeriod(a, GM_MOON);
  cout << "Satellite period: " << sat_period / SECS_HOUR << " minutes" << endl;

  // **************************************************************************
  // Case 0
  // **************************************************************************
  cout << endl << "*********** Case 0 ***********" << endl;

  // Time
  Real dt_total = sat_period;                           // [s] Total propagation time
  Real dt_step = 5 * SECS_MINUTE;                       // [s] Time step
  Real dt_prop = 20;                                    // [s] Propagation time step
  VecX tspan = arange(0, dt_total + dt_step, dt_step);  // [s] Time span
  VecX tfs = t0 + tspan.array();                        // [s] Final times
  int n_steps = tspan.size();
  cout << "Total duration    " << dt_total / SECS_HOUR << " hours" << endl;
  cout << "Time step         " << dt_step / SECS_MINUTE << " minutes" << endl;
  cout << "Propagation step  " << dt_prop << " seconds" << endl;
  cout << "Start epoch       " << Time2GregorianString(t0) << endl;
  cout << "End epoch         " << Time2GregorianString(t0 + dt_total) << endl;
  cout << "Number of steps   " << n_steps << endl;

  // Initial state
  Vec6 rv0_op = Classical2Cart(coe0_op, GM_MOON);
  Vec6 rv0_ci = ConvertFrame(t0, rv0_op, Frame::MOON_OP, Frame::MOON_CI);

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
    PlotBody(NaifId::MOON);
    p = Plot3(rv_case0_ci.col(0), rv_case0_ci.col(1), rv_case0_ci.col(2), "b-");
    p->line_width(2);
    p->marker_indices({0});
    SetLim(12e3);
    fig->draw();
  }

  // **************************************************************************
  // Case 1
  // **************************************************************************
  cout << endl << "*********** Case 1 ***********" << endl;

  Real moon_period = GetOrbitalPeriod(D_EARTH_MOON, GM_EARTH);  // [s] Moon period
  cout << "Moon period: " << moon_period / SECS_DAY << " days" << endl;

  // Time
  dt_total = 2 * DAYS_YEAR * SECS_DAY;             // [s] Total propagation time
  dt_step = 60 * SECS_MINUTE;                      // [s] Time step
  dt_prop = 1 * SECS_MINUTE;                       // [s] Propagation time step
  tspan = arange(0, dt_total + dt_step, dt_step);  // [s] Time span
  tfs = t0 + tspan.array();                        // [s] Final times
  n_steps = tspan.size();
  cout << "Total duration    " << dt_total / SECS_DAY << " days" << endl;
  cout << "Time step         " << dt_step / SECS_MINUTE << " minutes" << endl;
  cout << "Propagation step  " << dt_prop << " seconds" << endl;
  cout << "Start epoch       " << Time2GregorianString(t0) << endl;
  cout << "End epoch         " << Time2GregorianString(t0 + dt_total) << endl;
  cout << "Number of steps   " << n_steps << endl;

  // Initial state
  Vec6 rv0_moon_op = GetBodyPosVel(t0, NaifId::EARTH, NaifId::MOON, Frame::MOON_OP);
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
    PlotBody(NaifId::EARTH);
    PlotBody(NaifId::MOON, rv0_moon_op.head(3));
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
    rv_case1_op = dyn_3body_circ.Propagate(rv0_op, t0, tfs, true);
    dump(file_part1, "/rv_case1_op", rv_case1_op.cast<double>(), H5Easy::DumpMode::Overwrite);
  } else {
    rv_case1_op = H5Easy::load<MatX6d>(file_part1, "/rv_case1_op");
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
    for (int i = 0; i < n_plot; ++i) {
      int j = i * n_steps / n_plot;
      VecX tfs_plot = tfs[j] + tspan_plot.array();
      Vec6 rv0_plot = rv_case1_op.row(j);
      rv_plot_op.push_back(dyn_3body_circ.Propagate(rv0_plot, tfs[j], tfs_plot));
    }
    fig = figure(true);
    title("Case 1: Satellite orbit in OP frame");
    hold(true);
    PlotBody(NaifId::MOON);
    for (int i = 0; i < n_plot; ++i) {
      p = Plot3(rv_plot_op[i].col(0), rv_plot_op[i].col(1), rv_plot_op[i].col(2), "b-");
      p->line_width(2);
    }
    SetLim(12e3);
    fig->draw();
  }

  // **************************************************************************
  // Case 2
  // **************************************************************************
  cout << endl << "*********** Case 2 ***********" << endl;

  // Dynamics
  NBodyDynamics dyn_3body(IntegratorType::RK4);
  dyn_3body.AddBody(Body::Moon());
  dyn_3body.AddBody(Body::Earth());
  dyn_3body.SetTimeStep(dt_prop);
  dyn_3body.SetFrame(Frame::MOON_CI);

  // Propagate
  MatX6 rv_case2_ci;
  if (config.recompute_part1 || !file_part1.exist("/rv_case2_ci")) {
    rv_case2_ci = dyn_3body.Propagate(rv0_ci, t0, tfs, true);
    H5Easy::dump(file_part1, "/rv_case2_ci", rv_case2_ci.cast<double>(),
                 H5Easy::DumpMode::Overwrite);
  } else {
    rv_case2_ci = H5Easy::load<MatX6d>(file_part1, "/rv_case2_ci");
  }

  // Convert
  MatX6 rv_case2_op = ConvertFrame(tfs, rv_case2_ci, Frame::MOON_CI, Frame::MOON_OP);
  MatX6 coe_case2_op = Cart2Classical(rv_case2_op, GM_MOON);

  // **************************************************************************
  // Case 3
  // **************************************************************************
  cout << endl << "*********** Case 3 ***********" << endl;

  NBodyDynamics dyn_nbody(IntegratorType::RK4);
  dyn_nbody.AddBody(Body::Moon(7, 1));
  dyn_nbody.AddBody(Body::Earth());
  dyn_nbody.AddBody(Body::Sun());
  dyn_nbody.SetTimeStep(dt_prop);
  dyn_nbody.SetFrame(Frame::MOON_CI);

  // Propagate
  MatX6 rv_case3_ci;
  if (config.recompute_part1 || !file_part1.exist("/rv_case3_ci")) {
    rv_case3_ci = dyn_nbody.Propagate(rv0_ci, t0, tfs, true);
    H5Easy::dump(file_part1, "/rv_case3_ci", rv_case3_ci.cast<double>(),
                 H5Easy::DumpMode::Overwrite);
  } else {
    rv_case3_ci = H5Easy::load<MatX6d>(file_part1, "/rv_case3_ci");
  }

  MatX6 rv_case3_op = ConvertFrame(tfs, rv_case3_ci, Frame::MOON_CI, Frame::MOON_OP);
  MatX6 coe_case3_op = Cart2Classical(rv_case3_op, GM_MOON);

  // **************************************************************************
  // Case 4
  // **************************************************************************
  cout << endl << "*********** Case 4 ***********" << endl;
  MatX6 rv_case4_me = ConvertFrame(tfs, rv_case3_ci, Frame::MOON_CI, Frame::MOON_ME, true);
  MatX6 coe_case4_me = Cart2Classical(rv_case4_me, GM_MOON);

  // **************************************************************************
  // Case 5
  // **************************************************************************
  cout << endl << "*********** Case 5 ***********" << endl;

  NBodyDynamics dyn_nbody50(IntegratorType::RK4);
  dyn_nbody50.AddBody(Body::Moon(50, 50));
  dyn_nbody50.AddBody(Body::Earth());
  dyn_nbody50.AddBody(Body::Sun());
  dyn_nbody50.SetTimeStep(dt_prop);
  dyn_nbody50.SetFrame(Frame::MOON_CI);

  // Propagate
  MatX6 rv_case5_ci;
  if (config.recompute_part1 || !file_part1.exist("/rv_case5_ci")) {
    rv_case5_ci = dyn_nbody50.Propagate(rv0_ci, t0, tfs, true);
    H5Easy::dump(file_part1, "/rv_case5_ci", rv_case5_ci.cast<double>(),
                 H5Easy::DumpMode::Overwrite);
  } else {
    rv_case5_ci = H5Easy::load<MatX6d>(file_part1, "/rv_case5_ci");
  }

  MatX6 rv_case5_me = ConvertFrame(tfs, rv_case5_ci, Frame::MOON_CI, Frame::MOON_ME, true);
  MatX6 coe_case5_me = Cart2Classical(rv_case5_me, GM_MOON);

  // **************************************************************************
  // Case 6
  // **************************************************************************
  cout << endl << "*********** Case 6 ***********" << endl;

  // Time
  dt_total = 10 * DAYS_YEAR * SECS_DAY;            // [s] Total propagation time
  dt_step = 60 * SECS_MINUTE;                      // [s] Time step
  dt_prop = 1 * SECS_MINUTE;                       // [s] Propagation time step
  tspan = arange(0, dt_total + dt_step, dt_step);  // [s] Time span
  tfs = t0 + tspan.array();                        // [s] Final times
  n_steps = tspan.size();
  cout << "Total duration    " << dt_total / SECS_DAY << " days" << endl;
  cout << "Time step         " << dt_step / SECS_MINUTE << " minutes" << endl;
  cout << "Propagation step  " << dt_prop << " seconds" << endl;
  cout << "Start epoch       " << Time2GregorianString(t0) << endl;
  cout << "End epoch         " << Time2GregorianString(t0 + dt_total) << endl;
  cout << "Number of steps   " << n_steps << endl;

  // Propagate
  MatX6 rv_case6_ci;
  if (config.recompute_part1 || !file_part1.exist("/rv_case6_ci")) {
    rv_case6_ci = dyn_nbody.Propagate(rv0_ci, t0, tfs, true);
    H5Easy::dump(file_part1, "/rv_case6_ci", rv_case6_ci.cast<double>(),
                 H5Easy::DumpMode::Overwrite);
  } else {
    rv_case6_ci = H5Easy::load<MatX6d>(file_part1, "/rv_case6_ci");
  }
  MatX6 rv_case6_me = ConvertFrame(tfs, rv_case6_ci, Frame::MOON_CI, Frame::MOON_ME, true);
  MatX6 coe_case6_me = Cart2Classical(rv_case6_me, GM_MOON);

  // **************************************************************************
  // e-w plots
  // **************************************************************************

  std::vector<MatX6> coe_cases
      = {coe_case1_op, coe_case2_op, coe_case3_op, coe_case4_me, coe_case5_me, coe_case6_me};

  if (config.plot_ew) {
    // Plot
    fig = figure(true);
    title("e-w plots");
    for (size_t i = 0; i < coe_cases.size(); ++i) {
      fig->add_subplot(3, 2, i);
      hold(true);
      VecX e_vec = coe_cases[i].col(1);
      VecX w_vec = coe_cases[i].col(4) * DEG;
      Plot(w_vec, e_vec, "b");
      xlabel("w [deg]");
      ylabel("e [-]");
      grid(true);
      xlim({70, 110});
      ylim({0.5, 0.75});
      title("Case " + to_string(i + 1));
    }
    fig->draw();
  }

  // **************************************************************************
  // Constellation stability
  // **************************************************************************
  cout << endl << "*********** Constellation stability ***********" << endl;

  const int n_sat = 3;

  // Time
  dt_total = 2 * DAYS_YEAR * SECS_DAY;             // [s] Total propagation time
  dt_step = 15 * SECS_MINUTE;                      // [s] Time step
  dt_prop = 1 * SECS_MINUTE;                       // [s] Propagation time step
  tspan = arange(0, dt_total + dt_step, dt_step);  // [s] Time span
  tfs = t0 + tspan.array();                        // [s] Final times
  n_steps = tspan.size();
  cout << "Total duration    " << dt_total / SECS_DAY << " days" << endl;
  cout << "Time step         " << dt_step / SECS_MINUTE << " minutes" << endl;
  cout << "Propagation step  " << dt_prop << " seconds" << endl;
  cout << "Start epoch       " << Time2GregorianString(t0) << endl;
  cout << "End epoch         " << Time2GregorianString(t0 + dt_total) << endl;
  cout << "Number of steps   " << n_steps << endl;

  auto path2 = output_path / "data_part2.h5";
  if (config.recompute_part2 && filesystem::exists(path2)) filesystem::remove(path2);
  auto open_mode2 = (config.recompute_part2 || !filesystem::exists(path2))
                        ? H5Easy::File::OpenOrCreate
                        : H5Easy::File::ReadOnly;
  H5Easy::File file_part2(path2, open_mode2);

  vector<Vec6> coes0_op = {
      Vec6(6541.4, 0.6, 56.2 * RAD, 0, 90 * RAD, 0),
      Vec6(6541.4, 0.6, 56.2 * RAD, 0, 90 * RAD, 120 * RAD),
      Vec6(6541.4, 0.6, 56.2 * RAD, 0, 90 * RAD, 240 * RAD),
  };

  vector<MatX6> rvs_ci;
  for (int i = 0; i < n_sat; ++i) {
    if (config.recompute_part2 || !file_part2.exist("/rvs_ci" + to_string(i))) {
      Vec6 coe0_op_ = coes0_op[i];
      Vec6 rv0_op_ = Classical2Cart(coe0_op_, GM_MOON);
      Vec6 rv0_ci_ = ConvertFrame(t0, rv0_op_, Frame::MOON_OP, Frame::MOON_CI);
      rvs_ci.push_back(dyn_nbody.Propagate(rv0_ci_, t0, tfs, true));
      H5Easy::dump(file_part2, "/rvs_ci" + to_string(i), rvs_ci[i].cast<double>(),
                   H5Easy::DumpMode::Overwrite);
    } else {
      rvs_ci.push_back(H5Easy::load<MatX6d>(file_part2, "/rvs_ci" + to_string(i)));
    }
  }

  vector<MatX6> coes_ci;
  for (int i = 0; i < n_sat; ++i) coes_ci.push_back(Cart2Classical(rvs_ci[i], GM_MOON));

  if (config.plot_delta_M) {
    fig = figure(true);
    title("Constellation stability");
    for (int i = 0; i < 2; ++i) {
      fig->add_subplot(1, 2, i);
      hold(true);
      VecX delta_M_tmp = Wrap2Pi(coes_ci[i + 1].col(5) - coes_ci[0].col(5)) * DEG;
      VecX delta_M(delta_M_tmp.size());
      delta_M[0] = delta_M_tmp[0];
      for (int j = 1; j < delta_M.size(); ++j)
        delta_M[j] = delta_M[j - 1] + Wrap2Pi(delta_M_tmp[j] - delta_M_tmp[j - 1]);
      Plot(tspan / SECS_DAY, delta_M);
      xlabel("Time [days]");
      ylabel("M_" + to_string(i + 1) + " - M_0 [deg]");
      xlim({0, dt_total.val() / SECS_DAY});
      grid(true);
    }
    fig->draw();
  }

  // **************************************************************************
  // Constellation stability (phasing adjusted)
  // **************************************************************************
  cout << endl
       << endl
       << "*********** Constellation stability (phasing adjusted) ***********" << endl;

  vector<Vec6> coes0_op_adjusted = {
      Vec6(6541.400000, 0.6, 56.2 * RAD, 0, 90 * RAD, 0),
      Vec6(6541.623458, 0.6, 56.2 * RAD, 0, 90 * RAD, 120 * RAD),
      Vec6(6539.069348, 0.6, 56.2 * RAD, 0, 90 * RAD, 240 * RAD),
  };

  vector<MatX6> rvs_ci_adjusted;
  for (int i = 0; i < n_sat; ++i) {
    if (config.recompute_part2 || !file_part2.exist("/rvs_ci_adjusted" + to_string(i))) {
      Vec6 coe0_op_ = coes0_op_adjusted[i];
      Vec6 rv0_op_i = Classical2Cart(coe0_op_, GM_MOON);
      Vec6 rv0_ci_i = ConvertFrame(t0, rv0_op_i, Frame::MOON_OP, Frame::MOON_CI);
      rvs_ci_adjusted.push_back(dyn_nbody.Propagate(rv0_ci_i, t0, tfs, true));
      H5Easy::dump(file_part2, "/rvs_ci_adjusted" + to_string(i), rvs_ci_adjusted[i].cast<double>(),
                   H5Easy::DumpMode::Overwrite);
    } else {
      rvs_ci_adjusted.push_back(
          H5Easy::load<MatX6d>(file_part2, "/rvs_ci_adjusted" + to_string(i)));
    }
  }

  vector<MatX6> coes_ci_adjusted;
  for (int i = 0; i < n_sat; ++i) {
    coes_ci_adjusted.push_back(Cart2Classical(rvs_ci_adjusted[i], GM_MOON));
  }

  if (config.plot_delta_M) {
    fig = figure(true);
    title("Constellation stability (adjusted)");
    for (int i = 0; i < 2; ++i) {
      fig->add_subplot(1, 2, i);
      hold(true);
      VecX delta_M = Wrap2Pi(coes_ci_adjusted[i + 1].col(5) - coes_ci_adjusted[0].col(5)) * DEG;
      if (i == 1) {
        for (int j = 0; j < delta_M.size(); ++j) {
          if (delta_M[j] > 0) delta_M[j] -= 360;
        }
      }
      Plot(tspan / SECS_DAY, delta_M);
      xlabel("Time [days]");
      ylabel("M_" + to_string(i + 1) + " - M_0 [deg]");
      xlim({0, dt_total.val() / SECS_DAY});
      if (i == 0)
        ylim({110, 130});
      else
        ylim({-130, -110});
      grid(true);
    }
  }

  // **************************************************************************
  // South Pole Coverage
  // **************************************************************************
  cout << endl << "*********** Coverage ***********" << endl;

  Real min_elevation = 10 * RAD;  // [rad] Minimum elevation

  // Satellites in ME frame
  vector<MatX3> rs_me;
  for (int i = 0; i < n_sat; ++i)
    rs_me.push_back(
        ConvertFrame(tfs, rvs_ci_adjusted[i], Frame::MOON_CI, Frame::MOON_ME, true).leftCols(3));

  // Elevation
  Vec3 r_south_pole = LatLonAlt2Cart(Vec3(-90 * RAD, 0, 0), R_MOON);
  MatXi visibility = MatXi::Zero(n_sat, n_steps);
  for (int i = 0; i < n_sat; ++i) {
    VecX elev = Cart2AzElRange(rs_me[i], r_south_pole).col(1);
    for (int j = 0; j < n_steps; ++j)
      if (elev[j] >= min_elevation) visibility(i, j) = 1;
  }

  // Metrics
  Vec3 mean_pass, mean_gap, coverage;
  for (int i = 0; i < n_sat; ++i) {
    VecXi visibility_i = VecXi::Zero(n_steps + 2);
    visibility_i.segment(1, n_steps) = visibility.row(i);
    int n_passes = 0;
    for (int j = 1; j < visibility_i.size(); ++j) {
      if (visibility_i[j] - visibility_i[j - 1] > 0) {
        ++n_passes;
      }
    }
    coverage[i] = 100. * visibility.row(i).sum() / n_steps;
    mean_pass[i] = coverage[i] / 100. * dt_total / n_passes / SECS_HOUR;
    mean_gap[i] = (1 - coverage[i] / 100.) * dt_total / (n_passes - 1) / SECS_HOUR;
  }
  Real one_fold_coverage
      = 100. * (visibility.colwise().sum().array() >= 1).cast<int>().sum() / n_steps;
  Real two_fold_coverage
      = 100. * (visibility.colwise().sum().array() >= 2).cast<int>().sum() / n_steps;

  cout << "Mean pass duration  " << mean_pass.transpose() << " h" << endl;
  cout << "Mean gap duration   " << mean_gap.transpose() << " h" << endl;
  cout << "Coverage            " << coverage.transpose() << " %" << endl;
  cout << "One-fold coverage   " << one_fold_coverage << " %" << endl;
  cout << "Two-fold coverage   " << two_fold_coverage << " %" << endl;

  auto end = GetSystemTime();
  cout << "Total elapsed time: " << PrintDuration(end - begin) << endl;
  show();

  return 0;
}

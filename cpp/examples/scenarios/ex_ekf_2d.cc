// End-to-end 2D EKF scenario with synthetic range measurements and plot output.
#include <lupnt/lupnt.h>

#include "lupnt/core/constants.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/states/state.h"

using namespace lupnt;
using namespace matplot;

int main(int argc, char** argv) {
  //   RandomEngine::SetSeed(0);
  RandomEngine::SetSeed(112323);
  std::cout << SampleNormal(0.0, 1.0) << std::endl;
  std::cout << SampleNormal(0.0, 1.0) << std::endl;
  RandomEngine::SetSeed(112323);
  std::cout << SampleNormal(0.0, 1.0) << std::endl;

  // Parameters
  Real radius = 1.0e3;        // [m]
  Real tf = 1.0 * SECS_HOUR;  // [s]
  Real dt = 1.0;              // [s]

  Real v = 2.0 * PI * radius / tf;  // [m/s]
  Real w = v / radius;              // [rad/s]
  std::cout << "v: " << v << " [m/s]" << std::endl;
  std::cout << "w: " << w * DEG << " [deg/s]" << std::endl;

  // Error
  Real sigma_r = 2.0 * v;      // [m/s]
  Real sigma_theta = 5.0 * w;  // [rad/s]
  Mat3d Q = Vec3d{sigma_r, sigma_r, sigma_theta}.array().square().matrix().asDiagonal();

  // Control
  State control(Vec2{v, w}, {"v", "w"}, {"m/s", "rad/s"});

  // States
  SurfaceState2D state_gt({0.0, 0.0, 0.0});  // [m, m, rad]
  SurfaceState2D state_est = state_gt;
  state_est += SampleMvNormal(Vec3::Zero(), Q * dt * dt, 1).transpose();

  // Dynamics
  SurfaceDynamics2D dyn_gt, dyn_est;
  dyn_gt.SetQ(Q);

  // Propagate
  VecX times = Arange(Real(0.0), tf, dt);  // [s]
  MatX states_gt = dyn_gt.Propagate(state_gt, times, &control);

  // Filter
  EKF ekf;
  FilterDynamicsFunction f_dyn
      = [&dyn_est](const State& x, Real t0, Real tf, const State* u, MatXd* F) {
          return dyn_est.Propagate(x, t0, tf, u, F);
        };

  FilterMeasurementFunction f_meas = [&](const State& x, MatXd* H, MatXd* R) { return VecX(0); };
  ProcessNoiseFunction f_proc = [&](const State& x, Real t0, Real tf) {
    Real dt = tf - t0;
    return Q * dt * dt;
  };
  ekf.SetDynamicsFunction(f_dyn);
  ekf.SetMeasurementFunction(f_meas);
  ekf.SetProcessNoiseFunction(f_proc);
  ekf.SetState(state_est);
  ekf.SetCovariance(Q * dt * dt);

  // Propagate
  MatX states_est(times.size(), 3);
  for (int i = 0; i < times.size(); i++) {
    ekf.Predict(times(i), &control);
    states_est.row(i) = ekf.GetState().cast<double>();
  }

  // Plot
  Real dt_plt = 10.0;
  int step = dt_plt / dt;
  VecX t_plt = Subsample(times, step);
  MatX s_gt_plt = Subsample(states_gt, step, 1);
  MatX s_est_plt = Subsample(states_est, step, 1);

  // States
  auto fig = figure();
  fig->size(600, 600);
  fig->color("w");
  auto ax = subplot(2, 1, 0);
  hold(ax, on);
  ax->plot(t_plt, s_gt_plt.col(0).eval(), "r-");
  ax->plot(t_plt, s_gt_plt.col(1).eval(), "b-");
  ax->plot(t_plt, s_est_plt.col(0).eval(), "r--");
  ax->plot(t_plt, s_est_plt.col(1).eval(), "b--");
  ax->legend({"x", "y"});
  ax->xlabel("Time [s]");
  ax->ylabel("Position [m]");
  ax->grid(true);
  ax->title("Position");
  // Theta
  ax = subplot(2, 1, 1);
  hold(ax, on);
  ax->plot(t_plt, s_gt_plt.col(2).eval(), "g-");
  ax->plot(t_plt, s_est_plt.col(2).eval(), "g--");
  ax->legend({"x", "y", "theta"});
  ax->xlabel("Time [s]");
  ax->ylabel("Position [m] / Angle [rad]");
  ax->grid(true);
  ax->title("Position and Theta");
  fig->draw();

  // Trajectory
  fig = figure();
  fig->size(400, 400);
  fig->color("w");
  hold(on);
  plot(s_gt_plt.col(0).eval(), s_gt_plt.col(1).eval(), "r-");
  plot(s_est_plt.col(0).eval(), s_est_plt.col(1).eval(), "r--");
  axis("equal");
  grid(true);
  title("Trajectory");
  xlabel("x [m]");
  ylabel("y [m]");
  fig->draw();

  show();

  return 0;
}

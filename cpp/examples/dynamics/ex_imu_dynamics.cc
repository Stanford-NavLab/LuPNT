// Propagates an LN200S IMU bias model and plots gyro/accelerometer bias drift
// over one day.
#include <lupnt/lupnt.h>

#include "lupnt/core/constants.h"

using namespace lupnt;
using namespace matplot;

int main(int argc, char** argv) {
  ImuState imu_state_0;

  ImuDynamics imu_dynamics;
  imu_dynamics.SetModel(ImuModel::LN200S);

  Real tf = SECS_DAY;                      // [s]
  Real dt = 1.0;                           // [s]
  VecX times = Arange(Real(0.0), tf, dt);  // [s]
  MatX imu_states = imu_dynamics.Propagate(imu_state_0, times);

  // Units
  Real dt_plot = 60.0;  // [s]
  int step = dt_plot / dt;
  VecX times_plot = Subsample(times, step) / SECS_HOUR;
  MatX imu_states_plot = Subsample(imu_states, step, 1);
  for (int i = 0; i < 3; i++) imu_states_plot.col(i) *= DEG;
  for (int i = 3; i < 6; i++) imu_states_plot.col(i) *= M_KM;

  auto fig = figure();
  fig->size(600, 400);
  fig->color("w");

  auto ax = subplot(2, 1, 0);
  hold(ax, on);
  ax->plot(times_plot, imu_states_plot.col(0).eval(), "r-");
  ax->plot(times_plot, imu_states_plot.col(1).eval(), "b-");
  ax->plot(times_plot, imu_states_plot.col(2).eval(), "g-");
  legend(ax, {"b_w_x", "b_w_y", "b_w_z"});
  xlabel(ax, "Time [h]");
  ylabel(ax, "Angular velocity bias [deg/s]");

  ax = subplot(2, 1, 1);
  hold(ax, on);
  ax->plot(times_plot, imu_states_plot.col(3).eval(), "r-");
  ax->plot(times_plot, imu_states_plot.col(4).eval(), "b-");
  ax->plot(times_plot, imu_states_plot.col(5).eval(), "g-");
  legend(ax, {"b_a_x", "b_a_y", "b_a_z"});
  xlabel(ax, "Time [h]");
  ylabel(ax, "Acceleration bias [m/s^2]");

  show();
  return 0;
}

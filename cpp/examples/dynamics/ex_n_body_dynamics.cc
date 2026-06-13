// Propagates a lunar orbit with NBodyDynamics and then repeats propagation with
// autodiff enabled to inspect the state transition matrix.
#include <lupnt/lupnt.h>

using namespace lupnt;
using namespace matplot;

int main() {
  // Time
  Real t0 = GregorianToTime(2009, 7, 15, 1, 0, 0);  // [s] Start time (TAI)

  // Orbital elements
  Real a = 6541.4e3;    // [m] Semi-major axis
  Real e = 0.60;        // [-] Eccentricity
  Real i = 56.2 * RAD;  // [deg] Inclination
  Real O = 0.00 * RAD;  // [deg] Right ascension of the ascending node
  Real w = 90.0 * RAD;  // [deg] Argument of perigee
  Real M = 0.00 * RAD;  // [deg] Mean anomaly

  Vec6 coe0_op(a, e, i, O, w, M);
  Vec6 rv0_op = ClassicalToCart(coe0_op, GM_MOON);
  Vec6 rv0_ci = ConvertFrame(t0, rv0_op, Frame::MOON_OP, Frame::MOON_CI);

  // Time
  Real dt_total = 30 * SECS_DAY;                                // [s] Total propagation time
  Real dt_step = 5 * SECS_MINUTE;                               // [s] Time step
  Real dt_prop = 20;                                            // [s] Propagation time step
  VecX tspan = Arange(Real(0.0), dt_total + dt_step, dt_step);  // [s] Time span
  VecX tfs = t0 + tspan.array();                                // [s] Final times
  int n_steps = tspan.size();

  // Dynamics
  NBodyDynamics dyn;
  dyn.SetIntegrator(IntegratorType::RKF45);
  dyn.AddBody(Body::Moon(5, 5));
  // dyn.AddBody(Body::Earth());
  // dyn.AddBody(Body::Sun());
  dyn.SetTimeStep(dt_prop);
  dyn.SetFrame(Frame::MOON_CI);

  // Propagate
  dyn.SetAutodiff(false);
  MatX6 rv_ci = dyn.Propagate(rv0_ci, tspan);

  auto fig = figure(true);
  hold(true);
  grid(true);
  Plot3(rv_ci.col(0).eval(), rv_ci.col(1).eval(), rv_ci.col(2).eval());
  fig->show();

  // Propagate with state transition matrix.
  dyn.SetAutodiff(true);
  MatX6 rv_ci_ad(n_steps, 6);
  auto pbar = Logger::GetProgressBar(tfs.size(), "Propagating", "Propagation");
  rv_ci_ad.row(0) = rv0_ci;
  for (int i = 1; i < n_steps; i++) {
    Real t0_i = tfs(i - 1);
    Real tf_i = tfs(i);
    Vec6 rv0_i = rv_ci_ad.row(i - 1);
    MatXd stm;
    Vec6 rv_i = dyn.Propagate(rv0_i, t0_i, tf_i, nullptr, &stm);
    if (i == 1) std::cout << stm << std::endl;
    rv_ci_ad.row(i) = rv_i;
    pbar->Update(i);
  }
  pbar->Finish();
}

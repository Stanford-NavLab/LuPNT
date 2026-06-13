// Minimal Plotly helper example for producing an interactive orbit plot through
// the Python plotting bridge.
#include <lupnt/lupnt.h>
#include <matplot/matplot.h>

using namespace lupnt;

int main() {
  double t0_utc = GregorianToTime("2025-01-01T12:00:00");
  double t0_tai = UtcToTai(t0_utc);
  double tf_tai = t0_tai + 12 * SECS_HOUR;
  VecXd tspan = VecXd::LinSpaced(100, 0, tf_tai - t0_tai);

  // //
  // ************************************************************************************
  // // Moon Spacecraft

  // State coe_mop
  //     = ClassicalOE({6541.4, 0.6, 65.5 * RAD, 0 * RAD, 90 * RAD, 180 * RAD},
  //     Frame::MOON_OP);
  // State rv_mci = ConvertFrame(t0_tai, ClassicalToCart(coe_mop, GM_MOON),
  // Frame::MOON_CI);

  // auto sc_dyn = MakePtr<NBodyDynamics>();
  // sc_dyn->SetFrame(Frame::MOON_CI);
  // sc_dyn->AddBody(Body::Earth());
  // sc_dyn->AddBody(Body::Moon());
  // sc_dyn->SetTimeStep(10.0);  // [s]

  // GNSS spacecraft
  State coe = ClassicalOE({20200e3, 0.1, 55 * RAD, 0 * RAD, 90 * RAD, 180 * RAD}, Frame::ECI);
  State rv_eci = ConvertFrame(t0_tai, ClassicalToCart(coe, GM_EARTH), Frame::ECI);

  KeplerianDynamics<ClassicalOE> sc_dyn(GM_EARTH);
  MatX6 coes = sc_dyn.Propagate(coe, t0_tai + tspan.array());  // [s]
  MatX6d rvs = ClassicalToCart(coes, GM_EARTH);                // Convert to Cartesian

  auto fig = matplot::figure(true);
  fig->size(900, 500);
  matplot::hold(matplot::on);
  lupnt::Plot(tspan / SECS_DAY, rvs.col(0).eval())->display_name("x");
  lupnt::Plot(tspan / SECS_DAY, rvs.col(1).eval())->display_name("y");
  lupnt::Plot(tspan / SECS_DAY, rvs.col(2).eval())->display_name("z");
  matplot::xlabel("Time [days]");
  matplot::ylabel("Position [m]");
  matplot::legend();
  matplot::show();
  return 0;
}

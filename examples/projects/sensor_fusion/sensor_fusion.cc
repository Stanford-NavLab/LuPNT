#include <lupnt/lupnt.h>
#include <omp.h>

#include <cmath>

using namespace lupnt;
using namespace matplot;

int main() {
  // Parameters
  int N_sat = 3;
  int N_planes = 3;
  real tai0 = StringToTAI("2030/01/01 12:00:00.00 UTC");

  // Orbital Elements
  Vector3 sma{6142.4, 6142.4, 6142.4};
  Vector3 ecc{0.001, 0.6, 0.6};
  Vector3 imc{0.001, 57.7, 57.7};
  Vector3 raan{0, -90, 0};
  Vector3 aop{-45, -90, 90};

  real orbit_period = 2 * M_PI * pow(sma[0], 1.5) / sqrt(GM_MOON);
  int N_sat_plane = N_sat / N_planes;

  // Constellation
  std::vector<Vector6> rv0_mci(N_sat);
  for (int i = 0; i < N_planes; i++) {
    for (int j = 0; j < N_sat_plane; j++) {
      int idx = i * N_sat_plane + j;
      real ma = 360.0 / N_sat_plane * j;
      Vector6 coe0{sma[i], ecc[i], imc[i], raan[i], aop[i], ma};
      coe0.segment(2, 4) *= DEG_PER_RAD;
      Vector6 rv0_mop = ClassicalToCartesian(coe0, GM_MOON);
      rv0_mci[idx] = FrameConverter::Convert(tai0, rv0_mop, Frame::MOON_OP,
                                             Frame::MOON_CI);
    }
  }

  // Dynamics
  NBodyDynamics dyn_true;
  dyn_true.SetPrimaryBody(Body::Moon(5, 5));
  dyn_true.AddBody(Body::Earth());
  dyn_true.SetTimeStep(1.0);

  std::vector<NBodyDynamics> dyns(N_sat);
  for (int i = 0; i < N_sat; i++) {
    dyns[i].SetPrimaryBody(Body::Moon(5, 5));
    dyns[i].AddBody(Body::Earth());
    dyns[i].SetTimeStep(1.0);
  }

  // Propagate
  real Dt = 15.0;
  real tf = 1.0 * orbit_period;
  std::cout << "tf = " << tf / SECS_PER_HOUR << " h" << std::endl;
  int N_steps = (tf / Dt).val();
  VectorX tspan = VectorX::LinSpaced(N_steps, 0, tf);
  VectorX tai = tai0 + tspan.array();

  std::vector<MatrixX> rv_hist(N_sat);
  // compute time
  auto start = omp_get_wtime();
  // Print number of threads
  std::cout << "Number of threads: " << omp_get_max_threads() << std::endl;
#pragma omp parallel for
  for (int i = 0; i < N_sat; i++) {
    std::cout << "Thread " << i << " start." << std::endl;
    NBodyDynamics dyn_true;
    dyn_true.SetPrimaryBody(Body::Moon(5, 5));
    dyn_true.AddBody(Body::Earth());
    dyn_true.SetTimeStep(1.0);
    dyn_true.Propagate(rv0_mci[i], tai0, tai, Dt, true);
    std::cout << "Thread " << i << " end." << std::endl;
    // rv_hist[i] = dyns[i].Propagate(rv0_mci[i], tai0, tai, Dt, true);
  }
  auto end = omp_get_wtime();
  std::cout << "Time: " << end - start << " s" << std::endl;
  return 0;

  // Plotting
  auto fig = figure(true);
  title("Satellite Trajectories");
  xlabel("X (km)");
  ylabel("Y (km)");
  zlabel("Z (km)");

  hold(on);
  for (int i = 0; i < N_sat; i++) {
    VectorX x = rv_hist[i].col(0);
    VectorX y = rv_hist[i].col(1);
    VectorX z = rv_hist[i].col(2);
    lupnt::plot3(x, y, z);
  }

  // Generate sphere data
  std::vector<std::vector<double>> x, y, z;
  int n = 50;  // Number of points in each direction
  for (int i = 0; i <= n; ++i) {
    std::vector<double> x_row, y_row, z_row;
    for (int j = 0; j <= n; ++j) {
      double theta = i * pi / n;
      double phi = j * 2 * pi / n;
      x_row.push_back(R_MOON * sin(theta) * cos(phi));
      y_row.push_back(R_MOON * sin(theta) * sin(phi));
      z_row.push_back(R_MOON * cos(theta));
    }
    x.push_back(x_row);
    y.push_back(y_row);
    z.push_back(z_row);
  }
  surf(x, y, z);

  double lim = 10e3;
  xlim({-lim, lim});
  ylim({-lim, lim});
  zlim({-lim, lim});

  // fig->draw();

  // Create 3 2D plots for x, y, z vs time
  auto fig2 = figure(true);
  for (int j = 0; j < 3; j++) {
    subplot(3, 1, j);
    hold(on);
    for (int i = 0; i < N_sat; i++) {
      lupnt::plot(tspan / SECS_PER_HOUR, rv_hist[i].col(j));
    }
    grid(on);
    xlabel("Time [h]");
    ylabel("Position [km]");
  }

  // show();
}
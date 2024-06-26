#include <lupnt/lupnt.h>
#include <omp.h>
#include <unistd.h>  // For sleep function

#include <iostream>

using namespace lupnt;

int main() {
  int num_threads = 5;     // Number of threads to use
  int num_iterations = 5;  // Number of iterations in the loop

  omp_set_num_threads(num_threads);

  real t0_tai = StringToTAI("2030/01/01 12:00:00.00 UTC");
  real Dt = 5.;
  real tf = 1 * 12. * SECS_PER_HOUR;
  int N_steps = (tf / Dt).val();
  VectorX tspan = VectorX::LinSpaced(N_steps, 0, tf);
  VectorX t_tai = t0_tai + tspan.array();

  VectorX x0(6);
  x0 << -4556, 4003, 933, -0.584, -0.586, -0.340;

  // std::vector<MatrixX> x(num_iterations);
  auto start = std::chrono::high_resolution_clock::now();

  // #pragma omp parallel for
  for (int i = 0; i < num_iterations; i++) {
    int thread_id = omp_get_thread_num();
    NBodyDynamics dyn_true;
    dyn_true.SetPrimaryBody(Body::Moon(5, 5));
    // dyn_true.AddBody(Body::Earth());
    dyn_true.SetTimeStep(1.0);
    // CartesianTwoBodyDynamics dyn_true(MU_MOON);
    dyn_true.Propagate(x0, t0_tai, t_tai, Dt, true);
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Total time: " << elapsed.count() << " s" << std::endl;

  return 0;
}
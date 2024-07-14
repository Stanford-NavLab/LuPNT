
#include <lupnt/dynamics/propagator.h>

using namespace lupnt;

double OMEGA = 2.0 * M_PI;  // angular frequency

ODE HarmonicOscillator = [](const Real t, const VecX& x) {
  VecX dxdt(2);
  dxdt[0] = x[1];
  dxdt[1] = -OMEGA * OMEGA * x[0];
  return dxdt;
};

int main() {
  IntegratorParams params =
      IntegratorParams(20, 1e-6, 1e-6);  // max_iter, abstol, reltol
  NumericalPropagator prop_rk4 = NumericalPropagator("RK4");
  NumericalPropagator prop_rk8 = NumericalPropagator("RK8");
  NumericalPropagator prop_rkf45 = NumericalPropagator("RKF45", params);
  NumericalPropagator prop_rkf78 = NumericalPropagator("RKF78", params);

  Real dt = 0.01;  // time step

  std::vector<NumericalPropagator*> vec_prop = {&prop_rk4, &prop_rk8,
                                                &prop_rkf45, &prop_rkf78};
  std::vector<std::string> vec_method = {"RK4", "RK8", "RKF45", "RKF78"};
  int n_methods = vec_prop.size();
  std::vector<std::chrono::duration<double>> vec_elapsed_seconds(n_methods);
  std::vector<VecX> vec_x(n_methods);
  std::vector<std::vector<Real>> vec_t_history(n_methods);

  for (int j = 0; j < n_methods; j++) {
    // Time the computation time
    vec_prop[j]->SetLogHistory(true);
    auto start = std::chrono::high_resolution_clock::now();
    VecX x(2);
    x << 1.0, 0.0;  // initial position
    x = vec_prop[j]->Propagate(HarmonicOscillator, 0.0, 10.0, x, dt);
    auto end = std::chrono::high_resolution_clock::now();
    vec_elapsed_seconds[j] = end - start;
    vec_x[j] = x;
    vec_prop[j]->GetTimeHistory(vec_t_history[j]);
  }

  // After one period, the oscillator should be back to the initial position.
  std::cout << "Method    " << "Elapsed Time (ms)    " << "Final Position"
            << std::endl;
  std::cout << "----------------------------------------------------"
            << std::endl;
  for (int j = 0; j < n_methods; j++) {
    std::cout << vec_method[j] << "      "
              << vec_elapsed_seconds[j].count() * 1e3 << "    "
              << vec_x[j].transpose() << std::endl;
  }
}

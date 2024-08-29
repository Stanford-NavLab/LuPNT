
#include <lupnt/lupnt.h>

using namespace lupnt;
using namespace matplot;
using namespace std;

double OMEGA = 2.0 * PI;  // angular frequency

ODE HarmonicOscillator = [](const Real t, const VecX& x) {
  VecX dxdt(2);
  dxdt[0] = x[1];
  dxdt[1] = -OMEGA * OMEGA * x[0];
  return dxdt;
};

// Van der Pol oscillator with period 1
ODE VanderPol = [](const Real t, const VecX& x) {
  VecX dxdt(2);
  dxdt(0) = x(1);
  dxdt(1) = 0.1 * (1 - x(0) * x(0)) * x(1) - x(0);
  return dxdt;
};

ODE VanderPolStiff = [](const Real t, const VecX& x) {
  VecX dxdt(2);
  dxdt(0) = x(1);
  dxdt(1) = 10.0 * (1 - x(0) * x(0)) * x(1) - x(0);
  return dxdt;
};

int main() {
  // list of problems
  std::vector<ODE> vec_ode = {HarmonicOscillator, VanderPol, VanderPolStiff};
  std::vector<std::string> vec_ode_name
      = {"Harmonic Oscillator", "Van der Pol (non-stiff, mu=0.1)", "Van der Pol (stiff, mu=10.0)"};
  std::vector<std::string> vec_method = {"RK4", "RK8", "RKF45"};
  std::vector<std::vector<std::chrono::duration<double>>> vv_elapsed_seconds;
  std::vector<std::vector<VecX>> vv_x;
  int n_methods = vec_method.size();
  int n_problem = vec_ode.size();

  for (int p = 0; p < n_problem; p++) {
    IntegratorParams params = IntegratorParams(20, 1e-8, 1e-8);  // max_iter, abstol, reltol
    NumericalPropagator prop_rk4 = NumericalPropagator(IntegratorType::RK4);
    NumericalPropagator prop_rk8 = NumericalPropagator(IntegratorType::RK8);
    NumericalPropagator prop_rkf45 = NumericalPropagator(IntegratorType::RKF45, params);

    Real dt = 0.01;  // time step
    Real t0 = 0.0;   // initial time
    Real tf = 10.0;  // final time

    std::vector<NumericalPropagator*> vec_prop = {&prop_rk4, &prop_rk8, &prop_rkf45};
    std::vector<std::chrono::duration<double>> vec_elapsed_seconds(n_methods);
    std::vector<VecX> vec_x(n_methods);
    std::vector<std::vector<Real>> vec_t_history(n_methods);
    std::vector<std::vector<VecX>> vec_x_history(n_methods);

    for (int j = 0; j < n_methods; j++) {
      // Time the computation time
      vec_prop[j]->SetLogHistory(true);
      auto start = std::chrono::high_resolution_clock::now();
      VecX x(2);
      x << 1.0, 0.0;  // initial position
      x = vec_prop[j]->Propagate(vec_ode[p], t0, tf, x, dt);
      auto end = std::chrono::high_resolution_clock::now();
      vec_elapsed_seconds[j] = end - start;
      vec_x[j] = x;
      vec_prop[j]->GetTimeHistory(vec_t_history[j]);
      vec_prop[j]->GetStateHistory(vec_x_history[j]);
      cout << "Finished Method: " << vec_method[j] << " - " << vec_ode_name[p] << endl;
    }

    vv_x.push_back(vec_x);
    vv_elapsed_seconds.push_back(vec_elapsed_seconds);

    // Plot the results
    bool plot_results = true;

    if (plot_results) {
      figure(true);
      for (int i = 0; i < 2; i++) {
        subplot(2, 1, i + 1);
        for (int j = 0; j < n_methods; j++) {
          int tsize = vec_t_history[j].size();
          std::vector<double> t_plot(tsize);
          std::vector<double> x_plot(tsize);
          std::vector<double> v_plot(tsize);
          for (int k = 0; k < t_plot.size(); k++) {
            t_plot[k] = vec_t_history[j][k].val();
            x_plot[k] = vec_x_history[j][k](0).val();
            v_plot[k] = vec_x_history[j][k](1).val();
          }
          if (i == 0) {
            auto sc = plot(t_plot, x_plot, "o");
            sc->display_name(vec_method[j]);
          } else {
            auto sc = plot(t_plot, v_plot, "o");
            sc->display_name(vec_method[j]);
          }
          if (j == 0) hold(on);
        }  // end for methods
        grid(on);
        xlabel("Time (s)");
        if (i == 0)
          ylabel("Position");
        else
          ylabel("Velocity");
        title("Problem: " + vec_ode_name[p]);
        matplot::legend(true);
        hold(off);
      }  // end for position and velocity
      show();
    }  // end if plot_results
  }

  // Print Results
  for (int p = 0; p < n_problem; p++) {
    cout << "Problem: " << vec_ode_name[p] << endl;
    cout << "----------------------------------------------------" << endl;
    cout << setw(10) << "Method"
         << "  Elapsed Time (ms)"
         << "    Final Position" << endl;
    cout << "----------------------------------------------------" << endl;
    for (int j = 0; j < n_methods; j++) {
      cout << setw(10) << vec_method[j] << "      " << vv_elapsed_seconds[p][j].count() * 1e3
           << "       " << vv_x[p][j].transpose() << endl;
    }
    cout << "----------------------------------------------------" << endl;
    cout << " " << endl;
  }
}

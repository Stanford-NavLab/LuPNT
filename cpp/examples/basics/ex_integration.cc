// Compares LuPNT integrators on simple oscillator problems and plots the
// resulting state histories.
#include <lupnt/lupnt.h>

using namespace lupnt;
using namespace matplot;
using namespace std;

double OMEGA = 2.0 * PI;  // angular frequency

ODE HarmonicOscillator
    = [](const Real t, const VecX& x) { return Vec2(x(1), -OMEGA * OMEGA * x(0)); };

// Van der Pol oscillator with period 1
ODE VanderPol = [](const Real t, const VecX& x) {
  return Vec2(x(1), 0.10 * (1 - x(0) * x(0)) * x(1) - x(0));
};

ODE VanderPolStiff = [](const Real t, const VecX& x) {
  return Vec2(x(1), 10.0 * (1 - x(0) * x(0)) * x(1) - x(0));
};

int main() {
  // Problems span a simple harmonic oscillator and non-stiff/stiff Van der Pol
  // variants, making it easy to see where adaptive methods help.
  std::vector<ODE> odes = {HarmonicOscillator, VanderPol, VanderPolStiff};
  std::vector<std::string> names
      = {"Harmonic Oscillator", "Van der Pol (non-stiff, mu=0.1)", "Van der Pol (stiff, mu=10.0)"};
  std::vector<std::string> methods = {"RK4", "RK8", "RKF45", "PD45"};
  int n_methods = methods.size();
  int n_problems = odes.size();

  IntegratorParams params;
  params.max_iter = 40;
  params.abstol = 1e-8;
  params.reltol = 1e-8;

  RK4 rk4;
  RK8 rk8;
  RKF45 rkf45;
  rkf45.SetParams(params);
  PD45 pd45;
  pd45.SetParams(params);

  Real dt = 0.01;
  VecX tfs = Arange(Real(0.0), Real(10.0), Real(0.01));

  std::vector<Integrator*> integrators = {&rk4, &rk8, &rkf45, &pd45};

  for (int i = 0; i < n_problems; i++) {
    VecX elapsed_time(n_methods);
    std::vector<MatX2> solutions(n_methods);

    for (int j = 0; j < n_methods; j++) {
      solutions[j].resize(tfs.size(), 2);
      auto start = std::chrono::high_resolution_clock::now();
      Vec2 x{1.0, 0.0};
      for (int k = 1; k < tfs.size(); k++) {
        // cout << odes[i](tfs[k - 1], x).transpose() << endl;
        x = integrators[j]->Step(odes[i], tfs[k - 1], x, dt);
        solutions[j].row(k) = x;
      }
      auto end = std::chrono::high_resolution_clock::now();
      elapsed_time[j] = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

      cout << "Finished Method: " << methods[j] << " - " << names[i] << endl;
    }

    cout << "Problem: " << names[i] << endl;
    cout << string(50, '-') << endl << std::right;
    cout << setw(15) << "Method";
    cout << setw(15) << "Time [us]";
    cout << setw(20) << "Final Position" << endl;
    cout << string(50, '-') << endl << std::right;
    for (int j = 0; j < n_methods; j++) {
      cout << setw(15) << methods[j];
      cout << setw(15) << elapsed_time[j];
      cout << "  " << solutions[i].row(tfs.size() - 1) << endl;
    }
    cout << string(50, '-') << endl;

    bool plot_results = true;
    if (plot_results) {
      figure(true);
      for (int k = 0; k < 2; k++) {
        subplot(2, 1, k + 1);
        hold(on);
        for (int j = 0; j < n_methods; j++)
          Plot(tfs, solutions[j].col(k))->display_name(methods[j]);
        grid(on);
        xlabel("Time (s)");
        if (i == 0)
          ylabel("Position");
        else
          ylabel("Velocity");
        title("Problem: " + names[i]);
        matplot::legend(true);
        hold(off);
      }
      show();
    }
  }
}

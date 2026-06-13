// Monte Carlo propagation of LuPNT clock models. The plots compare the
// one-sigma growth for 2-state and 3-state clock dynamics over one day.
#include <lupnt/lupnt.h>
#include <matplot/backend/backend_interface.h>
#include <matplot/freestanding/plot.h>
#include <matplot/matplot.h>

using namespace lupnt;
using namespace matplot;

int main() {
  // Hardware classes are ordered roughly from lower to higher stability.
  std::vector<ClockModel> models = {
      ClockModel::OCXO,      ClockModel::USO,  ClockModel::CSAC,
      ClockModel::MINI_RAFS, ClockModel::RAFS, ClockModel::DSAC,
  };
  std::vector<std::string> model_names = {"OCXO", "USO", "CSAC", "MiniRAFS", "RAFS", "DSAC"};

  Real tf = SECS_DAY;                           // [s]
  Real dt = 1.0;                                // [s]
  VecX times = Arange(Real(0.0), tf + dt, dt);  // [s]
  int n_t = times.size();
  int n_runs = 100;
  // For each model, for each run, store the propagated states
  std::vector<std::vector<MatX2>> x_clk_2_runs(models.size(), std::vector<MatX2>(n_runs));
  std::vector<std::vector<MatX3>> x_clk_3_runs(models.size(), std::vector<MatX3>(n_runs));

  // For each model, for each time, store the sigma (stddev) for each state variable
  std::vector<MatX2> sigma_clk_2(models.size(), MatX2::Zero(n_t, 2));
  std::vector<MatX3> sigma_clk_3(models.size(), MatX3::Zero(n_t, 3));

  auto pbar = Logger::GetProgressBar(models.size() * n_runs, "Propagating clock states");
  for (size_t i = 0; i < models.size(); ++i) {
    ClockDynamics clk_dyn;
    clk_dyn.SetModel(models[i]);
#pragma omp parallel for
    for (int run = 0; run < n_runs; ++run) {
      x_clk_2_runs[i][run] = clk_dyn.Propagate(Vec2d::Zero(), times);
      x_clk_3_runs[i][run] = clk_dyn.Propagate(Vec3d::Zero(), times);
      pbar->Update();
    }
  }

  // Compute sigma for each time step and each state variable
  for (size_t i = 0; i < models.size(); ++i) {
    // 2-state
    for (int t = 0; t < n_t; ++t) {
      std::vector<double> vals_0, vals_1;
      for (int run = 0; run < n_runs; ++run) {
        vals_0.push_back(x_clk_2_runs[i][run](t, 0));
        vals_1.push_back(x_clk_2_runs[i][run](t, 1));
      }
      sigma_clk_2[i](t, 0) = stddev(vals_0);
      sigma_clk_2[i](t, 1) = stddev(vals_1);
    }
    // 3-state
    for (int t = 0; t < n_t; ++t) {
      std::vector<double> vals_0, vals_1, vals_2;
      for (int run = 0; run < n_runs; ++run) {
        vals_0.push_back(x_clk_3_runs[i][run](t, 0));
        vals_1.push_back(x_clk_3_runs[i][run](t, 1));
        vals_2.push_back(x_clk_3_runs[i][run](t, 2));
      }
      sigma_clk_3[i](t, 0) = stddev(vals_0);
      sigma_clk_3[i](t, 1) = stddev(vals_1);
      sigma_clk_3[i](t, 2) = stddev(vals_2);
    }
  }

  std::cout << "State size: 2 (sigma at final time step)" << std::endl;
  for (size_t i = 0; i < models.size(); ++i) {
    std::cout << fmt::format("- {}: {:.2g} km, {:.2g} m/s", model_names[i],
                             sigma_clk_2[i](n_t - 1, 0) * C, sigma_clk_2[i](n_t - 1, 1) * C * M_KM);
    std::cout << std::endl;
  }
  std::cout << "State size: 3 (sigma at final time step)" << std::endl;
  for (size_t i = 0; i < models.size(); ++i) {
    std::cout << fmt::format("- {}: {:.2g} km, {:.2g} m/s, {:.2g} mm/s^2", model_names[i],
                             sigma_clk_3[i](n_t - 1, 0) * C, sigma_clk_3[i](n_t - 1, 1) * C * M_KM,
                             sigma_clk_3[i](n_t - 1, 2) * C * M_KM * MM_M);
    std::cout << std::endl;
  }

  auto fig = figure(true);
  fig->size(1200, 600);
  fig->color("w");

  Real dt_plot = 60.0;  // [s]
  int step = static_cast<int>(dt_plot / dt);
  for (int i = 0; i < 3; ++i) {
    for (int k = 0; k < 2; ++k) {
      if (i == 3 && k == 1) continue;

      auto ax = subplot(3, 2, k + i * 2);
      hold(ax, on);
      for (size_t j = 0; j < models.size(); ++j) {
        // 2-state
        if (k == 0 && i < 2) {
          std::vector<double> state_2(n_t / step);
          for (int t = 0; t < n_t; t += step) state_2[t / step] = sigma_clk_2[j](t, i) * C * M_KM;
          semilogy(ax, times, state_2)->display_name(model_names[j] + " (2-state)");
        }
        // 3-state
        if (k == 1) {
          std::vector<double> state_3(n_t / step);
          for (int t = 0; t < n_t; t += step) state_3[t / step] = sigma_clk_3[j](t, i) * C * M_KM;
          semilogy(ax, times, state_3, "--")->display_name(model_names[j] + " (3-state)");
        }
      }

      xlabel(ax, "Time [s]");
      if (i == 0) {
        ylabel(ax, "Clock Bias [m]");
      } else if (i == 1) {
        ylabel(ax, "Clock Drift [m/s]");
      } else if (i == 2) {
        ylabel(ax, "Clock Drift Rate [m/s^2]");
      }
      if (i == 0) legend(ax);
    }
  }
  show();
  return 0;
}

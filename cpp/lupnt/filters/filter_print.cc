#include "lupnt/filters/filter_print.h"

#include <iomanip>

#include "lupnt/numerics/math_utils.h"

namespace lupnt {

  // Filter Error Messages
  void PrintEKFProgressHeaderPVC() {
    std::cout << "Run Simulation" << std::endl;
    std::cout << " " << std::endl;
    std::cout << " " << std::endl;
    std::cout << "Time [min]  | Pos Err [m] | Vel Err [mm/s] | Clk Bias Err [ms]" << std::endl;
    std::cout << "--------------------------------------------------------------" << std::endl;
  }

  void PrintEKFProgressHeaderPVC(int n_sat) {
    std::cout << "Run Simulation" << std::endl;
    std::cout << " " << std::endl;
    std::cout << " " << std::endl;

    //        |                  Sat1                     |                   Sat2 | Sat3  ...
    //  Time  | PosErr [m] | VelErr [mm/s] | ClkBErr [ms] | PosErr [m] | VelErr [mm/s] | ClkBErr
    //  [ms] | ...
    // -------------------------------------------------------------------------------------------------------

    std::cout << "            |";
    // first row, plot each sat
    for (int i = 0; i < n_sat; i++) {
      std::cout << "Sat" << i + 1 << "                                             | ";
    }
    std::cout << std::endl;

    // second row, time and pos, vel and clk bias error for each sat
    std::cout << "Time [min]  |";
    for (int i = 0; i < n_sat; i++) {
      std::cout << "Pos Err [m] | Vel Err [mm/s] | Clk Bias Err [ms] | ";
    }
    std::cout << std::endl;

    for (int i = 0; i < n_sat; i++) {
      std::cout << "----------------------------------------------------------";
    }
    std::cout << std::endl;
  }

  void PrintEKFProgressPVC(double t, double x_pos_err, double x_vel_err, double x_clk_bias_err) {
    std::cout.precision(5);
    std::cout << std::left << std::setw(12) << t / 60 << " " << std::left << std::setw(12)
              << x_pos_err << "  " << std::left << std::setw(14) << x_vel_err << "   " << std::left
              << std::setw(16) << x_clk_bias_err << std::endl;
  };

  void PrintEKFProgressPVC(double t, const VecXd& est_err, int n_sat) {
    std::cout.precision(5);
    std::cout << std::left << std::setw(12) << t / 60 << " ";
    // for each satellite
    for (int i = 0; i < n_sat; i++) {
      std::cout << std::left << std::setw(12) << est_err(i * 4) << "  " << std::left
                << std::setw(14) << est_err(i * 4 + 1) << "   " << std::left << std::setw(16)
                << est_err(i * 4 + 2) << "  | ";
    }
    std::cout << std::endl;
  };

  void PrintEstimationStatistics(const VecXd& num_meas, const MatXd& error_mat, double ratio,
                                 int n_sat) {
    int n_time = num_meas.size();
    Vec4d rms, means, stds, p68, p95, p99;

    if (error_mat.rows() != 4 * n_sat) {
      std::cout << "Wrong Mat Size, Error Mat size must be (4 x timestep)" << std::endl;
      return;
    }

    // extract statistics range data
    int start_idx = (int)((1.0 - ratio) * n_time);
    int end_idx = n_time - 1;
    int n_range = end_idx - start_idx;

    VecXd num_meas_range(n_range);
    MatXd error_mat_range(4 * n_sat, n_range);

    num_meas_range = num_meas.segment(start_idx, n_range);
    error_mat_range = error_mat.block(0, start_idx, 4 * n_sat, n_range);

    // reshape error_mat_range from (4xnsat, n_range) to (4, n_sat*n_range)
    MatXd error_mat_range_reshaped(4, n_sat * n_range);
    for (int i = 0; i < n_sat; i++) {
      error_mat_range_reshaped.block(0, i * n_range, 4, n_range)
          = error_mat_range.block(i * 4, 0, 4, n_range);
    }

    // compute statistics ----------------------------------------------------
    // rms
    for (int i = 0; i < 4; i++) {
      rms(i) = RootMeanSquare(error_mat_range_reshaped.row(i));
      means(i) = error_mat_range_reshaped.row(i).mean();
      stds(i) = Std(error_mat_range_reshaped.row(i));
      p68(i) = Percentile(error_mat_range_reshaped.row(i), 0.68);
      p95(i) = Percentile(error_mat_range_reshaped.row(i), 0.95);
      p99(i) = Percentile(error_mat_range_reshaped.row(i), 0.99);
    }

    std::cout << " " << std::endl;
    std::cout << " " << std::endl;
    std::cout << "< Simulation Statistics (Last " << ratio * 100 << "%)>" << std::endl;
    std::cout << " " << std::endl;
    std::cout << "Statistics  | Position [m]  | Velocity [mm/s] | Clock Bias [ns] "
                 "| Clk Drift [ns/s] "
              << std::endl;
    std::cout << "---------------------------------------------------------------"
                 "-----------------------"
              << std::endl;

    std::cout.precision(5);
    std::cout << "RMS         | " << std::left << std::setw(16) << rms(0) << "  " << std::left
              << std::setw(16) << rms(1) << "   " << std::left << std::setw(16) << rms(2)
              << std::left << std::setw(16) << rms(3) << std::endl;

    std::cout << "Mean+-Std   | " << std::left << means(0) << "+-" << std::left << stds(0) << "  "
              << std::left << means(1) << "+-" << std::left << stds(1) << "   " << std::left
              << means(2) << "+-" << std::left << stds(2) << "   " << std::left << means(3) << "+-"
              << std::left << stds(3) << std::endl;

    std::cout << "68%         | " << std::left << std::setw(16) << p68(0) << "  " << std::left
              << std::setw(16) << p68(1) << "   " << std::left << std::setw(16) << p68(2)
              << std::left << std::setw(16) << p68(3) << std::endl;

    std::cout << "95%         | " << std::left << std::setw(16) << p95(0) << "  " << std::left
              << std::setw(16) << p95(1) << "   " << std::left << std::setw(16) << p95(2)
              << std::left << std::setw(16) << p95(3) << std::endl;

    std::cout << "99%         | " << std::left << std::setw(16) << p99(0) << "  " << std::left
              << std::setw(16) << p99(1) << "   " << std::left << std::setw(16) << p99(2)
              << std::left << std::setw(16) << p99(3) << std::endl;

    std::cout << "  " << std::endl;
  }

}  // namespace lupnt

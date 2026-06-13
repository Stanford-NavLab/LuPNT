/**
 * @file sim_r12.cpp
 * @brief Example C++ code to run GCPM v2.4 simulations for various Rz12 values.
 * This code computes equatorial and meridianal slices of electron density
 * for different solar activity levels (Rz12) and saves the results to CSV
 * files.
 * @author Keidai Iiyama
 * @date 2024-06-20
 */

#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "lupnt/environment/plasma/plasma.h"

using namespace std;
using namespace pecsim;

void run_equatorial_slice(const DateTime& datetime, double akp,
                          const std::filesystem::path& output_dir_cpp) {
  std::vector<double> outn(8);
  std::vector<double> outnf(4);
  // Output density array
  double den_cpp[201][201] = {0};
  double alatr = 0.0;  // Magnetic latitude in radians

  std::cout << "Computing equatorial slice..." << std::endl;

  for (int i = 0; i <= 200; ++i) {
    double x = (i - 100) / 10.0;

    if (i % 2 == 0) {
      std::cout << "Progress: " << i << " / " << 200 << " (" << (i * 1.0 / 200 * 100) << " %)"
                << std::endl;
    }

    for (int j = 0; j <= 200; ++j) {
      double y = (j - 100) / 10.0;
      double r = std::sqrt(x * x + y * y);
      double along = std::atan2(y, x);
      double amlt = along / M_PI * 12.0 + 12.0;
      if (amlt > 24.0) {
        amlt -= 24.0;
      }

      // Call the gcpm_interface function
      outn = gcpm_v24(datetime, r, amlt, alatr, akp);
      den_cpp[i][j] = outn[0];
    }
  }

  std::cout << "Writing equatorial slice data to CSV..." << std::endl;

  std::string csv_path_cpp = (output_dir_cpp / "test_equatorial.csv").string();
  std::ofstream csvFile(csv_path_cpp);

  for (int i = 0; i <= 200; ++i) {
    for (int j = 0; j <= 200; ++j) {
      csvFile << den_cpp[i][j] << ",";
    }
    csvFile << std::endl;
  }

  return;
}

void run_meridianal_slice(const DateTime& datetime, double akp,
                          const std::filesystem::path& output_dir_cpp) {
  double r, alatr, amlt;
  std::vector<double> outn(8);

  const int grid_size = 201;
  const float cell_resolution = 0.1f;

  // Variables
  std::vector<std::vector<double>> density_cpp(grid_size, std::vector<double>(grid_size, 0.0));

  double zmlt = 0.0;  // Zenithal magnetic local time
  double mlt_n = zmlt;
  double mlt_p = zmlt + 12.0;
  if (mlt_p >= 24.0) mlt_p -= 24.0;

  std::cout << "Computing meridianal slice..." << std::endl;

  for (int i = 0; i < grid_size; ++i) {
    double x = (i - 100) * cell_resolution;
    amlt = (x >= 0.0) ? mlt_p : mlt_n;

    for (int j = 0; j < grid_size; ++j) {
      double z = (j - 100) * cell_resolution;
      r = std::sqrt(x * x + z * z);

      if (r > 1.0) {
        alatr = atan2(z, std::abs(x));

        // C++ call
        outn = gcpm_v24(datetime, r, amlt, alatr, akp);
        density_cpp[i][j] = outn[0];  // Total electron density
      } else {
        density_cpp[i][j] = 0.0;
      }
    }

    if ((i + 1) % 20 == 0) {
      // Print progress every 10 iterations
      std::cout << "  Progress: " << i + 1 << " / " << grid_size << " (" << std::fixed
                << std::setprecision(2) << (static_cast<double>(i + 1) / grid_size * 100.0) << " %)"
                << std::endl;
    }
  }

  std::cout << "Writing meridianal slice data to CSV..." << std::endl;

  std::string filename = "gcpm_v24_meridian_" + std::to_string(static_cast<int>(mlt_n)) + "h_"
                         + std::to_string(static_cast<int>(mlt_p)) + "h_kp1p0.csv";

  std::string csv_path_cpp = (output_dir_cpp / filename).string();
  std::ofstream output_csv_cpp(csv_path_cpp);
  for (int i = 0; i < grid_size; ++i) {
    for (int j = 0; j < grid_size; ++j) {
      output_csv_cpp << density_cpp[i][j] << ",";
    }
    output_csv_cpp << std::endl;
  }
  output_csv_cpp.close();
}

void run_sim_rz12(double rz12) {
  // Parameters
  DateTime datetime;  // Date and time for the computation
  double akp;
  std::filesystem::path output_dir_cpp;

  // custom date and kp
  datetime.year = 2025;  // Year
  int month = 3;         // Month (1-12)
  int day = 1;           // Day of the month (1-31)
  datetime.hour = 12;    // Hour (0-23)
  datetime.min = 0;      // Minute (0-59)
  datetime.sec = 0.0;    // Second (0.0-59.

  akp = 3.0;  // KP index <====== CHANGE HERE

  int doy;
  mmdd_to_doy(datetime.year, month, day, doy);  // Convert to day of year
  datetime.doy = doy;                           // Day of the year (1-366)

  std::string sim_str = std::to_string(datetime.year) + "_" + std::to_string(month) + "_"
                        + std::to_string(day) + "_" + std::to_string(datetime.hour) + "_kp"
                        + std::to_string(static_cast<int>(akp * 1000));  //

  sim_str += "_rz" + std::to_string(static_cast<int>(rz12));

  output_dir_cpp = std::filesystem::path(get_base_path()) / "output" / sim_str / "cpp" / "csv";

  // setup IRI model ----------------------------------------------
  std::string version = "iri2007";
  set_iri_model(IRIModel::IRI_2007);

  IRI2007Option iri_option_2007 = IRI2007Option();
  iri_option_2007.R12 = rz12;
  set_iri2007_option(iri_option_2007);

  // Create the output directory ----------------------------------------------
  if (!std::filesystem::exists(output_dir_cpp)) {
    std::filesystem::create_directories(output_dir_cpp);
    std::cout << "Created output directory: " << output_dir_cpp << std::endl;
  }

  // Loop over a grid of x and y values (equatorial slice)
  run_equatorial_slice(datetime, akp, output_dir_cpp);

  // Meridianal slice
  run_meridianal_slice(datetime, akp, output_dir_cpp);

  return;
}

int main() {
  std::vector<double> rz12_values = {10.0, 50.0, 100.0, 150.0, 200.0};

  for (double rz12 : rz12_values) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Running simulation for Rz12 = " << rz12 << std::endl;
    run_sim_rz12(rz12);
    std::cout << "----------------------------------------" << std::endl;
  }

  std::cout << "All simulations completed." << std::endl;
}

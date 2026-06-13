#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "lupnt/environment/plasma/plasma.h"

using namespace pecsim;

int main() {
  // Fixed-date version of the equatorial-slice example. Use
  // ex_plasma_equatorial_slice_custom for a configurable date/Kp/R12 run.
  double den_cpp[201][201] = {0};
  double den_fortran[201][201] = {0};

  DateTime datetime;
  datetime.year = 2002;
  datetime.doy = 185;
  datetime.hour = 12;
  datetime.min = 0;
  datetime.sec = 0.0;

  const double akp = 0.7;
  const double rz12 = -1.0;
  const bool run_cpp = true;
  const bool run_fortran = true;
  const double alatr = 0.0;

  std::filesystem::path output_dir_cpp
      = std::filesystem::path(get_base_path()) / "output" / "test" / "cpp" / "csv";
  std::filesystem::path output_dir_fortran
      = std::filesystem::path(get_base_path()) / "output" / "test" / "fortran" / "csv";

  set_iri_model(IRIModel::IRI_2007);

  if (rz12 > 0) {
    IRI2007Option iri_option_2007;
    iri_option_2007.R12 = rz12;
    set_iri2007_option(iri_option_2007);
  }

  if (run_cpp && !std::filesystem::exists(output_dir_cpp)) {
    std::filesystem::create_directories(output_dir_cpp);
    std::cout << "Created output directory: " << output_dir_cpp << std::endl;
  }

  if (run_fortran && !std::filesystem::exists(output_dir_fortran)) {
    std::filesystem::create_directories(output_dir_fortran);
    std::cout << "Created output directory: " << output_dir_fortran << std::endl;
  }

  std::vector<double> outn(8);
  std::vector<double> outnf(4);

  for (int i = 0; i <= 200; ++i) {
    const double x = (i - 100) / 10.0;

    if (i % 2 == 0) {
      std::cout << "Progress: " << i << " / 200 (" << (i * 100.0 / 200) << " %)" << std::endl;
    }

    for (int j = 0; j <= 200; ++j) {
      const double y = (j - 100) / 10.0;
      const double r = std::sqrt(x * x + y * y);
      const double along = std::atan2(y, x);
      double amlt = along / M_PI * 12.0 + 12.0;
      if (amlt > 24.0) {
        amlt -= 24.0;
      }

      if (run_cpp) {
        outn = gcpm_v24(datetime, r, amlt, alatr, akp);
        den_cpp[i][j] = outn[0];
      }
      if (run_fortran) {
        outnf = gcpm_v24_fortran(datetime, r, amlt, alatr, akp);
        den_fortran[i][j] = outnf[0];

        if (outnf[0] < 0.0f || outnf[0] > 1.0e7f) {
          std::cerr << "Error: out of range value detected." << std::endl;
          std::cerr << "outn: " << outnf[0] << ", x: " << x << ", y: " << y << ", r: " << r
                    << ", amlt: " << amlt << ", alatr: " << alatr << std::endl;
          return 1;
        }
      }
    }
  }

  if (run_fortran) {
    std::ofstream csv_file_fortran((output_dir_fortran / "test_equatorial.csv").string());
    for (int i = 0; i <= 200; ++i) {
      for (int j = 0; j <= 200; ++j) {
        csv_file_fortran << den_fortran[i][j] << ",";
      }
      csv_file_fortran << std::endl;
    }
  }

  if (run_cpp) {
    std::ofstream csv_file_cpp((output_dir_cpp / "test_equatorial.csv").string());
    for (int i = 0; i <= 200; ++i) {
      for (int j = 0; j <= 200; ++j) {
        csv_file_cpp << den_cpp[i][j] << ",";
      }
      csv_file_cpp << std::endl;
    }
  }

  std::cout << "Computation completed successfully." << std::endl;
  return 0;
}

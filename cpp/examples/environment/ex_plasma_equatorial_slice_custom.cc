#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "example_config.h"
#include "lupnt/environment/plasma/plasma.h"

using namespace pecsim;

namespace {

  // Keep the example runnable both from the LuPNT repository root and from the
  // original gcpm_wrap-style plasma base directory.
  std::filesystem::path default_config_path() {
    const std::filesystem::path repo_path = std::filesystem::path("cpp") / "examples"
                                            / "environment" / "config"
                                            / "equatorial_slice_custom.cfg";
    if (std::filesystem::exists(repo_path)) {
      return repo_path;
    }

    return std::filesystem::path(get_base_path()) / "examples" / "config"
           / "equatorial_slice_custom.cfg";
  }

}  // namespace

int main(int argc, char** argv) {
  // This example samples the GCPM electron density in the magnetic equatorial
  // plane. The output CSV is useful for quick contour plots and regression
  // checks against the original Fortran model.
  double den_cpp[201][201] = {0};
  double den_fortran[201][201] = {0};

  const std::filesystem::path default_config = default_config_path();
  const std::filesystem::path config_path
      = argc > 1 ? std::filesystem::path(argv[1]) : default_config;

  const examples::ExampleConfig config = examples::load_example_config(config_path);
  std::cout << "Loaded config: " << config_path << std::endl;

  DateTime datetime;
  datetime.year = config.year;
  const int month = config.month;
  const int day = config.day;
  datetime.hour = config.hour;
  datetime.min = config.minute;
  datetime.sec = config.second;

  double akp = config.kp;
  double rz12 = config.rz12;
  bool run_cpp = config.run_cpp;
  bool run_fortran = config.run_fortran;
  const double alatr = 0.0;

  int doy;
  mmdd_to_doy(datetime.year, month, day, doy);
  datetime.doy = doy;

  std::string sim_str = std::to_string(datetime.year) + "_" + std::to_string(month) + "_"
                        + std::to_string(day) + "_" + std::to_string(datetime.hour) + "_kp"
                        + std::to_string(static_cast<int>(akp * 1000));

  if (rz12 > 0.0) {
    sim_str += "_rz" + std::to_string(static_cast<int>(rz12));
  } else {
    sim_str += "_rzauto";
  }

  // Plasma example outputs live below the configured plasma data directory, so
  // they do not pollute the source tree.
  std::filesystem::path output_dir_cpp
      = std::filesystem::path(get_base_path()) / "output" / sim_str / "cpp" / "csv";
  std::filesystem::path output_dir_fortran
      = std::filesystem::path(get_base_path()) / "output" / sim_str / "fortran" / "csv";

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

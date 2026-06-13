#include <cmath>
#include <fstream>
#include <iomanip>
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
                                            / "meridianal_slice_custom.cfg";
    if (std::filesystem::exists(repo_path)) {
      return repo_path;
    }

    return std::filesystem::path(get_base_path()) / "examples" / "config"
           / "meridianal_slice_custom.cfg";
  }

}  // namespace

int main(int argc, char** argv) {
  // This example samples a meridianal GCPM density slice. By default it runs
  // only one magnetic-local-time plane; set only_xz=false in the config for all
  // 12 MLT pairs.
  const int grid_size = 201;
  const float cell_resolution = 0.1f;

  std::vector<std::vector<double>> density_cpp(grid_size, std::vector<double>(grid_size, 0.0));
  std::vector<std::vector<double>> density_fortran(grid_size, std::vector<double>(grid_size, 0.0));

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
  double r12 = config.rz12;
  bool run_cpp = config.run_cpp;
  bool run_fortran = config.run_fortran;
  bool only_xz = config.only_xz;

  int doy;
  mmdd_to_doy(datetime.year, month, day, doy);
  datetime.doy = doy;

  std::string sim_str = std::to_string(datetime.year) + "_" + std::to_string(month) + "_"
                        + std::to_string(day) + "_" + std::to_string(datetime.hour) + "_kp"
                        + std::to_string(static_cast<int>(akp * 1000));

  if (r12 > 0.0) {
    sim_str += "_rz" + std::to_string(static_cast<int>(r12));
  } else {
    sim_str += "_rzauto";
  }

  std::filesystem::path output_dir_cpp
      = std::filesystem::path(get_base_path()) / "output" / sim_str / "cpp" / "csv";
  std::filesystem::path output_dir_fortran
      = std::filesystem::path(get_base_path()) / "output" / sim_str / "fortran" / "csv";

  std::cout << "Output directory (C++): " << output_dir_cpp << std::endl;
  std::cout << "Output directory (Fortran): " << output_dir_fortran << std::endl;

  set_iri_model(IRIModel::IRI_2007);
  if (r12 > 0.0) {
    IRI2007Option iri_option_2007;
    iri_option_2007.R12 = r12;
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

  for (double zmlt = 0.0; zmlt < 12.0; zmlt += 1.0) {
    double mlt_n = zmlt;
    double mlt_p = zmlt + 12.0;
    if (mlt_p >= 24.0) {
      mlt_p -= 24.0;
    }

    std::cout << "Computing for MLT range: " << mlt_n << "h - " << mlt_p << "h ..." << std::endl;

    for (int i = 0; i < grid_size; ++i) {
      const double x = (i - 100) * cell_resolution;
      const double amlt = (x >= 0.0) ? mlt_p : mlt_n;

      for (int j = 0; j < grid_size; ++j) {
        const double z = (j - 100) * cell_resolution;
        const double r = std::sqrt(x * x + z * z);

        if (r > 1.0) {
          const double alatr = atan2(z, std::abs(x));

          if (run_cpp) {
            outn = gcpm_v24(datetime, r, amlt, alatr, akp);
            density_cpp[i][j] = outn[0];
          }

          if (run_fortran) {
            outnf = gcpm_v24_fortran(datetime, r, amlt, alatr, akp);
            density_fortran[i][j] = outnf[0];
          }
        } else {
          density_cpp[i][j] = 0.0;
          density_fortran[i][j] = 0.0;
        }
      }

      if ((i + 1) % 20 == 0) {
        std::cout << "  Progress: " << i + 1 << " / " << grid_size << " (" << std::fixed
                  << std::setprecision(2) << (static_cast<double>(i + 1) / grid_size * 100.0)
                  << " %)" << std::endl;
      }
    }

    std::cout << "Finished computing for MLT range: " << mlt_n << "h - " << mlt_p << "h"
              << std::endl;
    std::cout << "   " << std::endl;

    const std::string filename = "gcpm_v24_meridian_" + std::to_string(static_cast<int>(mlt_n))
                                 + "h_" + std::to_string(static_cast<int>(mlt_p)) + "h_kp1p0.csv";

    if (run_fortran) {
      std::ofstream output_csv_fortran((output_dir_fortran / filename).string());
      for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
          output_csv_fortran << density_fortran[i][j] << ",";
        }
        output_csv_fortran << std::endl;
      }
    }

    if (run_cpp) {
      std::ofstream output_csv_cpp((output_dir_cpp / filename).string());
      for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
          output_csv_cpp << density_cpp[i][j] << ",";
        }
        output_csv_cpp << std::endl;
      }
    }

    if (only_xz) {
      break;
    }
  }

  return 0;
}

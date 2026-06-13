#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "lupnt/environment/plasma/plasma.h"

using namespace std;
using namespace pecsim;

int main() {
  // Parameters
  DateTime datetime;
  datetime.year = 2002;  // Calendar year
  datetime.doy = 185;    // Day of the year
  datetime.hour = 12;    // Hour of the day
  datetime.min = 0;      // Minute of the hour
  datetime.sec = 0;      // Second of the minute

  double akp = 0.7;  // Planetary Kp index
  double amlt = 23.74;
  double al = 1.107226364;
  double alat_max = acos(sqrt(1.014 / al));

  // Create the output directory
  std::string version = "iri2007";
  set_iri_model(IRIModel::IRI_2007);

  std::filesystem::path output_dir_cpp
      = std::filesystem::path(get_base_path()) / "output" / "cpp" / "csv";
  std::filesystem::path output_dir_fortran
      = std::filesystem::path(get_base_path()) / "output" / "fortran" / "csv";

  if (!std::filesystem::exists(output_dir_cpp)) {
    std::filesystem::create_directories(output_dir_cpp);
    std::cout << "Created output directory: " << output_dir_cpp << std::endl;
  }

  if (!std::filesystem::exists(output_dir_fortran)) {
    std::filesystem::create_directories(output_dir_fortran);
    std::cout << "Created output directory: " << output_dir_fortran << std::endl;
  }

  std::string csv_path_cpp = output_dir_cpp / "test_fieldaligned.csv";
  std::string csv_path_fortran = output_dir_fortran / "test_fieldaligned.csv";
  std::ofstream csv_file_cpp(csv_path_cpp);
  std::ofstream csv_file_fortran(csv_path_fortran);

  // Loop over magnetic latitude from (-alat_max, alat_max) in steps of 0.001 rad.
  double alatr, alatd, r;
  alatr = -alat_max;
  while (alatr <= alat_max) {
    alatd = alatr / 3.1415927 * 180.0;
    r = al * pow(cos(alatr), 2.0);
    std::vector<double> outn = gcpm_v24(datetime, r, amlt, alatr, akp);
    std::vector<double> outnf = gcpm_v24_fortran(datetime, r, amlt, alatr, akp);
    alatr += 0.001;
    // output to csv (alatd, r, outn[0])
    csv_file_cpp << alatd << "," << r << "," << outn[0] << std::endl;
    csv_file_fortran << alatd << "," << r << "," << outnf[0] << std::endl;
  }

  std::cout << "Finished computing for field-aligned slice" << std::endl;
  csv_file_cpp.close();
  csv_file_fortran.close();

  return 0;
}

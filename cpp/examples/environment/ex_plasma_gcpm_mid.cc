#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "lupnt/environment/plasma/plasma.h"

using namespace std;
using namespace pecsim;
using namespace std::chrono;

int main(int argc, char *argv[]) {
  // Parameters
  DateTime datetime;
  datetime.year = 2002;  // Calendar year
  datetime.doy = 185;    // Day of the year
  datetime.hour = 12;    // Hour of the day
  datetime.min = 0.;     // Minute of the hour
  datetime.sec = 0.;     // Second of the minute

  // Create the output directory
  set_iri_model(IRIModel::IRI_2007);

  std::vector<double> outn;
  std::vector<double> outn2;
  std::vector<double> outnf;

  double akp, amlt, al, alatr, alatd;
  double x, z, along;

  akp = 0.7;  // Planetary Kp index

  if (argc < 4) {
    cerr << "Usage: " << argv[0] << " <ni(int)> <latd0(double)> <r(double)>" << endl;
    return 1;  // Exit with error if no argument is provided
  }

  // error case: r = 2  lat_d0=64
  int ni = stoi(argv[1]);         // Number of latitude samples to run
  double lat_d0 = stod(argv[2]);  // Starting magnetic latitude in degrees
  double r = stod(argv[3]);       // Geocentric distance in Earth radii

  double dlat = PI / ni;            // Step size in radians for magnetic latitude
  double latr0 = lat_d0 * DEG2RAD;  // Starting magnetic latitude in

  std::cout << "Running test_gcpm_repeat with ni = " << ni << endl;

  std::cout << " " << endl;
  std::cout << "-------------------------------------------------------------" << endl;
  std::cout << "r  lat   C++_Density (2007)   Fortran_Density (2007)   "
               "C++_Density (2020)      TYPE"
            << std::endl;
  std::cout << "--------------------------------------------------------------" << std::endl;

  double dur_cpp = 0.0;       // Duration for C++ computation
  double dur_fortran = 0.0;   // Duration for Fortran computation
  double dur_cpp_2020 = 0.0;  // Duration for C++ computation with IRI-2020

  double mlt_p = 12.0;  // Magnetic local time in hours for positive longitude
  double mlt_n = 0.0;   // Magnetic local time in hours for negative longitude

  for (int i = 0; i < ni; i++) {
    // Case 2: Change the location
    alatr = i * dlat + latr0;                  // Magnetic latitude in radians
    alatd = alatr / M_PI * 180.0;              // Magnetic latitude in degrees
    x = r * std::cos(alatr);                   // Calculate the x coordinate
    z = r * std::sin(alatr);                   // Calculate the z coordinate
    amlt = (x >= 0.0) ? mlt_p + 12.0 : mlt_n;  // Magnetic local time in hours

    // C++ 2007
    set_iri_model(IRIModel::IRI_2007);  // Set the IRI model to IRI-2007
    auto start_time = high_resolution_clock::now();
    outn = gcpm_v24(datetime, r, amlt, alatr, akp);
    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end_time - start_time);
    dur_cpp += static_cast<double>(duration.count());

    // Fortran 2007
    auto start_time_f = high_resolution_clock::now();
    outnf = gcpm_v24_fortran(datetime, r, amlt, alatr, akp);
    auto end_time_f = high_resolution_clock::now();
    auto duration_f = duration_cast<microseconds>(end_time_f - start_time_f);
    dur_fortran += static_cast<double>(duration_f.count());

    // C++ 2020
    set_iri_model(IRIModel::IRI_2020);
    auto start_time_cpp_2020 = high_resolution_clock::now();
    outn2 = gcpm_v24(datetime, r, amlt, alatr, akp);
    auto end_time_cpp_2020 = high_resolution_clock::now();
    auto duration_cpp_2020 = duration_cast<microseconds>(end_time_cpp_2020 - start_time_cpp_2020);
    dur_cpp_2020 += static_cast<double>(duration_cpp_2020.count());

    int type = outn[7];
    std::string type_str;
    if (type == 0) {
      type_str = "low lat";
    } else if (type == 1) {
      type_str = "mid lat";
    } else if (type == 2) {
      type_str = "polar";
    }

    cout << r << "  " << int(alatd) << "     " << outn[0] << " (" << duration.count()
         << " us)       " << outnf[0] << " (" << duration_f.count() << " us)       " << outn2[0]
         << " (" << duration_cpp_2020.count() << " us)    " << type_str << endl;
  }

  std::cout << "-------------------------------------" << std::endl;
  std::cout << "Total C++ Duration (2007): " << dur_cpp / 1e6 << " s" << std::endl;
  std::cout << "Total C++ Duration (2020): " << dur_cpp_2020 / 1e6 << " s" << std::endl;
  std::cout << "Total Fortran Duration: " << dur_fortran / 1e6 << " s" << std::endl;
}

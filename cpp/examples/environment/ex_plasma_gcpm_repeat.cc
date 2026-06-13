#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "lupnt/environment/plasma/plasma.h"

using namespace std;
using namespace pecsim;

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

  std::vector<double> outn, outn1;
  std::vector<double> outnf, outnf1;

  double akp, amlt, al, r, alatr, alatd;
  double x, y, along;

  akp = 0.7;  // Planetary Kp index

  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " <run_as_first (0 or 1)>" << endl;
    return 1;  // Exit with error if no argument is provided
  }

  int run_as_first = stoi(argv[1]);  // Flag to run the first case

  // First call
  if (run_as_first == 0) {
    amlt = 23.74;
    al = 1.107226364;
    r = 1.014;
    alatr = acos(sqrt(r / al));
    alatd = alatr / 3.1415927 * 180.0;

    // Case 1
    outn1 = gcpm_v24(datetime, r, amlt, alatr, akp);
    outnf1 = gcpm_v24_fortran(datetime, r, amlt, alatr, akp);
  }

  // Case 2: Change the location
  x = 0.7;                            // Geocentric distance in Earth radii
  y = 2.5;                            // Magnetic latitude in degrees
  r = std::sqrt(x * x + y * y);       // Calculate the geocentric distance
  along = std::atan2(y, x);           // Calculate the magnetic local time
  amlt = along / M_PI * 12.0 + 12.0;  // Convert to hours
  if (amlt > 24.0) {
    amlt -= 24.0;  // Adjust for values greater than 24 hours
  }
  alatr = 0.0;                   // Magnetic latitude in radians
  alatd = alatr / M_PI * 180.0;  // Magnetic latitude in degrees

  outn = gcpm_v24(datetime, r, amlt, alatr, akp);
  outnf = gcpm_v24_fortran(datetime, r, amlt, alatr, akp);
  cout << "-----------------------------------------------------" << endl;
  cout << "First time: " << run_as_first << endl;
  cout << "  IRI Model: " << get_iri_model_str() << endl;
  cout << "  Year : " << datetime.year << " Day of Year: " << datetime.doy
       << " Hour: " << datetime.hour << endl;
  cout << "  Magnetic Latitude (rad) : " << alatr << endl;
  cout << "  Magnetic Latitude (deg) : " << alatd << endl;
  cout << "  Magnetic Local Time (hours): " << amlt << endl;
  cout << "  Geocentric Distance (Re): " << r << endl;
  cout << " " << endl;
  cout << "  f107 : " << outn[4] << endl;
  cout << "  rz12 : " << outn[5] << endl;
  cout << "  hmf2_km : " << outn[6] << endl;
  cout << " " << endl;
  cout << "  Density (e, H+, He+, O+) (C++)    : " << outn[0] << " " << outn[1] << " " << outn[2]
       << " " << outn[3] << endl;
  cout << "  Density (e, H+, He+, O+) (Fortran): " << outnf[0] << " " << outnf[1] << " " << outnf[2]
       << " " << outnf[3] << endl;
  cout << "-----------------------------------------------------" << endl;
  cout << " " << endl;
}

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
  datetime.min = 0.;     // Minute of the hour
  datetime.sec = 0.;     // Second of the minute

  double akp = 0.7;  // Planetary Kp index
  double amlt = 23.74;
  double al = 1.107226364;
  double r = 1.014;
  double alatr = acos(sqrt(r / al));
  double alatd = alatr / 3.1415927 * 180.0;

  double mjd = datetime_to_mjd(datetime);
  double tj2000 = mjd_to_tj2000(mjd);

  // Create the output directory
  set_iri_model(IRIModel::IRI_2007);

  std::vector<double> outn;
  std::vector<double> outnf;
  Vec3d iono_params;

  // Case 1
  outn = gcpm_v24(datetime, r, amlt, alatr, akp);
  outnf = gcpm_v24_fortran(datetime, r, amlt, alatr, akp);
  iono_params = get_iono_params(tj2000);
  cout << "-----------------------------------------------------" << endl;
  cout << "Case 1: " << endl;
  cout << "  IRI Model: " << get_iri_model_str() << endl;
  cout << "  Year : " << datetime.year << " Day of Year: " << datetime.doy
       << " Hour: " << datetime.hour << endl;
  cout << "  Magnetic Latitude (rad) : " << alatr << endl;
  cout << "  Magnetic Latitude (deg) : " << alatd << endl;
  cout << "  Geocentric Distance (Re): " << r << endl;
  cout << " " << endl;
  cout << "  f107 : " << outn[4] << " "
       << "(from iono_params:" << iono_params[0] << ")" << endl;
  cout << "  rz12 : " << outn[5] << " "
       << "(from iono_params:" << iono_params[1] << ")" << endl;
  cout << "  hmf2_km : " << outn[6] << " "
       << "(from iono_params:" << iono_params[2] << ")" << endl;
  cout << " " << endl;
  cout << "  Density (e, H+, He+, O+) (C++)    : " << outn[0] << " " << outn[1] << " " << outn[2]
       << " " << outn[3] << endl;
  cout << "  Density (e, H+, He+, O+) (Fortran): " << outnf[0] << " " << outnf[1] << " " << outnf[2]
       << " " << outnf[3] << endl;
  cout << "-----------------------------------------------------" << endl;
  cout << " " << endl;

  // Case 2: Change the location
  double x = 0.7;                     // Geocentric distance in Earth radii
  double y = 2.5;                     // Magnetic latitude in degrees
  r = std::sqrt(x * x + y * y);       // Calculate the geocentric distance
  double along = std::atan2(y, x);    // Calculate the magnetic local time
  amlt = along / M_PI * 12.0 + 12.0;  // Convert to hours
  if (amlt > 24.0) {
    amlt -= 24.0;  // Adjust for values greater than 24 hours
  }
  alatr = 0.0;                   // Magnetic latitude in radians
  alatd = alatr / M_PI * 180.0;  // Magnetic latitude in degrees

  outn = gcpm_v24(datetime, r, amlt, alatr, akp);
  outnf = gcpm_v24_fortran(datetime, r, amlt, alatr, akp);
  cout << "-----------------------------------------------------" << endl;
  cout << "Case 2: " << endl;
  cout << "  IRI Model: " << get_iri_model_str() << endl;
  cout << "  Year : " << datetime.year << " Day of Year: " << datetime.doy
       << " Hour: " << datetime.hour << endl;
  cout << "  Magnetic Latitude (rad) : " << alatr << endl;
  cout << "  Magnetic Latitude (deg) : " << alatd << endl;
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

  // Case 3: Change the date (2024)
  datetime.year = 2024;
  datetime.doy = 10;
  datetime.hour = 3;
  datetime.min = 0;  // Minute of the hour
  datetime.sec = 0;  // Second of the minute
  r = 2.2;

  akp = 6.5;  // Kp index

  outn = gcpm_v24(datetime, r, amlt, alatr, akp);
  outnf = gcpm_v24_fortran(datetime, r, amlt, alatr, akp);
  cout << "-----------------------------------------------------" << endl;
  cout << "Case 3: (Kp fixed)" << endl;
  cout << "  IRI Model: " << get_iri_model_str() << endl;
  cout << "  Year : " << datetime.year << " Day of Year: " << datetime.doy
       << " Hour: " << datetime.hour << endl;
  cout << "  Kp Index: " << akp << endl;
  cout << "  Magnetic Latitude (rad) : " << alatr << endl;
  cout << "  Magnetic Latitude (deg) : " << alatd << endl;
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

  // Case 4: Free kp index
  akp = get_kp_index(datetime);
  outn = gcpm_v24(datetime, r, amlt, alatr);           // no kp input
  outnf = gcpm_v24_fortran(datetime, r, amlt, alatr);  // no kp input

  cout << "-----------------------------------------------------" << endl;
  cout << "Case 4: (Kp index from data)" << endl;
  cout << "  IRI Model: " << get_iri_model_str() << endl;
  cout << "  Year : " << datetime.year << " Day of Year: " << datetime.doy
       << " Hour: " << datetime.hour << endl;
  cout << "  Kp Index: " << akp << endl;
  cout << "  Magnetic Latitude (rad) : " << alatr << endl;
  cout << "  Magnetic Latitude (deg) : " << alatd << endl;
  cout << "  Geocentric Distance (Re): " << r << endl;
  cout << " " << endl;
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

  // Case 5: Use IRI model 2020
  r = 1.014;
  set_iri_model(IRIModel::IRI_2020);
  outn = gcpm_v24(datetime, r, amlt, alatr);           // no kp input
  outnf = gcpm_v24_fortran(datetime, r, amlt, alatr);  // no kp input

  IRI2020Option op = IRI2020Option();
  op.compute_teti = false;     // Disable TETI computatio
  op.compute_ni = false;       // Disable NI computation
  op.output_messages = true;   // Disable output messages
  op.output_to_text = true;    // Output to text file
  op.plasma_model = 0;         // Use Ozhogin plasma model
  op.without_plasmapause = 1;  // Without plasmapause
  op.R12 = -1.0;               // Use historical or projected R12

  // set_iri2020_option(op);  // Set the IRI options

  std::vector<double> outni = iri_2020(datetime, r, amlt, alatr, akp, op);
  cout << "-----------------------------------------------------" << endl;
  cout << "Case 5: (Using IRI 2020)" << endl;
  cout << "  IRI Model: " << get_iri_model_str() << endl;
  cout << "  Year : " << datetime.year << " Day of Year: " << datetime.doy
       << " Hour: " << datetime.hour << endl;
  cout << "  Kp Index: " << akp << endl;
  cout << "  Magnetic Latitude (rad) : " << alatr << endl;
  cout << "  Magnetic Latitude (deg) : " << alatd << endl;
  cout << "  Geocentric Distance (Re): " << r << endl;
  cout << " " << endl;
  cout << " " << endl;
  cout << "  f107 : " << outn[4] << endl;
  cout << "  rz12 : " << outn[5] << endl;
  cout << "  hmf2_km : " << outn[6] << endl;
  cout << " " << endl;
  cout << "  Density (e, H+, He+, O+) (C++)    : " << outn[0] << " " << outn[1] << " " << outn[2]
       << " " << outn[3] << endl;
  cout << "  Density (e, H+, He+, O+) (Fortran): " << outnf[0] << " " << outnf[1] << " " << outnf[2]
       << " " << outnf[3] << endl;
  cout << "  Density (e, H+, He+, O+) (IRI2020): " << outni[0] << " " << outni[1] << " " << outni[2]
       << " " << outni[3] << endl;
  cout << "-----------------------------------------------------" << endl;
  cout << " " << endl;
}

/**
 * @file ex_plasma_iri.cc
 * @author Keidai Iiyama
 * @brief Compare the IRI-2007 and IRI-2020 interfaces used by GCPM.
 * @version 0.1
 * @date 2025-02-15
 *
 * @copyright Copyright (c) 2025
 *
 */

#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "lupnt/environment/plasma/plasma.h"

using namespace pecsim;
using namespace std;

int main() {
  // The legacy IRI interface uses a compact mmdd argument. Negative values
  // indicate day-of-year, matching the original Fortran convention.
  int jmag = 0;
  double alati = 45.;
  double along = -80.;
  int iyyyy = 2002;
  int iyyyy2 = 2024;
  int mmdd = -10;
  int mmdd2 = -100;
  double dhour = 12. + 25.;  // use universal time
  double dhour2 = 3. + 25.;  // use universal time
  double heibeg = 250.;
  double heiend = 250.;
  double heistp = 1.;

  /* *****************************************
   * IRI Model 2007
   *****************************************/

  // Test 2007 IRI Model with date 2002
  set_iri_model(IRIModel::IRI_2007);

  IRIParams params = iri_sub(jmag, alati, along, iyyyy, mmdd, dhour, heibeg, heiend, heistp);
  cout << " " << endl;
  cout << "--------------------------------" << endl;
  cout << "2007 IRI Model" << endl;
  cout << "Year: " << iyyyy << " Day of Year: " << -mmdd << " Hour: " << dhour - 25 << endl;
  cout << "Latitude: " << alati << " Longitude: " << along << endl;
  cout << "R12: " << params.rz12 << endl;
  cout << "F107: " << params.f107 << endl;
  cout << "F2 Peak Height: " << params.hmf2_km << endl;
  cout << "Total Electron Density: " << params.neiri << endl;
  cout << "--------------------------------" << endl;
  cout << " " << endl;

  // Test 2007 IRI Model with date 2024
  IRIParams params2 = iri_sub(jmag, alati, along, iyyyy2, mmdd2, dhour2, heibeg, heiend, heistp);
  cout << " " << endl;
  cout << "--------------------------------" << endl;
  cout << "2007 IRI Model" << endl;
  cout << "Year: " << iyyyy2 << " Day of Year: " << -mmdd2 << " Hour: " << dhour2 - 25 << endl;
  cout << "Latitude: " << alati << " Longitude: " << along << endl;
  cout << "R12: " << params2.rz12 << endl;
  cout << "F107: " << params2.f107 << endl;
  cout << "F2 Peak Height: " << params2.hmf2_km << endl;
  cout << "Total Electron Density: " << params2.neiri << endl;
  cout << "--------------------------------" << endl;
  cout << " " << endl;

  // Test 2007 IRI model with custom R12 value
  double r12 = 50.0;  // Example R12 value
  IRI2007Option iri_option_2007 = IRI2007Option();
  iri_option_2007.R12 = r12;
  set_iri2007_option(iri_option_2007);

  // Re-test 2007 IRI Model with custom R12 value
  IRIParams params3 = iri_sub(jmag, alati, along, iyyyy2, mmdd2, dhour2, heibeg, heiend, heistp);
  cout << " " << endl;
  cout << "--------------------------------" << endl;
  cout << "2007 IRI Model with Custom R12" << endl;
  cout << "Year: " << iyyyy2 << " Day of Year: " << -mmdd2 << " Hour: " << dhour2 - 25 << endl;
  cout << "Latitude: " << alati << " Longitude: " << along << endl;
  cout << "Custom R12: " << params3.rz12 << endl;
  cout << "F107: " << params3.f107 << endl;
  cout << "F2 Peak Height: " << params3.hmf2_km << endl;
  cout << "Total Electron Density: " << params3.neiri << endl;
  cout << "--------------------------------" << endl;
  cout << " " << endl;

  /* *****************************************
   * IRI Model 2020
   *****************************************/

  // Test 2020 IRI Model with date 2024
  set_iri_model(IRIModel::IRI_2020);

  IRIParams params4 = iri_sub(jmag, alati, along, iyyyy2, mmdd2, dhour2, heibeg, heiend, heistp);
  cout << " " << endl;
  cout << "--------------------------------" << endl;
  cout << "2020 IRI Model" << endl;
  cout << "Year: " << iyyyy2 << " Day of Year: " << -mmdd2 << " Hour: " << dhour2 - 25 << endl;
  cout << "Latitude: " << alati << " Longitude: " << along << endl;
  cout << "Rz12: " << params4.rz12 << endl;
  cout << "F107: " << params4.f107 << endl;
  cout << "F2 Peak Height: " << params4.hmf2_km << endl;
  cout << "Total Electron Density: " << params4.neiri << endl;
  cout << "--------------------------------" << endl;
  cout << " " << endl;

  // Test 2020 IRI Model with custom R12 value
  double r12_2020 = 10.0;  // Example R12 value
  IRI2020Option iri_option_2020 = IRI2020Option();
  iri_option_2020.R12 = r12_2020;
  set_iri2020_option(iri_option_2020);

  // Re-test 2020 IRI Model with custom R12 value
  IRIParams params5 = iri_sub(jmag, alati, along, iyyyy2, mmdd2, dhour2, heibeg, heiend, heistp);
  cout << " " << endl;
  cout << "--------------------------------" << endl;
  cout << "2020 IRI Model with Custom R12" << endl;
  cout << "Year: " << iyyyy2 << " Day of Year: " << -mmdd2 << " Hour: " << dhour2 - 25 << endl;
  cout << "Latitude: " << alati << " Longitude: " << along << endl;
  cout << "Custom R12: " << params5.rz12 << endl << "F107: " << params5.f107 << endl;
  cout << "F2 Peak Height: " << params5.hmf2_km << endl;
  cout << "Total Electron Density: " << params5.neiri << endl;
  cout << "--------------------------------" << endl;
  cout << " " << endl;
}

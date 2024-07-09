#include <lupnt/lupnt.h>

#include <filesystem>
#include <iostream>
using namespace lupnt;
using namespace std;
using namespace jpl_de;

int main() {
  EphemerisData data;
  ReadEphemerisFile(ASCII_KERNEL_DIR / "de440" / "header.440", data);

  // Output the parsed data for verification
  cout << "KSIZE: " << data.KSIZE << "\nNCOEFF: " << data.NCOEFF << endl;
  cout << "Start Info: " << data.start_info
       << "\nFinal Info: " << data.final_info << endl;
  cout << "Start Epoch: " << data.jd_start_tbd
       << "\nFinal Epoch: " << data.jd_end_tbd << "\nStep: " << data.step
       << endl;

  cout << "Constant Names:" << endl;
  int i = 0;
  for (const auto& name : data.constant_names) {
    cout << name << " ";
    if (++i > 10) break;
  }
  cout << endl;

  cout << "Constant Values:" << endl;
  i = 0;
  for (const auto& value : data.constant_values) {
    cout << value << " ";
    if (++i > 10) break;
  }
  cout << endl;

  cout << "Group 1050 Data:" << endl;
  cout << "Start Locations: ";
  for (const auto& loc : data.coeff_start_locations) {
    cout << loc << " ";
  }
  cout << "\nNumber of Coefficients: ";
  for (const auto& num : data.num_coefficients) {
    cout << num << " ";
  }
  cout << "\nNumber of Properties: ";
  for (const auto& sets : data.num_properties) {
    cout << sets << " ";
  }
  cout << endl;

  return 0;
}
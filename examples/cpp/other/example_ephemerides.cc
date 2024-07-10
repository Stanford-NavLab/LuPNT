#include <lupnt/lupnt.h>
#include <omp.h>

#include <filesystem>
#include <iostream>

using namespace std;
using namespace lupnt;
using namespace jpl_de;

int main() {
  EphemerisHeaderData data;
  ReadEphemerisHeaderFile(ASCII_KERNEL_DIR / "de440" / "header.440", data);

  // Output the parsed data for verification
  cout << "KSIZE: " << data.KSIZE << "\nNCOEFF: " << data.NCOEFF << endl;
  cout << "Start Info: " << data.start_info
       << "\nFinal Info: " << data.final_info << endl;
  cout << "Start Epoch: " << data.jd_tdb_start
       << "\nFinal Epoch: " << data.jd_tdb_end << "\nStep: " << data.step
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
  for (const auto& loc : data.coeff_offset) {
    cout << loc << " ";
  }
  cout << "\nNumber of Coefficients: ";
  for (const auto& num : data.n_coeffs) {
    cout << num << " ";
  }
  cout << "\nNumber of Properties: ";
  for (const auto& sets : data.n_subintervals) {
    cout << sets << " ";
  }
  cout << endl;

  ReadEphemerisCoefficientsFile(ASCII_KERNEL_DIR / "de440" / "ascp01950.440",
                                data);

  double jd_tdb = 2458832.6;
  Real t_tdb = JDtoTime(jd_tdb);
  Real t_tai = ConvertTime(t_tdb, "TDB", "TAI");
  Vec6 rv_mercury_ = GetBodyPosVel(t_tai, NaifId::SSB, NaifId::MERCURY);
  Vec3 r_mercury__ =
      GetBodyPosSpice(NaifId::MERCURY, t_tai, Frame::GCRF, NaifId::SSB, "NONE");
  Vec6 rv_mercury = GetBodyPosVel(t_tai, EphemID::MERCURY);
  auto fmt = Eigen::IOFormat(16, 0, ", ", ", ", "", "", "[", "]");
  cout << "Cheby (new)" << rv_mercury.head(3).transpose().format(fmt) << endl;
  cout << "Cheby (old)" << rv_mercury_.head(3).transpose().format(fmt) << endl;
  cout << "Spice      " << r_mercury__.head(3).transpose().format(fmt) << endl;
  cout << "Cheby (new)" << rv_mercury_.tail(3).transpose().format(fmt) << endl;
  cout << "Cheby (old)" << rv_mercury.tail(3).transpose().format(fmt) << endl;

  // Benchmark the ephemerides
  int n = 5'000'000;
  auto start = chrono::high_resolution_clock::now();
#pragma omp parallel for
  for (int i = 0; i < n; i++) {
    GetBodyPosVel(t_tai, EphemID::MERCURY);
  }
  auto end = chrono::high_resolution_clock::now();
  auto duration =
      chrono::duration_cast<chrono::milliseconds>(end - start).count();
  cout << "Total cheby (new, multiproc.) " << duration << " ms" << endl;

  start = chrono::high_resolution_clock::now();
  for (int i = 0; i < n; i++) {
    GetBodyPosVel(t_tai, EphemID::MERCURY);
  }
  end = chrono::high_resolution_clock::now();
  duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();
  cout << "Total cheby (new) " << duration << " ms" << endl;

  start = chrono::high_resolution_clock::now();
  for (int i = 0; i < n; i++) {
    GetBodyPosVel(t_tai, NaifId::SSB, NaifId::MERCURY);
  }
  end = chrono::high_resolution_clock::now();
  duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();
  cout << "Total cheby (old) " << duration << " ms" << endl;

  start = chrono::high_resolution_clock::now();
  for (int i = 0; i < n; i++) {
    GetBodyPosSpice(NaifId::MERCURY, t_tai, Frame::GCRF, NaifId::SSB, "NONE");
    // GetFrameConversionMat(t_tai, Frame::MOON_PA, Frame::ITRF);
  }
  end = chrono::high_resolution_clock::now();
  duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();
  cout << "Total spice " << duration << " ms" << endl;

  return 0;
}
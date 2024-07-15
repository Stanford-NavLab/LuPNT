#include <lupnt/lupnt.h>
#include <omp.h>

#include <ctime>

#include <filesystem>
#include <iostream>

using namespace std;
using namespace lupnt;
using namespace std::chrono;

int main() {
  double jd_tdb = 2458850.5;
  Real t_tdb = JDtoTime(jd_tdb);
  Real t_tai = ConvertTime(t_tdb, "TDB", "TAI");
  Vec6 rv_old = spice::GetBodyPosVel(t_tai, NaifId::SSB, NaifId::MERCURY);
  Vec6 rv_spi =
    spice::GetBodyPosVelSpice(t_tai, NaifId::SSB, NaifId::MERCURY);
  Vec6 rv_new = GetBodyPosVel(t_tai, NaifId::SSB, NaifId::MERCURY);
  auto fmt = Eigen::IOFormat(14, Eigen::DontAlignCols, ", ", ", ", "", "", "[", "]");
  cout << "Position and velocity" << endl;
  cout << "LuPNT       " << rv_new.transpose().format(fmt) << endl;
  cout << "LuPNT+Spice " << rv_old.transpose().format(fmt) << endl;
  cout << "Spice       " << rv_spi.transpose().format(fmt) << endl;

  // Benchmark 
  int n = 2'000'000;
  int w = 20;

  cout << endl << "Benchmark" << endl;
  auto start = high_resolution_clock::now();
  for (int i = 0; i < n; i++)
    GetBodyPosVel(t_tai, NaifId::SSB, NaifId::MERCURY);
  auto end = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(end - start).count();
  cout << setw(w) << left << "LuPNT" << duration << " ms" << endl;

  start = high_resolution_clock::now();
  int n_threads = omp_get_num_procs();
#pragma omp parallel for
  for (int i = 0; i < n; i++) {
    GetBodyPosVel(t_tai, NaifId::SSB, NaifId::MERCURY);
  }
  end = high_resolution_clock::now();
  duration = duration_cast<milliseconds>(end - start).count();
  cout << setw(w) << left << ("LuPNT (" + to_string(n_threads) + " threads)") << duration << " ms" << endl;

  start = high_resolution_clock::now();
  for (int i = 0; i < n; i++)
    GetBodyPosVel(t_tai, NaifId::SSB, NaifId::MERCURY);
  end = high_resolution_clock::now();
  duration = duration_cast<milliseconds>(end - start).count();
  cout << setw(w) << left << "LuPNT+Spice" << duration << " ms" << endl;

  start = high_resolution_clock::now();
  for (int i = 0; i < n; i++)
    spice::GetBodyPosSpice(t_tai, NaifId::SSB, NaifId::MERCURY);
  end = high_resolution_clock::now();
  duration = duration_cast<milliseconds>(end - start).count();
  cout << setw(w) << left << "Spice" << duration << " ms" << endl;

  return 0;
}
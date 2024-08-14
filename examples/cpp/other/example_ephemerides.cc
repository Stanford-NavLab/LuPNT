#include <lupnt/lupnt.h>
#include <omp.h>

#include <ctime>
#include <filesystem>
#include <iostream>

using namespace std;
using namespace lupnt;
using namespace std::chrono;

int main() {
  double jd_tdb = 2458832.6;
  Real t_tdb = JDtoTime(jd_tdb);
  Real t_tai = ConvertTime(t_tdb, "TDB", "TAI");
  NaifId ids[] = {NaifId::SUN,
                  NaifId::MERCURY_BARYCENTER,
                  NaifId::VENUS_BARYCENTER,
                  NaifId::EMB,
                  NaifId::EARTH,
                  NaifId::MOON,
                  NaifId::MARS_BARYCENTER,
                  NaifId::JUPITER_BARYCENTER,
                  NaifId::SATURN_BARYCENTER,
                  NaifId::URANUS_BARYCENTER,
                  NaifId::NEPTUNE_BARYCENTER};
  NaifId center = NaifId::SSB;
  NaifId target = NaifId::MOON;

  for (auto id : ids) {
    Vec6 rv_old = spice::GetBodyPosVel(t_tai, center, id);
    Vec6 rv_spi = spice::GetBodyPosVelSpice(t_tai, center, id);
    Vec6 rv_new = GetBodyPosVel(t_tai, center, id);
    auto fmt = Eigen::IOFormat(16, Eigen::AutoAlign, ", ", ", ", "", "", "[", "]");
    cout << toString(center) << " -> " << toString(id) << endl;
    cout << "LuPNT       " << rv_new.transpose().format(fmt) << endl;
    cout << "LuPNT+Spice " << rv_old.transpose().format(fmt) << endl;
    cout << "Spice       " << rv_spi.transpose().format(fmt) << endl;
    cout << endl;
  }

  // Benchmark
  int n = 1'000'000;
  int w = 20;

  cout << endl << "Benchmark" << endl;
  auto start = high_resolution_clock::now();
  for (int i = 0; i < n; i++) GetBodyPosVel(t_tai, center, target);
  auto end = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(end - start).count();
  cout << setw(w) << left << "LuPNT" << duration << " ms" << endl;

  start = high_resolution_clock::now();
  int n_threads = omp_get_num_procs();
#pragma omp parallel for
  for (int i = 0; i < n; i++) GetBodyPosVel(t_tai, center, target);
  end = high_resolution_clock::now();
  duration = duration_cast<milliseconds>(end - start).count();
  cout << setw(w) << left << ("LuPNT (" + to_string(n_threads) + " threads)") << duration << " ms"
       << endl;

  start = high_resolution_clock::now();
  for (int i = 0; i < n; i++) spice::GetBodyPosVel(t_tai, center, target);
  end = high_resolution_clock::now();
  duration = duration_cast<milliseconds>(end - start).count();
  cout << setw(w) << left << "LuPNT+Spice" << duration << " ms" << endl;

  start = high_resolution_clock::now();
  for (int i = 0; i < n; i++) spice::GetBodyPosVelSpice(t_tai, center, target);
  end = high_resolution_clock::now();
  duration = duration_cast<milliseconds>(end - start).count();
  cout << setw(w) << left << "Spice" << duration << " ms" << endl;
  return 0;
}

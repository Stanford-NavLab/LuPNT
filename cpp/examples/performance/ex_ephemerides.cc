// Benchmarks repeated ephemeris lookups for common solar-system bodies.
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
  Real t_tdb = JdToTime(jd_tdb);
  BodyId ids[] = {BodyId::SUN,
                  BodyId::MERCURY_BARYCENTER,
                  BodyId::VENUS_BARYCENTER,
                  BodyId::EMB,
                  BodyId::EARTH,
                  BodyId::MOON,
                  BodyId::MARS_BARYCENTER,
                  BodyId::JUPITER_BARYCENTER,
                  BodyId::SATURN_BARYCENTER,
                  BodyId::URANUS_BARYCENTER,
                  BodyId::NEPTUNE_BARYCENTER};
  BodyId center = BodyId::SSB;
  BodyId target = BodyId::MOON;
  Frame frame = Frame::ICRF;

  for (auto id : ids) {
    cout << int(center) << " -> " << int(id) << endl;
    Vec6 rv_old = spice::GetBodyPosVel(t_tdb, center, id);
    Vec6 rv_spi = spice::GetBodyPosVelSpice(t_tdb, center, id);
    Vec6 rv_new = GetBodyPosVel(t_tdb, center, id, frame);
    auto fmt = Eigen::IOFormat(16, Eigen::AutoAlign, ", ", ", ", "", "", "[", "]");
    cout << "LuPNT       " << rv_new.transpose().format(fmt) << endl;
    cout << "LuPNT+Spice " << rv_old.transpose().format(fmt) << endl;
    cout << "Spice       " << rv_spi.transpose().format(fmt) << endl;
    cout << endl;
  }

  // Benchmark
  int n = 20;
  int w = 20;

  Real dt = 10;
  VecX ts_tdb = Arange(t_tdb, t_tdb + 7 * SECS_DAY, dt);
  cout << endl << "Benchmark" << endl;

  // LuPNT ********************************************************************
  auto bar = Logger::GetProgressBar(n, "", "Ephemerides");
  auto start = high_resolution_clock::now();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < ts_tdb.size(); j++) {
      GetBodyPosVel(ts_tdb(j), center, target, frame);
    }
    bar->Update(i);
  }
  bar->Finish();
  auto end = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(end - start).count();
  cout << setw(w) << left << "LuPNT" << duration << " ms" << endl;

  bar = Logger::GetProgressBar(ts_tdb.size(), "", "Ephemerides");
  start = high_resolution_clock::now();
  for (int j = 0; j < ts_tdb.size(); j++) {
    for (int i = 0; i < n; i++) {
      GetBodyPosVel(ts_tdb(j), center, target, frame);
    }
    bar->Update(j);
  }
  bar->Finish();
  end = high_resolution_clock::now();
  duration = duration_cast<milliseconds>(end - start).count();
  cout << setw(w) << left << "LuPNT" << duration << " ms" << endl;

  // LuPNT ********************************************************************
  start = high_resolution_clock::now();
  int n_threads = omp_get_num_procs();
#pragma omp parallel for
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < ts_tdb.size(); j++) {
      GetBodyPosVel(ts_tdb(j), center, target, frame);
    }
  }
  bar->Finish();
  end = high_resolution_clock::now();
  duration = duration_cast<milliseconds>(end - start).count();
  cout << setw(w) << left << ("LuPNT (" + to_string(n_threads) + " threads)") << duration << " ms"
       << endl;

#pragma omp parallel for
  for (int j = 0; j < ts_tdb.size(); j++) {
    for (int i = 0; i < n; i++) {
      GetBodyPosVel(ts_tdb(j), center, target, frame);
    }
  }
  bar->Finish();
  end = high_resolution_clock::now();
  duration = duration_cast<milliseconds>(end - start).count();
  cout << setw(w) << left << ("LuPNT (" + to_string(n_threads) + " threads)") << duration << " ms"
       << endl;

  // LuPNT+Spice ***************************************************************
  bar = Logger::GetProgressBar(n, "", "Ephemerides");
  start = high_resolution_clock::now();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < ts_tdb.size(); j++) {
      spice::GetBodyPosVel(ts_tdb(j), center, target);
    }
    bar->Update(i);
  }
  bar->Finish();
  end = high_resolution_clock::now();
  duration = duration_cast<milliseconds>(end - start).count();
  cout << setw(w) << left << "LuPNT+Spice" << duration << " ms" << endl;

  bar = Logger::GetProgressBar(ts_tdb.size(), "", "Ephemerides");
  start = high_resolution_clock::now();
  for (int j = 0; j < ts_tdb.size(); j++) {
    for (int i = 0; i < n; i++) {
      spice::GetBodyPosVel(ts_tdb(j), center, target);
    }
    bar->Update(j);
  }
  bar->Finish();
  end = high_resolution_clock::now();
  duration = duration_cast<milliseconds>(end - start).count();
  cout << setw(w) << left << "LuPNT+Spice" << duration << " ms" << endl;

  // Spice ********************************************************************
  bar = Logger::GetProgressBar(n, "", "Ephemerides");
  start = high_resolution_clock::now();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < ts_tdb.size(); j++) {
      spice::GetBodyPosVelSpice(ts_tdb(j), center, target);
    }
    bar->Update(i);
  }
  bar->Finish();
  end = high_resolution_clock::now();
  duration = duration_cast<milliseconds>(end - start).count();
  cout << setw(w) << left << "Spice" << duration << " ms" << endl;

  bar = Logger::GetProgressBar(ts_tdb.size(), "", "Ephemerides");
  start = high_resolution_clock::now();
  for (int j = 0; j < ts_tdb.size(); j++) {
    for (int i = 0; i < n; i++) {
      spice::GetBodyPosVelSpice(ts_tdb(j), center, target);
    }
    bar->Update(j);
  }
  bar->Finish();
  end = high_resolution_clock::now();
  duration = duration_cast<milliseconds>(end - start).count();
  cout << setw(w) << left << "Spice" << duration << " ms" << endl;
  return 0;
}

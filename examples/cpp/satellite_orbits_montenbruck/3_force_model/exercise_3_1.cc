#include <lupnt/lupnt.h>
#include <omp.h>

#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>

using namespace lupnt;
using namespace std;
using namespace matplot;

// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
// applications. Berlin : New York: Springer, 2000.
// doi: 10.1007/978-3-642-58351-3.
int main() {
  int N_Step = 2'000'000;  // Recommended for 0.01 sec timer (Linux)
  int n_max = 50;
  int N_threads = 12;
  Vec3 r(6525.919, 1710.416, 2508.886);  // Position [km]

  bool normalized = true;
  std::string filename = "JGM3.cof";
  GravityField<Real> grav = ReadHarmonicGravityField<Real>(filename, n_max, n_max, normalized);

  cout << "Exercise 3-1: Gravity Field Computation " << endl << endl;
  cout << " Order   CPU Time [s]" << endl << endl;

  // Real
  vector<double> ns;
  vector<double> times;
  omp_set_num_threads(N_threads);
  for (int n = 0; n <= n_max; n += 10) {
    double start = omp_get_wtime();
#pragma omp parallel for
    for (int i = 0; i < N_Step; i++) {
      Vec3 a = AccelarationGravityField(r, grav.GM, grav.R, grav.CS, n, n);
    }
    double end = omp_get_wtime();
    cout << setw(4) << n << setprecision(2) << fixed << setw(13) << (end - start) << endl;
    ns.push_back(n);
    times.push_back(end - start);
  }

  // Double
  GravityField<double> grav_d
      = ReadHarmonicGravityField<double>(filename, n_max, n_max, normalized);
  Vec3d r_d = r.cast<double>();

  cout << endl << " Order   CPU Time [s]" << endl << endl;

  vector<double> times_d;
  for (int n = 0; n <= n_max; n += 10) {
    double start = omp_get_wtime();
#pragma omp parallel for
    for (int i = 0; i < N_Step; i++) {
      Vec3d a_d = AccelarationGravityField(r_d, grav_d.GM, grav_d.R, grav_d.CS, n, n);
    }
    double end = omp_get_wtime();
    cout << setw(4) << n << setprecision(2) << fixed << setw(13) << (end - start) << endl;
    times_d.push_back(end - start);
  }

  // Plot
  figure();
  plot(ns, times, "-o");
  hold(on);
  plot(ns, times_d, "-o");
  xlabel("Degree");
  ylabel("CPU Time [s]");
  title("Gravity Field Computation (" + to_string(N_Step) + " evaluations)");
  legend({"Real", "Double"});
  grid(on);
  show();

  return 0;
}

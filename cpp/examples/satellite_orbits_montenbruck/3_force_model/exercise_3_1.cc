// Montenbruck Exercise 3-1: benchmark gravity-field acceleration as harmonic
// degree/order increases.
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
int main(int argc, char** argv) {
  int N_Step = argc > 1 ? std::stoi(argv[1]) : 200'000;
  int n_max = 100;
  int n_step = 20;
  int N_threads = 12;
  Vec3 r(6525.919e3, 1710.416e3, 2508.886e3);  // Position [m]

  std::string filename = "grgm900c.cof";

  // Load through the Body factory so the gravity field matches the Body/Real
  // scalar type used by the dynamics code.
  Body moon = Body::Moon(n_max, n_max, filename);
  const GravityField<Real>& grav = moon.gravity_field;
  vector<double> ns;
  vector<double> times;
  omp_set_num_threads(N_threads);
  cout << "Exercise 3-1: Gravity Field Computation " << endl << endl;
  cout << " Order   CPU Time [s]" << endl << endl;
  for (int n = 0; n <= n_max; n += n_step) {
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
  //   GravityField<double> grav_d
  //       = ReadHarmonicGravityField<double>(filename, n_max, n_max,
  //       normalized);
  //   Vec3d r_d = r.cast<double>();
  //   vector<double> times_d;
  //   cout << endl << " Order   CPU Time [s]" << endl << endl;
  //   for (int n = 0; n <= n_max; n += n_step) {
  //     double start = omp_get_wtime();
  // #pragma omp parallel for
  //     for (int i = 0; i < N_Step; i++) {
  //       Vec3d a_d = AccelarationGravityField(r_d, grav_d.GM, grav_d.R,
  //       grav_d.CS, n, n);
  //     }
  //     double end = omp_get_wtime();
  //     cout << setw(4) << n << setprecision(2) << fixed << setw(13) << (end -
  //     start) << endl; times_d.push_back(end - start);
  //   }

  //   // Plot
  //   figure();
  //   plot(ns, times, "-o");
  //   hold(on);
  //   plot(ns, times_d, "-o");
  //   xlabel("Degree");
  //   ylabel("CPU Time [s]");
  //   title("Gravity Field Computation (" + to_string(N_Step) + "
  //   evaluations)"); matplot::legend({"Real", "Double"}); grid(on); show();

  return 0;
}

#include <lupnt/lupnt.h>
#include <omp.h>

#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>

using namespace lupnt;
using namespace std;

// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
// applications. Berlin : New York: Springer, 2000.
// doi: 10.1007/978-3-642-58351-3.
int main() {
  int N_Step = 2'000'000;  // Recommended for 0.01 sec timer (Linux)
  int n_max = 20;
  int N_threads = 20;
  Vec3 r(6525.919, 1710.416, 2508.886);  // Position [km]

  bool normalized = true;
  std::string filename = "JGM3.cof";
  auto fmt = Eigen::IOFormat(10, 0, ", ", "\n", "[", "]");
  GravityField grav = ReadHarmonicGravityField(filename, n_max, n_max, normalized);

  cout << "Exercise 3-1: Gravity Field Computation " << endl << endl;
  cout << " Order   CPU Time [s]" << endl << endl;

  omp_set_num_threads(N_threads);
  int n;
  for (n = 2; n <= n_max; n += 2) {
    double start = omp_get_wtime();
#pragma omp parallel for
    for (int i = 0; i < N_Step; i++) {
      Vec3 a = AccelarationGravityField(r, grav.GM, grav.R, grav.CS, n, n);
    }
    double end = omp_get_wtime();
    cout << setw(4) << n << setprecision(2) << fixed << setw(13) << (end - start) << endl;
  }
  return 0;
}

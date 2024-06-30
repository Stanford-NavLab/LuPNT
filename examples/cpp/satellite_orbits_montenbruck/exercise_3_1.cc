#include <lupnt/lupnt.h>

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
  const int N_Step = 10000;  // Recommended for 0.01 sec timer
  const int n_max = 5;
  const int m_max = 5;
  const bool normalized = true;
  const std::string filename = "JGM3.cof";
  const Vec3 r(6525.919e3, 1710.416e3, 2508.886e3);  // Position [m]

  auto fmt = Eigen::IOFormat(7, 0, ", ", "\n", "[", "]");
  GravityField grav =
      ReadHarmonicGravityField(filename, n_max, m_max, normalized);
  cout << grav.R << endl;
  cout << grav.GM << endl;
  cout << grav.CS.format(fmt) << endl;

  // int i, n;            // Loop counters
  // clock_t start, end;  // Processor time at start and end
  // double duration;
  // Vec3 a;

  // cout << "Exercise 3-1: Gravity Field Computation " << endl << endl;

  // cout << " Order   CPU Time [s]" << endl << endl;

  // // Outer loop [2,4,...,n_max]

  // for (int n = 2; n <= n_max; n += 2) {
  //   // Start timing
  //   start = clock();

  //   // Evaluate gravitational acceleration N_Step times
  //   for (int i = 0; i <= N_Step; i++)
  //     a = AccelHarmonic(r, Id(3), Grav.GM, Grav.R_ref, Grav.CS, n, n);

  //   // Stop CPU time measurement
  //   end = clock();

  //   duration = (double)(end - start) / (double)(CLOCKS_PER_SEC);

  //   cout << setw(4) << n << setprecision(2) << fixed << setw(13) << duration
  //        << endl;
  // };

  return 0;
}

#include <lupnt/lupnt.h>
#include <omp.h>

#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>

using namespace lupnt;
using namespace std;

Vec3d AccelHarmonic(const Vec3d& r, const Mat3d& E, double GM, double R_ref,
                    const MatXd& CS, int n_max, int m_max) {
  // Local variables

  int n, m;                // Loop counters
  double r_sqr, rho, Fac;  // Auxiliary quantities
  double x0, y0, z0;       // Normalized coordinates
  double ax, ay, az;       // Acceleration vector
  double C, S;             // Gravitational coefficients
  Vec3d r_bf;              // Body-fixed position
  Vec3d a_bf;              // Body-fixed acceleration

  MatXd V(n_max + 2, n_max + 2);  // Harmonic functions
  MatXd W(n_max + 2, n_max + 2);  // work array (0..n_max+1,0..n_max+1)

  // Body-fixed position

  r_bf = E * r;

  // Auxiliary quantities

  r_sqr = r_bf.squaredNorm();
  rho = R_ref * R_ref / r_sqr;

  x0 = R_ref * r_bf(0) / r_sqr;  // Normalized
  y0 = R_ref * r_bf(1) / r_sqr;  // coordinates
  z0 = R_ref * r_bf(2) / r_sqr;

  //
  // Evaluate harmonic functions
  //   V_nm = (R_ref/r)^(n+1) * P_nm(sin(phi)) * cos(m*lambda)
  // and
  //   W_nm = (R_ref/r)^(n+1) * P_nm(sin(phi)) * sin(m*lambda)
  // up to degree and order n_max+1
  //

  // Calculate zonal terms V(n,0); set W(n,0)=0.0

  V(0, 0) = R_ref / sqrt(r_sqr);
  W(0, 0) = 0.0;

  V(1, 0) = z0 * V(0, 0);
  W(1, 0) = 0.0;

  for (n = 2; n <= n_max + 1; n++) {
    V(n, 0) =
        ((2 * n - 1) * z0 * V(n - 1, 0) - (n - 1) * rho * V(n - 2, 0)) / n;
    W(n, 0) = 0.0;
  };

  // Calculate tesseral and sectorial terms

  for (m = 1; m <= m_max + 1; m++) {
    // Calculate V(m,m) .. V(n_max+1,m)

    V(m, m) = (2 * m - 1) * (x0 * V(m - 1, m - 1) - y0 * W(m - 1, m - 1));
    W(m, m) = (2 * m - 1) * (x0 * W(m - 1, m - 1) + y0 * V(m - 1, m - 1));

    if (m <= n_max) {
      V(m + 1, m) = (2 * m + 1) * z0 * V(m, m);
      W(m + 1, m) = (2 * m + 1) * z0 * W(m, m);
    };

    for (n = m + 2; n <= n_max + 1; n++) {
      V(n, m) =
          ((2 * n - 1) * z0 * V(n - 1, m) - (n + m - 1) * rho * V(n - 2, m)) /
          (n - m);
      W(n, m) =
          ((2 * n - 1) * z0 * W(n - 1, m) - (n + m - 1) * rho * W(n - 2, m)) /
          (n - m);
    };
  };

  //
  // Calculate accelerations ax,ay,az
  //

  ax = ay = az = 0.0;

  for (m = 0; m <= m_max; m++)
    for (n = m; n <= n_max; n++)
      if (m == 0) {
        C = CS(n, 0);  // = C_n,0
        ax -= C * V(n + 1, 1);
        ay -= C * W(n + 1, 1);
        az -= (n + 1) * C * V(n + 1, 0);
      } else {
        C = CS(n, m);      // = C_n,m
        S = CS(m - 1, n);  // = S_n,m
        Fac = 0.5 * (n - m + 1) * (n - m + 2);
        ax += +0.5 * (-C * V(n + 1, m + 1) - S * W(n + 1, m + 1)) +
              Fac * (+C * V(n + 1, m - 1) + S * W(n + 1, m - 1));
        ay += +0.5 * (-C * W(n + 1, m + 1) + S * V(n + 1, m + 1)) +
              Fac * (-C * W(n + 1, m - 1) + S * V(n + 1, m - 1));
        az += (n - m + 1) * (-C * V(n + 1, m) - S * W(n + 1, m));
      }

  // Body-fixed acceleration

  a_bf = (GM / (R_ref * R_ref)) * Vec3d(ax, ay, az);

  // Inertial acceleration

  return E.transpose() * a_bf;
}

// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
// applications. Berlin : New York: Springer, 2000.
// doi: 10.1007/978-3-642-58351-3.
int main() {
  // Constants

  const int N_Step = 500'000;  // Recommended for 0.01 sec timer (Linux)
  const int n_max = 20;
  const Vec3d r(6525.919, 1710.416, 2508.886);  // Position [m]

  // Variables
  const bool normalized = true;
  const std::string filename = "JGM3.cof";
  auto fmt = Eigen::IOFormat(10, 0, ", ", "\n", "[", "]");
  GravityField grav =
      ReadHarmonicGravityField(filename, n_max, n_max, normalized);

  int i, n;            // Loop counters
  clock_t start, end;  // Processor time at start and end
  double duration;
  Vec3d a;
  Mat3d E = Mat3d::Identity();

  // Header

  cout << "Exercise 3-1: Gravity Field Computation " << endl << endl;

  cout << " Order   CPU Time [s]" << endl << endl;

  // Outer loop [2,4,...,n_max]
  MatXd CS(n_max + 1, n_max + 1);
  for (int i = 0; i < n_max + 1; i++) {
    for (int j = 0; j < n_max + 1; j++) {
      CS(i, j) = grav.CS(i, j).val();
    }
  }
  double r_ref = grav.R.val();
  double gm = grav.GM.val();

  for (n = n_max; n <= n_max; n += 2) {
    // Start timing
    start = clock();

    // Evaluate gravitational acceleration N_Step times
    for (int i = 0; i < N_Step; i++) {
      a = AccelHarmonic(r, E, gm, r_ref, CS, n, n) / 1e3;
    }

    // Stop CPU time measurement
    end = clock();

    cout << fixed << setprecision(10);
    cout << a << endl;

    duration = (double)(end - start) / (double)(CLOCKS_PER_SEC);

    cout << setw(4) << n << setprecision(2) << fixed << setw(13) << duration
         << endl;
  };

  return 0;
}

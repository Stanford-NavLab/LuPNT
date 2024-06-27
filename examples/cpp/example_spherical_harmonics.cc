#include <lupnt/lupnt.h>
#include <omp.h>

#include <cmath>
#include <iostream>
#include <vector>

using namespace lupnt;

Vec3 SphericalHarmonicsGravity(const Vec3 &r, const MatX &C, const MatX &S,
                               int n_max, int m_max, Real R_body, Real GM,
                               bool normalized = true) {
  // Intermediate computations
  Real r_sqr = r.squaredNorm();
  Real r_sqr_inv = (r_sqr > EPS) ? 1.0 / r_sqr : 0.0;  // Safe division
  Real rho = R_body * R_body * r_sqr_inv;
  Real x0 = R_body * r[0] * r_sqr_inv;
  Real y0 = R_body * r[1] * r_sqr_inv;
  Real z0 = R_body * r[2] * r_sqr_inv;

  // Initialize Intermediary Matrices
  MatX V = MatX::Zero(n_max + 2, n_max + 2);
  MatX W = MatX::Zero(n_max + 2, n_max + 2);

  // Calculate zonal terms V(n, 0). Set W(n,0)=0.0
  V(0, 0) = R_body / sqrt(r_sqr);
  W(0, 0) = 0.0;

  V(1, 0) = z0 * V(0, 0);
  W(1, 0) = 0.0;

  for (int n = 2; n <= n_max + 1; ++n) {
    V(n, 0) =
        ((2 * n - 1) * z0 * V(n - 1, 0) - (n - 1) * rho * V(n - 2, 0)) / n;
    W(n, 0) = 0.0;
  }

  // Calculate tesseral and sectoral terms
  for (int m = 1; m <= m_max + 1; ++m) {
    // Calculate V(m,m) to V(n_max+1,m)
    V(m, m) = (2 * m - 1) * (x0 * V(m - 1, m - 1) - y0 * W(m - 1, m - 1));
    W(m, m) = (2 * m - 1) * (x0 * W(m - 1, m - 1) + y0 * V(m - 1, m - 1));

    if (m <= m_max) {
      V(m + 1, m) = (2 * m + 1) * z0 * V(m, m);
      W(m + 1, m) = (2 * m + 1) * z0 * W(m, m);
    }

    for (int n = m + 2; n <= n_max + 1; ++n) {
      V(n, m) =
          ((2 * n - 1) * z0 * V(n - 1, m) - (n + m - 1) * rho * V(n - 2, m)) /
          (n - m);
      W(n, m) =
          ((2 * n - 1) * z0 * W(n - 1, m) - (n + m - 1) * rho * W(n - 2, m)) /
          (n - m);
    }
  }

  // Calculate accelerations
  Real ax = 0.0;
  Real ay = 0.0;
  Real az = 0.0;

  for (int m = 0; m <= m_max; ++m) {
    for (int n = m; n <= n_max; ++n) {
      Real Cnm = 0.0;
      Real Snm = 0.0;
      if (m == 0) {
        // Denormalize Coefficients
        if (normalized) {
          Real N = sqrt(2 * n + 1);
          Cnm = N * C(n, 0);
        } else {
          Cnm = C(n, 0);
        }

        ax -= Cnm * V(n + 1, 1);
        ay -= Cnm * W(n + 1, 1);
        az -= (n + 1) * Cnm * V(n + 1, 0);

      } else {
        if (normalized) {
          Real N = sqrt((2 - ((m == 0) ? 1 : 0)) * (2 * n + 1) *
                        tgamma(n - m + 1) / tgamma(n + m + 1));
          Cnm = N * C(n, m);
          Snm = N * S(n, m);
        } else {
          Cnm = C(n, m);
          Snm = S(n, m);
        }

        Real fac = 0.5 * (n - m + 1) * (n - m + 2);
        ax += 0.5 * (-Cnm * V(n + 1, m + 1) - Snm * W(n + 1, m + 1)) +
              fac * (Cnm * V(n + 1, m - 1) + Snm * W(n + 1, m - 1));
        ay += 0.5 * (-Cnm * W(n + 1, m + 1) + Snm * V(n + 1, m + 1)) +
              fac * (-Cnm * W(n + 1, m - 1) + Snm * V(n + 1, m - 1));
        az += (n - m + 1) * (-Cnm * V(n + 1, m) - Snm * W(n + 1, m));
      }
    }
  }

  Vec3 a = (GM / (R_body * R_body)) * Vec3(ax, ay, az);

  return a;
}

int main() {
  // Example coefficients (truncated for simplicity)
  int n = 5;
  int m = 5;
  Body moon = Body::Moon(n, m);

  // Example positions
  std::vector<Vec3> positions = {{36e3, 0, 0}, {0, 2e3, 0}, {0, 0, -400}};

  // Compute accelerations for multiple positions
  std::vector<Vec3> accelerations(positions.size());

  omp_set_num_threads(1);
#pragma omp parallel for
  for (size_t i = 0; i < positions.size(); ++i) {
    accelerations[i] =
        SphericalHarmonicsGravity(positions[i], moon.Cnm, moon.Snm, moon.n_max,
                                  moon.m_max, moon.R, moon.mu, false);
  }

  for (size_t i = 0; i < positions.size(); ++i) {
    std::cout << "Position: " << positions[i].transpose()
              << " Acceleration: " << accelerations[i].transpose() << std::endl;
  }

  std::vector<Vec3> accelerations_2(positions.size());
  NBodyDynamics dyn_true;
  dyn_true.SetPrimaryBody(moon);
  for (size_t i = 0; i < positions.size(); ++i) {
    accelerations_2[i] = dyn_true.ComputeNBodyGravity(0, positions[i]);
  }

  for (size_t i = 0; i < positions.size(); ++i) {
    std::cout << "Position: " << positions[i].transpose()
              << " Acceleration: " << accelerations_2[i].transpose()
              << std::endl;
  }
  return 0;
}
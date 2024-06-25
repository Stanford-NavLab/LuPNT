#include <lupnt/core/constants.h>

#include <cmath>
#include <iomanip>
#include <iostream>

// Function to compute the eccentric anomaly
double EccAnom(double M, double e) {
  const double tolerance = 1e-10;
  double E = M;
  double delta = 1.0;

  while (std::abs(delta) > tolerance) {
    delta = (E - e * std::sin(E) - M) / (1 - e * std::cos(E));
    E -= delta;
  }
  return E;
}

int main() {
  Eigen::Vector2d MeanAnom(4, 50);  // [deg]

  for (int iCase = 0; iCase < 2; ++iCase) {
    double M = MeanAnom(iCase) * Rad;
    double e = 0.72;
    double E_ref = EccAnom(M, e);

    std::cout << "Exercise 2-2: Solution of Kepler's equation \n\n";
    std::cout << "  M = " << std::fixed << std::setprecision(11) << M << "\n";
    std::cout << "  e = " << std::fixed << std::setprecision(11) << e << "\n";
    std::cout << "  E = " << std::fixed << std::setprecision(11) << E_ref
              << "\n\n";

    // Newton's iteration
    std::cout << "  a) Newton's iteration \n\n";
    std::cout << "  i          E          Accuracy   sin/cos \n";

    double E = M;
    int i = 0;

    while (std::abs(E - E_ref) > 1e-10) {
      i++;
      E = E - (E - e * std::sin(E) - M) / (1 - e * std::cos(E));
      std::cout << std::setw(3) << i << " ";
      std::cout << std::fixed << std::setprecision(11) << E << " ";
      std::cout << std::scientific << std::setprecision(2)
                << std::abs(E - E_ref) << " ";
      std::cout << std::fixed << 2 * i << "\n";
    }

    std::cout << "\n";

    // Fixed point iteration
    std::cout << "  b) Fixed point iteration \n\n";
    std::cout << "  i          E          Accuracy   sin/cos \n";

    E = M;
    i = 0;

    while (std::abs(E - E_ref) > 1e-10) {
      i++;
      E = M + e * std::sin(E);
      std::cout << std::setw(3) << i << " ";
      std::cout << std::fixed << std::setprecision(11) << E << " ";
      std::cout << std::scientific << std::setprecision(2)
                << std::abs(E - E_ref) << " ";
      std::cout << std::fixed << i << "\n";
    }
  }

  return 0;
}
#include <lupnt/lupnt.h>

using namespace lupnt;
using std::cout, std::endl, std::setw, std::fixed, std::setprecision,
    std::scientific;

// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
// applications. Berlin : New York: Springer, 2000.
// doi: 10.1007/978-3-642-58351-3.

int main() {
  real e = 0.72;
  Vec2 mean_anomalies(4, 50);  // [deg]

  cout << "Exercise 2-2 (Kepler’s equation)";

  for (int i = 0; i < 2; ++i) {
    real M = mean_anomalies(i) * RAD_PER_DEG;
    real E_ref = MeanToEccentricAnomaly(M, e);

    cout << "\n\nExample " << i + 1 << "\n\n";
    cout << fixed << setprecision(11);
    cout << "  M = " << M << "\n";
    cout << "  e = " << e << "\n";
    cout << "  E = " << E_ref << "\n\n";

    // Newton's iteration
    cout << "a) Newton's iteration \n\n";
    cout << setw(5) << "i";
    cout << setw(15) << "E";
    cout << setw(10) << "ΔE";
    cout << setw(10) << "n_trig";
    cout << "\n";

    real E = M;
    int it = 0;
    while (abs(E - E_ref) > 1e-10) {
      it++;
      E = E - (E - e * sin(E) - M) / (1 - e * cos(E));
      cout << setw(5) << it;
      cout << setw(15) << fixed << setprecision(11) << E;
      cout << setw(10) << scientific << setprecision(2) << abs(E - E_ref);
      cout << setw(10) << fixed << 2 * it << " ";
      cout << "\n";
    }
    cout << "\n";

    // Fixed point iteration
    cout << "b) Fixed point iteration \n\n";
    cout << setw(5) << "i";
    cout << setw(15) << "E";
    cout << setw(10) << "ΔE";
    cout << setw(10) << "n_trig";
    cout << "\n";

    E = M;
    it = 0;
    while (abs(E - E_ref) > 1e-10) {
      it++;
      E = M + e * sin(E);
      cout << setw(5) << it;
      cout << setw(15) << fixed << setprecision(11) << E;
      cout << setw(10) << scientific << setprecision(2) << abs(E - E_ref);
      cout << setw(10) << fixed << it;
      cout << "\n";
    }
  }

  return 0;
}
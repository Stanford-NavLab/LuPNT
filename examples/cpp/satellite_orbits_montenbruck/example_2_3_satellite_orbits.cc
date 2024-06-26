#include <lupnt/lupnt.h>

using namespace lupnt;
using std::cout, std::endl, std::setw, std::fixed, std::setprecision,
    std::scientific;

// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
// applications. Berlin : New York: Springer, 2000.
// doi: 10.1007/978-3-642-58351-3.

int main() {
  // Position and velocity
  Vector3 r(10e3, 40e3, -5e3);  // [km]
  Vector3 v(-1.5, 1, -0.1);     // [km/s]

  Vector6 rv(6);
  rv << r, v;

  // Compute orbital elements
  Vector6 coe = CartesianToClassical(rv, GM_EARTH);

  // Output
  cout << "Exercise 2-3 (Osculating elements)\n\n";
  cout << fixed << setprecision(3);
  cout << "  r = " << r.transpose() << "  [km]\n";
  cout << "  v = " << v.transpose() << "  [km/s]\n";

  int w = 9;
  cout << "\nOrbital elements:\n\n";

  cout << fixed << setprecision(3);
  cout << "  Semimajor axis   a = ";
  cout << setw(w) << coe(0) << " [km]\n";
  cout << fixed << setprecision(7);
  cout << "  Eccentricity     e = ";
  cout << setw(w) << coe(1) << " [-]\n";
  cout << fixed << setprecision(3);
  cout << "  Inclination      i = ";
  cout << setw(w) << coe(2) * DEG_PER_RAD << " [deg]\n";
  cout << "  RA ascend. node  Ω = ";
  cout << setw(w) << coe(3) * DEG_PER_RAD << " [deg]\n";
  cout << "  Arg. of perigee  ω = ";
  cout << setw(w) << coe(4) * DEG_PER_RAD << " [deg]\n";
  cout << "  Mean anomaly     M = ";
  cout << setw(w) << coe(5) * DEG_PER_RAD << " [deg]\n";

  return 0;
}
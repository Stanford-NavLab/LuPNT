// Montenbruck Exercise 2-3: convert Cartesian position/velocity to osculating
// Keplerian elements.
#include <lupnt/lupnt.h>

using namespace lupnt;
using std::cout, std::endl, std::setw, std::fixed, std::setprecision, std::scientific;

// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
// applications. Berlin : New York: Springer, 2000.
// doi: 10.1007/978-3-642-58351-3.

int main() {
  const Vec3 r(10e6, 40e6, -5e6);        // [m]
  const Vec3 v(-1.5e3, 1.0e3, -0.10e3);  // [m/s]

  Vec6 y;
  y << r, v;
  Vec6 coe = CartToClassical(y, GM_EARTH);  // [m, -, rad, rad, rad, rad]

  cout << endl << "Exercise 2-3: Osculating elements" << endl << endl;

  cout << "State vector:" << endl << endl;
  cout << "  Position       " << fixed << setprecision(3);
  for (int i = 0; i < 3; i++) {
    cout << setw(12) << r(i) / 1000.0;
  };
  cout << "  [m]" << endl;
  cout << "  Velocity       " << setprecision(6);
  for (int i = 0; i < 3; i++) {
    cout << setw(12) << v(i) / 1000.0;
  };
  cout << "  [km/s]" << endl;
  cout << endl;

  cout << "Orbital elements:" << endl << endl;
  cout << setprecision(3);
  cout << "  Semimajor axis   " << setw(10) << coe(0) << " m" << endl;
  cout << setprecision(7);
  cout << "  Eccentricity     " << setw(10) << coe(1) << endl;
  cout << setprecision(3);
  cout << "  Inclination      " << setw(10) << coe(2) * DEG << " deg" << endl;
  cout << "  RA ascend. node  " << setw(10) << coe(3) * DEG << " deg" << endl;
  cout << "  Arg. of perigee  " << setw(10) << coe(4) * DEG << " deg" << endl;
  cout << "  Mean anomaly     " << setw(10) << coe(5) * DEG << " deg" << endl;
  cout << endl;
  return 0;
}

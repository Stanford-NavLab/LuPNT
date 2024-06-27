#include <lupnt/lupnt.h>

using namespace lupnt;
using namespace std;

// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
// applications. Berlin : New York: Springer, 2000.
// doi: 10.1007/978-3-642-58351-3.

int main() {
  const Real Ms[] = {4.0, 50.0};  // [deg]
  int i;
  Real M, E, e;

  cout << endl << "Exercise 2-2: Solution of Kepler's equation" << endl << endl;

  for (int i = 0; i <= 1; i++) {
    M = Ms[i] * RAD;
    e = 0.72;
    E = Mean2EccAnomaly(M, e);

    cout << fixed << setprecision(11);
    cout << "  M =" << setw(14) << M << endl;
    cout << "  e =" << setw(14) << e << endl;
    cout << "  E =" << setw(14) << E << endl << endl;
  }
  return 0;
}

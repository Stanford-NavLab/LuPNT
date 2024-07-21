#include <lupnt/lupnt.h>

#include <iostream>

using namespace lupnt;
namespace ad = autodiff;

Real f(Real x) { return 1 + x + x * x + 1 / x + log(x); }

int main() {
  Real x = 1.0;
  Real u = f(x);
  double dudx = ad::derivative(f, wrt(x), at(x));

  std::cout << "u = " << u << std::endl;
  std::cout << "du/dx = " << dudx << std::endl;
}

// Demonstrates LuPNT's autodiff wrappers for scalar and vector Jacobians.
// The first example differentiates a Vec3 -> Vec3 function; the second shows
// how to differentiate with respect to multiple vector arguments.
#include <lupnt/lupnt.h>
using namespace lupnt;

Vec3 func(Vec3 x, bool flag) { return flag ? x : x.array().square(); }

Vec1 func2(Vec3 a, Vec3 b) { return a.transpose() * b; }

int main() {
  Vec3 x(1, 2, 3);
  bool flag = false;

  Vec3 y;
  MatX dydx = jacobian(func, wrt(x), autodiff::at(x, flag), y);
  std::cout << dydx << std::endl;

  Vec3 a(1, 2, 3);
  Vec3 b(4, 5, 6);
  Vec1 rho;
  MatX H;
  jacobian(func2, wrt(a, b), autodiff::at(a, b), rho, H);
  std::cout << "rho: " << rho << std::endl;
  std::cout << "H: " << H << std::endl;
}

// Output:
// 2 0 0
// 0 4 0
// 0 0 6

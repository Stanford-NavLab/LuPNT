#include <lupnt/lupnt.h>
using namespace lupnt;

Vec3 func(Vec3 x, bool flag) { return flag ? x : x.array().square(); }

int main() {
  Vec3 x(1, 2, 3);
  bool flag = false;

  Vec3 y;
  MatX dydx = jacobian(func, wrt(x), at(x, flag), y);
  std::cout << dydx << std::endl;
}

// Output:
// 2 0 0
// 0 4 0
// 0 0 6

#include <lupnt/core/constants.h>

using namespace lupnt;

int main() {
  Vector3real a{1, 2, 3};
  Vector3real b(a);
  b(2) += 1;
  auto c = a.transpose() * b;
  auto d = b.norm();
  std::cout << a << std::endl;
  std::cout << b << std::endl;
  std::cout << c << std::endl;
  std::cout << d << std::endl;
}
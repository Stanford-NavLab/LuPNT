// Minimal Eigen matrix construction example using LuPNT's Vec/Mat aliases.
// The two assignments show the difference between column-wise and row-wise
// insertion with the Eigen comma initializer.
#include <lupnt/lupnt.h>

using namespace lupnt;

int main() {
  Vec3 a{1, 2, 3};
  Vec3 b{4, 5, 6};
  Vec3 c{7, 8, 9};
  Mat3 R = Mat3::Zero();
  R << a, b, c;
  std::cout << R << std::endl << std::endl;

  R << a.transpose(), b.transpose(), c.transpose();
  std::cout << R << std::endl << std::endl;
}

#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/orbit_state_utils.h>

using namespace lupnt;

int main() {
  Vector6 coe{9750.5,        0.7,          deg2rad(63.5),
              deg2rad(90.0), deg2rad(0.0), deg2rad(30.0)};
  auto cart = ClassicalToCartesian(coe, GM_MOON);
  auto coe2 = CartesianToClassical(cart, GM_MOON);
  auto cart2 = ClassicalToCartesian(coe2, GM_MOON);
  auto coe3 = CartesianToClassical(cart2, GM_MOON);
  auto cart3 = ClassicalToCartesian(coe3, GM_MOON);

  // Formatter
  auto fmt = Eigen::IOFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, " ",
                             " ", "", "", "", "");

  std::cout << coe.transpose().format(fmt) << std::endl;
  std::cout << coe2.transpose().format(fmt) << std::endl;
  std::cout << coe3.transpose().format(fmt) << std::endl;
  std::cout << cart.transpose().format(fmt) << std::endl;
  std::cout << cart2.transpose().format(fmt) << std::endl;
  std::cout << cart3.transpose().format(fmt) << std::endl;
}
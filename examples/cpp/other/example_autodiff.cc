#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/orbit_state.h>

using namespace lupnt;

int main() {
  Vec6 coe{9750.5, 0.7, DEG * 63.5, DEG * 90.0, 0.0, DEG * 30.0};
  auto cart = Classical2Cart(coe, GM_MOON);
  auto coe2 = Cart2Classical(cart, GM_MOON);
  auto cart2 = Classical2Cart(coe2, GM_MOON);
  auto coe3 = Cart2Classical(cart2, GM_MOON);
  auto cart3 = Classical2Cart(coe3, GM_MOON);

  // Formatter
  auto fmt
      = Eigen::IOFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ", "", "", "", "");

  std::cout << coe.transpose().format(fmt) << std::endl;
  std::cout << coe2.transpose().format(fmt) << std::endl;
  std::cout << coe3.transpose().format(fmt) << std::endl;
  std::cout << cart.transpose().format(fmt) << std::endl;
  std::cout << cart2.transpose().format(fmt) << std::endl;
  std::cout << cart3.transpose().format(fmt) << std::endl;
}

#include <lupnt/core/constants.h>
#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/orbit_state_utils.h>

using namespace lupnt;

int main() {
  double a = 6541.4;
  double e = 0.6;
  double i = 65.5 * RAD_PER_DEG;
  double Omega = 0.0 * RAD_PER_DEG;
  double w = 90.0 * RAD_PER_DEG;
  double M = 0.0 * RAD_PER_DEG;

  ad::Vector6real coe;
  coe << a, e, i, Omega, w, M;

  Eigen::Vector6d coe_d;
  coe_d << a, e, i, Omega, w, M;

  auto cart = CoeToCart(coe, MU_EARTH);
  auto cart_d = CoeToCart(coe_d, MU_EARTH);

  // Print type
  std::cout << "coe type = " << typeid(coe).name() << std::endl;
  std::cout << "coe_d type = " << typeid(coe_d).name() << std::endl;

  std::cout << "cart = " << cart.transpose() << std::endl;
  std::cout << "cart_d = " << cart_d.transpose() << std::endl;

  auto rv = CartesianOrbitState(cart);
  std::cout << "type = " << typeid(rv.r()).name() << std::endl;
}
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

  Vector6 coe;
  coe << a, e, i, Omega, w, M;

  Vector6 coe_d{a, e, i, Omega, w, M};

  auto cart = ClassicalToCartesian(coe, GM_EARTH);
  auto cart_d = ClassicalToCartesian(coe_d, GM_EARTH);

  // Print type
  std::cout << "coe type = " << typeid(coe).name() << std::endl;
  std::cout << "coe_d type = " << typeid(coe_d).name() << std::endl;

  std::cout << "cart = " << cart.transpose() << std::endl;
  std::cout << "cart_d = " << cart_d.transpose() << std::endl;

  auto rv = CartesianOrbitState(cart);
  std::cout << "type = " << typeid(rv.r()).name() << std::endl;

  Vector3 r_cart_ref{0, 0, 0};
  Vector3 r_aer{1e-3, 1e-3, 5e3};
  Vector3 r = AzimuthElevationRangeToCartesian(r_cart_ref, r_aer);
  std::cout << "r = " << r.transpose() << std::endl;
}
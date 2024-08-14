#include <lupnt/core/constants.h>
#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/coordinates.h>
#include <lupnt/physics/orbit_state.h>

using namespace lupnt;

int main() {
  double a = 6541.4;
  double e = 0.6;
  double i = 65.5 * RAD;
  double Omega = 0.0 * RAD;
  double w = 90.0 * RAD;
  double M = 0.0 * RAD;

  Vec6 coe;
  coe << a, e, i, Omega, w, M;

  Vec6 coe_d{a, e, i, Omega, w, M};

  auto cart = Classical2Cart(coe, GM_EARTH);
  auto cart_d = Classical2Cart(coe_d, GM_EARTH);

  // Print type
  std::cout << "coe type = " << typeid(coe).name() << std::endl;
  std::cout << "coe_d type = " << typeid(coe_d).name() << std::endl;

  std::cout << "cart = " << cart.transpose() << std::endl;
  std::cout << "cart_d = " << cart_d.transpose() << std::endl;

  auto rv = CartesianOrbitState(cart);
  std::cout << "type = " << typeid(rv.r()).name() << std::endl;

  Vec3 r_cart_ref{1, 1, 1};
  Vec3 r_aer{30 * RAD, 60 * RAD, 1};
  Vec3 r = AzElRange2Cart(r_aer, r_cart_ref);
  std::cout << "r = " << r.transpose() << std::endl;
}

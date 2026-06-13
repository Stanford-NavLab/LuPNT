// Basic coordinate-conversion examples: classical orbital elements to
// Cartesian state, and local azimuth/elevation/range to Cartesian position.
#include <lupnt/conversions/coordinate_conversions.h>
#include <lupnt/conversions/state_conversions.h>
#include <lupnt/core/constants.h>
#include <lupnt/numerics/math_utils.h>

#include "lupnt/states/state.h"

using namespace lupnt;

int main() {
  double a = 6541.4;
  double e = 0.6;
  double i = 65.5 * RAD;
  double Omega = 0.0 * RAD;
  double w = 90.0 * RAD;
  double M = 0.0 * RAD;

  Vec6 coe{a, e, i, Omega, w, M};
  State coe_s = ClassicalOE({a, e, i, Omega, w, M});

  auto rv = ClassicalToCart(coe, GM_EARTH);
  auto rv_s = ClassicalToCart(coe_s, GM_EARTH);

  // Print type
  std::cout << "cart   = " << rv.transpose() << std::endl;
  std::cout << "cart_s = " << rv_s.transpose() << std::endl;

  Vec3 r_ref{1, 1, 1};
  Vec3 r_aer{30 * RAD, 60 * RAD, 1};
  Vec3 r = AzElRangeToCart(r_aer, r_ref);
  std::cout << "r = " << r.transpose() << std::endl;
}

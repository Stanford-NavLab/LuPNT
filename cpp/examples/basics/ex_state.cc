// Shows how LuPNT State objects retain names, units, and frame metadata while
// still behaving like Eigen vectors.
#include <lupnt/lupnt.h>

using namespace lupnt;

int main() {
  Vec6 coe_vec = {6500, 0.001, 30 * RAD, 60 * RAD, 90 * RAD, 45 * RAD};
  std::cout << "coe_vec " << coe_vec.transpose().format(FMT_CLEAN) << std::endl;

  State coe = ClassicalOE({6500, 0.001, 30 * RAD, 60 * RAD, 90 * RAD, 45 * RAD}, Frame::MOON_CI);
  std::cout << "coe     " << coe.transpose().format(FMT_CLEAN) << std::endl;
  std::cout << "names   " << coe.GetNames() << std::endl;
  std::cout << "units   " << coe.GetUnits() << std::endl;
  std::cout << "frame   " << coe.GetFrame() << std::endl;

  State rv = ClassicalToCart(coe, GM_MOON);
  std::cout << "rv      " << rv.transpose().format(FMT_CLEAN) << std::endl;
  std::cout << "names   " << rv.GetNames() << std::endl;
  std::cout << "units   " << rv.GetUnits() << std::endl;
  std::cout << "frame   " << rv.GetFrame() << std::endl;
}

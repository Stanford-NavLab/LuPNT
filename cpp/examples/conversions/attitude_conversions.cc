// Small attitude-conversion sandbox for Euler angles, quaternions, and rotation
// matrices. It is intended for inspecting convention/sign behavior.
#include "lupnt/conversions/attitude_conversions.h"

using namespace lupnt;
int main() {
  Mat3 R;
  Vec3 e1(0.969742, 0.230052, -0.0817073);
  Vec3 e2(0.237688, -0.966084, 0.100927);
  Vec3 e3(-0.0557177, -0.117294, -0.991533);
  R << e1, e2, e3;

  std::cout << "Rot:\n" << R << std::endl;
  std::cout << "Quat:\n" << RotToQuat(R).transpose() << std::endl;
  std::cout << "Rot:\n" << QuatToRot(RotToQuat(R)) << std::endl;
  return 0;
}

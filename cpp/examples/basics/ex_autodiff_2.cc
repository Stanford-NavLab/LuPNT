// Demonstrates differentiating rotation matrices and quaternions with LuPNT's
// autodiff scalar type. It compares finite differences with autodiff and then
// maps the rotation derivative into angular velocity.
#include <lupnt/lupnt.h>

using namespace lupnt;

using Eigen::AngleAxisd;
using Quat = Eigen::Quaternion<Real>;

#include <unsupported/Eigen/EulerAngles>

using AngleAxis = Eigen::AngleAxis<Real>;
int main() {
  Real t = 0.1;
  Real eps = 1e-6;

  Vec3 ea(sin(t), cos(t), 3 * t * t);
  Mat3 R = RotX(ea(0)) * RotY(ea(1)) * RotZ(ea(2));
  Quat q = Quat(R).normalized();
  Quat aa = (AngleAxis(-ea(0), Vec3::UnitX()) * AngleAxis(-ea(1), Vec3::UnitY())
             * AngleAxis(-ea(2), Vec3::UnitZ()));

  std::cout << "Rotation matrix:\n\n";
  std::cout << R << "\n\n";
  std::cout << q.toRotationMatrix() << "\n\n";
  std::cout << aa.toRotationMatrix() << "\n\n";

  std::cout << "Quaternion:\n\n";
  std::cout << ScalarFirstToLast(q.coeffs()).transpose() << "\n\n";
  std::cout << q.w() << " " << q.vec().transpose() << "\n\n";

  std::cout << "Euler Angles:\n\n";
  std::cout << WrapToPi(ea).transpose() << "\n\n";
  std::cout << WrapToPi(R.transpose().eulerAngles(2, 1, 0).reverse()).transpose() << "\n\n";

  std::cout << "Differentiation:\n\n";
  auto GetR
      = [](Real t) { return (Vec9)(RotX(sin(t)) * RotY(cos(t)) * RotZ(3 * t * t)).reshaped(); };
  auto GetQuat = [](Real t) {
    return ScalarLastToFirst(Quat(RotX(sin(t)) * RotY(cos(t)) * RotZ(3 * t * t)).coeffs());
  };

  std::cout << (GetR(t + eps) - GetR(t - eps)).reshaped(3, 3) / (2 * eps) << "\n\n";
  std::cout << derivative(GetR, wrt(t), at(t)).reshaped(3, 3) << "\n\n";

  std::cout << "Angular velocity:\n\n";
  Mat3 R_dot = derivative(GetR, wrt(t), at(t)).reshaped(3, 3);
  Mat3 Omega = -R_dot * R.transpose();
  Vec3 omega(Omega(2, 1), Omega(0, 2), Omega(1, 0));
  std::cout << Omega << "\n\n";
  std::cout << omega.transpose() << "\n\n";

  std::cout << "Quaternion:\n\n";
  std::cout << (GetQuat(t + eps) - GetQuat(t - eps)).transpose() / (2 * eps) << "\n\n";
  std::cout << derivative(GetQuat, wrt(t), at(t)).transpose() << "\n\n";
  std::cout << 0.5 * (Quat(0, omega(0), omega(1), omega(2)) * Quat(GetQuat(t))).coeffs() << "\n\n";
}

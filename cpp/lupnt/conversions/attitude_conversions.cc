/**
 * @file attitude_conversions.cc
 * @author Stanford NAV Lab
 * @brief
 * @version 0.1
 * @date 2024-08-11
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "lupnt/conversions/attitude_conversions.h"

namespace lupnt {

  Vec3 RotToAngularVelocity(const Mat3& rot_b2w, const Mat3& rot_dot_b2w) {
    Mat3 Omega = rot_dot_b2w * rot_b2w.transpose();
    Omega = 0.5 * (Omega - Omega.transpose());
    return Vec3(Omega(2, 1), Omega(0, 2), Omega(1, 0));
  }

  State ScalarFirstToLast(const State& quat) {
    return NormalizeQuat(State(Vec4(quat(3), quat(0), quat(1), quat(2))));
  }
  State ScalarLastToFirst(const State& quat) {
    return NormalizeQuat(State(Vec4(quat(1), quat(2), quat(3), quat(0))));
  }

  State NormalizeQuat(const State& quat) {
    State q = quat.normalized();
    if (q(0) < 0) q = -q;
    return q;  // [qw, qx, qy, qz] = [scalar, vector]
  }

  Mat3 QuatToRot(const State& quat) {
    Eigen::Quaternion<Real> q(quat(0), quat(1), quat(2), quat(3));  // [w, x, y, z]
    return q.toRotationMatrix();
  }

  State RotToQuat(const Mat3& rot) {
    Eigen::Quaternion<Real> q(rot);
    return NormalizeQuat(State(Vec4(q.w(), q.x(), q.y(), q.z())));
  }

  Mat3 RollPitchYawToRot(const Vec3& rpy) {
    Eigen::AngleAxis<Real> roll(rpy(0), Vec3::UnitX());
    Eigen::AngleAxis<Real> pitch(rpy(1), Vec3::UnitY());
    Eigen::AngleAxis<Real> yaw(rpy(2), Vec3::UnitZ());
    return (yaw * pitch * roll).toRotationMatrix();
  }

  Vec3 RotToRollPitchYaw(const Mat3& rot) {
    auto ypr = rot.eulerAngles(2, 1, 0);
    return Vec3(ypr[2], ypr[1], ypr[0]);
  }

  // Vector definitions
  VEC_IMP_VECTOR(ScalarFirstToLast, 4);
  VEC_IMP_VECTOR(ScalarLastToFirst, 4);
  VEC_IMP_VECTOR(NormalizeQuat, 4);

}  // namespace lupnt

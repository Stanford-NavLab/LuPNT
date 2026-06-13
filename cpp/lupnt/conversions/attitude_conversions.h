/**
 * @file attitude_conversions.h
 * @author Stanford NAV Lab
 * @brief
 * @version 0.1
 * @date 2024-08-11
 *
 * @copyright Copyright (c) 2024
 *
 */

#pragma once

#include "lupnt/core/constants.h"
#include "lupnt/numerics/vector_macros.h"
#include "lupnt/states/state.h"

namespace lupnt {

  Vec3 RotToAngularVelocity(const Mat3& rot_b2w, const Mat3& rot_dot_b2w);

  State ScalarFirstToLast(const State& quat);
  State ScalarLastToFirst(const State& quat);
  State NormalizeQuat(const State& quat);

  Mat3 QuatToRot(const State& quat);
  State RotToQuat(const Mat3& rot);

  Mat3 RollPitchYawToRot(const Vec3& rpy);
  Vec3 RotToRollPitchYaw(const Mat3& rot);

  // Vector definitions
  VEC_DEF_VECTOR(ScalarFirstToLast, 4);
  VEC_DEF_VECTOR(ScalarLastToFirst, 4);
  VEC_DEF_VECTOR(NormalizeQuat, 4);

}  // namespace lupnt

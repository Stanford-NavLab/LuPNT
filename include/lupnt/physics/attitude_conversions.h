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

namespace lupnt {

  /**
   * @brief Convert a quaternion to a direction cosine matrix
   *
   * @param quat quaternion [q0, q1, q2, q3] = [scalar, vector]
   * @return Mat3 direction cosine matrix
   */
  Mat3 QuatCoeffToDCM(Vec4 quat);

  /**
   * @brief Convert a direction cosine matrix to a quaternion
   *
   * @param dcm direction cosine matrix (3x3)
   * @return Vec4 quaternion [q0, q1, q2, q3] = [scalar, vector]
   */
  Vec4 DCMToQuatCoeff(Mat3 dcm);

  /**
   * @brief Convert a quaternion to euler angles
   * Reference:
   * https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
   *
   * @param quat quaternion [q0, q1, q2, q3] = [scalar, vector]
   * @return Vec3 euler angles [phi, theta, psi] = [roll, pitch, yaw] (radians)
   * sequence (XYZ)
   *
   */
  Vec3 QuatCoeffToEulerAngles(Vec4 quat);

  /**
   * @brief Convert euler angles to a quaternion
   *
   * @param euler Vec3 euler angles [phi, theta, psi] = [roll, pitch, yaw]
   * (radians) sequence (XYZ)
   * @return Vec4 quaternion [q0, q1, q2, q3] = [scalar, vector]
   */
  Vec4 EulerAnglesToQuatCoeff(Vec3 euler);

  /**
   * @brief Convert a direction cosine matrix to euler angles
   *
   * @param dcm direction cosine matrix (3x3)
   * @return Vec3 euler angles [phi, theta, psi] (radians) sequence (1, 2, 3)
   */
  Vec3 DCMToEulerAngles(Mat3 dcm);

  /**
   * @brief Convert euler angles to a direction cosine matrix
   *
   * @param euler euler angles [phi, theta, psi] (radians) sequence (1, 2, 3)
   * @return Mat3 direction cosine matrix
   */
  Mat3 EulerAnglesToDCM(Vec3 euler);

}  // namespace lupnt

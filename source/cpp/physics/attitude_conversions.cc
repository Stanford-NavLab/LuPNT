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

#include "lupnt/physics/attitude_conversions.h"

namespace lupnt {

  /**
   * @brief Convert a quaternion to a direction cosine matrix
   *  Reference:
   *   https://stevendumble.com/attitude-representations-understanding-direct-cosine-matrices-euler-angles-and-quaternions/
   *
   * @param quat quaternion [q0, q1, q2, q3] = [scalar, vector]
   * @return Mat3 direction cosine matrix
   *
   */
  Mat3 QuatCoeffToDCM(Vec4 quat) {
    Real q0 = quat(0);
    Real q1 = quat(1);
    Real q2 = quat(2);
    Real q3 = quat(3);
    Mat3 dcm;
    dcm(0, 0) = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
    dcm(0, 1) = 2 * (q1 * q2 + q0 * q3);
    dcm(0, 2) = 2 * (q1 * q3 - q0 * q2);
    dcm(1, 0) = 2 * (q1 * q2 - q0 * q3);
    dcm(1, 1) = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
    dcm(1, 2) = 2 * (q2 * q3 + q0 * q1);
    dcm(2, 0) = 2 * (q1 * q3 + q0 * q2);
    dcm(2, 1) = 2 * (q2 * q3 - q0 * q1);
    dcm(2, 2) = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
    return dcm;
  }

  Vec4 DCMToQuatCoeff(Mat3 dcm) {
    Real q0 = 0.5 * sqrt(1 + dcm(0, 0) + dcm(1, 1) + dcm(2, 2));
    Real q1 = 0.25 * (dcm(1, 2) - dcm(2, 1)) / q0;
    Real q2 = 0.25 * (dcm(2, 0) - dcm(0, 2)) / q0;
    Real q3 = 0.25 * (dcm(0, 1) - dcm(1, 0)) / q0;
    return Vec4(q0, q1, q2, q3);
  }

  Vec3 QuatCoeffToEulerAngles(Vec4 quat) {
    Real q0 = quat(0);
    Real q1 = quat(1);
    Real q2 = quat(2);
    Real q3 = quat(3);
    Real phi = atan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1 * q1 + q2 * q2));
    Real theta = asin(2 * (q0 * q2 - q3 * q1));
    Real psi = atan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2 * q2 + q3 * q3));

    return Vec3(phi, theta, psi);
  }

  /**
   * @brief Convert euler angles to a quaternion
   *
   * @param euler Vec3 euler angles [phi, theta, psi] = [roll, pitch, yaw]
   * (radians) sequence (XYZ)
   * @return Vec4 quaternion [q0, q1, q2, q3] = [scalar, vector]
   */
  Vec4 EulerAnglesToQuatCoeff(Vec3 euler) {
    Real phi = euler(0);
    Real theta = euler(1);
    Real psi = euler(2);
    Real q0 = cos(phi / 2) * cos(theta / 2) * cos(psi / 2)
              + sin(phi / 2) * sin(theta / 2) * sin(psi / 2);
    Real q1 = sin(phi / 2) * cos(theta / 2) * cos(psi / 2)
              - cos(phi / 2) * sin(theta / 2) * sin(psi / 2);
    Real q2 = cos(phi / 2) * sin(theta / 2) * cos(psi / 2)
              + sin(phi / 2) * cos(theta / 2) * sin(psi / 2);
    Real q3 = cos(phi / 2) * cos(theta / 2) * sin(psi / 2)
              - sin(phi / 2) * sin(theta / 2) * cos(psi / 2);
    return Vec4(q0, q1, q2, q3);
  }

  Vec3 DCMToEulerAngles(Mat3 dcm) {
    Vec4 quat = DCMToQuatCoeff(dcm);
    return QuatCoeffToEulerAngles(quat);
  }

  Mat3 EulerAnglesToDCM(Vec3 euler) {
    Vec4 quat = EulerAnglesToQuatCoeff(euler);
    return QuatCoeffToDCM(quat);
  }

}  // namespace lupnt

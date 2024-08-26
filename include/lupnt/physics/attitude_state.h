#pragma once

#include "lupnt/core/constants.h"
#include "lupnt/physics/attitude_conversions.h"
#include "lupnt/physics/frame_converter.h"
#include "lupnt/physics/state.h"

namespace lupnt {

  class AttitudeState : public IState {
  private:
    VecX x_;  // Either quaternion or quaternion + angular velocity

    bool is_attached_to_frame_ = false;  // true if the attitude is wrt to a orbit
                                         // frame

    bool with_angular_velocity_ = false;  // With Angular Velocity

    Frame frame_ = Frame::NONE;
    AttitudeState *attached_to_ = nullptr;

  public:
    AttitudeState(const VecX &x) : x_(x) {
      if (x.size() == 7) {
        with_angular_velocity_ = true;
      } else if (x.size() != 4) {
        with_angular_velocity_ = false;
        throw std::invalid_argument("Invalid size for attitude state");
      }
    }

    // Only set attitude
    AttitudeState(const Quat &x) {
      x_.resize(4);
      x_ = x.coeffs();
      with_angular_velocity_ = false;
    }

    AttitudeState(const Mat3 &dcm) {
      x_.resize(4);
      x_.head(4) = DCMToQuatCoeff(dcm);
      with_angular_velocity_ = false;
    }

    AttitudeState(const Vec3 &euler_angles) {
      x_.resize(4);
      x_.head(4) = EulerAnglesToQuatCoeff(euler_angles);
      with_angular_velocity_ = false;
    }

    AttitudeState(const Vec3 &e_x, const Vec3 &e_y, const Vec3 &e_z) {
      x_.resize(4);
      SetQuatCoeff(e_x, e_y, e_z);
      with_angular_velocity_ = false;
    }

    AttitudeState(const Vec3 &e_x, const Vec3 &e_y, const Vec3 &e_z, Frame frame) {
      x_.resize(4);
      SetQuatCoeff(e_x, e_y, e_z);
      with_angular_velocity_ = false;
      frame_ = frame;
      is_attached_to_frame_ = true;
    }

    VecX GetVec() const { return x_; }
    inline int GetSize() const { return x_.size(); }
    inline Real GetValue(int i) const { return x_(i); }
    inline void SetValue(int i, Real val) { x_(i) = val; }
    inline void SetVecX(const VecX &x) { x_ = x; }

    // Angular Velocity
    inline bool HasAngularVelocity() const { return with_angular_velocity_; }

    // Frames
    inline Frame GetFrame() const { return frame_; }
    inline void AttachToFrame(Frame frame) {
      frame_ = frame;
      is_attached_to_frame_ = true;
    }
    inline bool IsWrtToFrame() const { return is_attached_to_frame_; }

    // Attitude Setters
    inline void SetQuatCoeff(const Vec4 &q) {
      if (q.size() != 4) {
        throw std::invalid_argument("Invalid size for quaternion");
      }
      x_.head(4) = q;
    }

    inline void SetQuatCoeff(const Vec3 &e_x, const Vec3 &e_y, const Vec3 &e_z) {
      Mat3 DCM;
      DCM.row(0) = e_x;
      DCM.row(1) = e_y;
      DCM.row(2) = e_z;
      x_.head(4) = DCMToQuatCoeff(DCM);
    }

    inline void SetQuatCoeff(const Mat3 &dcm) { x_.head(4) = DCMToQuatCoeff(dcm); }

    inline void SetEulerAngles(const Vec3 &euler_angles) {
      x_.head(4) = EulerAnglesToQuatCoeff(euler_angles);
    }

    inline void SetAngularVelocity(const Vec3 &w) {
      if (!with_angular_velocity_) {
        // store quaternion
        Vec4 q = x_.head(4);
        x_.resize(7);
        x_.head(4) = q;
        x_.tail(3) = w;
        with_angular_velocity_ = true;
      } else {
        x_.tail(3) = w;
      }
    }

    // Attitude Getters
    inline Vec4 GetQuatCoeff() const { return x_.head(4); }
    inline Vec3 GetAngularVelocity() const {
      if (!with_angular_velocity_) {
        throw std::invalid_argument("No angular velocity in the state");
      }
      return x_.tail(3);
    }

    Mat3 GetDCM() const { return QuatCoeffToDCM(GetQuatCoeff()); }
    Vec3 GetEulerAngles() const { return QuatCoeffToEulerAngles(GetQuatCoeff()); }
  };

}  // namespace lupnt

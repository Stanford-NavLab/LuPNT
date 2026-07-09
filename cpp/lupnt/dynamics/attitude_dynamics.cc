#include "lupnt/dynamics/attitude_dynamics.h"

#include "lupnt/conversions/attitude_conversions.h"
#include "lupnt/conversions/frame_converter.h"
#include "lupnt/core/asset_factory.h"
#include "lupnt/interfaces/kernels.h"
#include "lupnt/states/state.h"

namespace lupnt {

  // AttitudeDynamics *************************************************************
  AttitudeDynamics::AttitudeDynamics(Config& config) : Dynamics(config) {}

  State AttitudeDynamics::Propagate(const State& x0, Real t0, Real tf, const State& rv) {
    (void)rv;
    return Propagate(x0, t0, tf, nullptr);
  }

  // FixedAttitudeDynamics *********************************************************
  FixedAttitudeDynamics::FixedAttitudeDynamics(Config& config) : AttitudeDynamics(config) {}

  State FixedAttitudeDynamics::Propagate(const State& x0, Real t0, Real tf, const State* u) {
    (void)t0;
    (void)u;
    (void)tf;
    return x0;
  }

  REGISTER_FACTORY_CLASS(AttitudeDynamics, FixedAttitudeDynamics)

  // FixedPointingDynamics ********************************************************
  FixedPointingDynamics::FixedPointingDynamics(Config& config) {
    primary_axis_ = enum_cast<Axis>(config["primary_axis"].as<std::string>()).value();
    primary_body_id_ = enum_cast<BodyId>(config["primary_body_id"].as<std::string>()).value();
    secondary_axis_ = enum_cast<Axis>(config["secondary_axis"].as<std::string>()).value();
    secondary_body_id_ = enum_cast<BodyId>(config["secondary_body_id"].as<std::string>()).value();
  }

  State FixedPointingDynamics::Propagate(const State& qw, Real t0, Real tf, const State* u) {
    LUPNT_CHECK(false, "Position and velocity required for fixed pointing dynamics",
                "FixedPointingDynamics");
    (void)t0;
    (void)u;
    (void)tf;
    return qw;
  }

  State FixedPointingDynamics::Propagate(const State& qw, Real t0, Real tf, const State& rv) {
    (void)t0;
    Frame frame = qw.GetFrame();
    LUPNT_CHECK(frame == rv.GetFrame(),
                fmt::format("Frame mismatch: {} {}", enum_name(frame), enum_name(rv.GetFrame())),
                "FixedPointingDynamics");
    LUPNT_CHECK(rv.GetType() == Cart6::TYPE,
                fmt::format("rv must be a Cart6 state: {}", rv.GetType()), "FixedPointingDynamics");
    Vec3 r = rv.head(3);

    auto func = [&](Real t) {
      Vec3 dir1 = -r + GetBodyPos(GetLupntEpoch() + t, primary_body_id_, frame);
      Vec3 dir2 = -r + GetBodyPos(GetLupntEpoch() + t, secondary_body_id_, frame);
      int i0 = static_cast<int>(primary_axis_);
      int i1 = static_cast<int>(secondary_axis_);
      int i2 = 3 - i0 - i1;

      Vec3 e0, e1, e2;
      e0 = dir1.normalized();
      e1 = dir2.normalized();
      e2 = (i1 == (i0 + 1) % 3) ? e0.cross(e1) : e1.cross(e0);
      e1 = (i2 == (i0 + 1) % 3) ? e0.cross(e2) : e2.cross(e0);

      Mat3 R_b2w = Mat3::Zero();
      R_b2w.col(i0) = e0.normalized();
      R_b2w.col(i1) = e1.normalized();
      R_b2w.col(i2) = e2.normalized();
      return (Vec9)R_b2w.reshaped();
    };

    Real t = tf.val();
    Mat3 R_b2w = func(tf).reshaped(3, 3);
    Mat3 R_b2w_dot = derivative(func, wrt(t), at(t)).reshaped(3, 3);
    Vec3 w_w = RotToAngularVelocity(R_b2w, R_b2w_dot);
    Vec4 q_b2w = RotToQuat(R_b2w);
    return Attitude(q_b2w, w_w, qw.GetFrame());
  }

  REGISTER_FACTORY_CLASS(AttitudeDynamics, FixedPointingDynamics)

  // Define the GetRegistry function for this specialization (must come before explicit
  // instantiation)
  template <> std::unordered_map<std::string, AssetFactory<AttitudeDynamics, Config&>::Creator>&
  AssetFactory<AttitudeDynamics, Config&>::GetRegistry() {
    return Registry();
  }

  // Explicit template instantiation to ensure single registry across library boundaries
  template class AssetFactory<AttitudeDynamics, Config&>;

}  // namespace lupnt

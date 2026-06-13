#pragma once
#include "lupnt/dynamics/dynamics.h"
#include "lupnt/interfaces/yaml.h"
#include "lupnt/states/state.h"

namespace lupnt {
  /// @brief Base class for attitude (quaternion + angular velocity) propagation
  /// models.
  ///
  /// `AttitudeDynamics` subclasses advance an `Attitude` state (`Attitude::TYPE`:
  /// quaternion `q0..q3` plus body angular velocity `wx,wy,wz`) from `t0` to `tf`.
  /// Used by agents/ to maintain a satellite's/rover's orientation (e.g. for
  /// pointing-dependent device models in devices/ such as cameras and antennas).
  /// The two-argument `Propagate(qw, t0, tf, rv)` overload additionally takes the
  /// object's translational (position/velocity) state when the attitude solution
  /// depends on it (e.g. nadir/sun pointing).
  class AttitudeDynamics : public Dynamics {
  public:
    AttitudeDynamics() = default;

    /// @brief Construct from a YAML configuration node (forwards to
    /// `Dynamics(config)`).
    AttitudeDynamics(Config& config);

    /// @brief Return `Attitude::TYPE`, the quaternion + angular-velocity state type
    /// propagated by attitude dynamics.
    StateType GetStateType() const override { return Attitude::TYPE; }

    using Dynamics::Propagate;

    /// @brief Propagate the attitude state, given the object's translational state.
    ///
    /// Default implementation ignores `rv` and forwards to
    /// `Propagate(qw, t0, tf, nullptr)`; pointing-dependent subclasses (e.g.
    /// FixedPointingDynamics) override this to compute the attitude from `rv`.
    ///
    /// @param qw  Initial attitude state (`Attitude`: quaternion + angular
    ///            velocity) at time `t0`.
    /// @param t0  Initial epoch [s, TDB since J2000].
    /// @param tf  Final epoch [s, TDB since J2000].
    /// @param rv  Object's Cartesian position/velocity state (`Cart6`), used by
    ///            pointing laws that depend on the object's location (e.g. nadir,
    ///            sun, or target pointing).
    /// @return    Propagated attitude state at time `tf`.
    virtual State Propagate(const State& qw, Real t0, Real tf, const State& rv);
  };

  /// @brief Attitude dynamics that holds the attitude state constant (no rotation).
  ///
  /// Used as a placeholder/default attitude model for objects whose orientation is
  /// not tracked or is irrelevant to the simulation (e.g. point-mass agents).
  class FixedAttitudeDynamics : public AttitudeDynamics {
  public:
    FixedAttitudeDynamics() = default;

    /// @brief Construct from a YAML configuration node (forwards to
    /// `AttitudeDynamics(config)`).
    FixedAttitudeDynamics(Config& config);

    using AttitudeDynamics::Propagate;

    /// @brief Return the attitude state unchanged (`xf = x0`).
    State Propagate(const State& qw, Real t0, Real tf, const State* u = nullptr) override;

    /// @brief Return `Attitude::TYPE`.
    StateType GetStateType() const override { return Attitude::TYPE; }
  };

  /// @brief Attitude dynamics that orients a body frame so two configured axes point
  /// at two configured target bodies ("two-vector" / fixed-pointing attitude).
  ///
  /// Builds a body-to-world rotation whose `primary_axis_` points toward
  /// `primary_body_id_` and whose `secondary_axis_` is constrained toward
  /// `secondary_body_id_` (e.g. antenna boresight at the target body, solar panel
  /// normal at the Sun). The resulting quaternion and angular velocity (obtained by
  /// differentiating the rotation in time via autodiff) are returned as an
  /// `Attitude` state. Used by agents/ (e.g. constellation satellites, see
  /// `Constellation` construction in agents/constellation.cc) to drive
  /// pointing-dependent devices.
  class FixedPointingDynamics : public AttitudeDynamics {
  private:
    Axis primary_axis_, secondary_axis_;
    BodyId primary_body_id_, secondary_body_id_;

  public:
    FixedPointingDynamics() = default;

    /// @brief Construct from a YAML configuration node specifying the primary and
    /// secondary pointing axes and target bodies.
    ///
    /// Reads `primary_axis`, `primary_body_id`, `secondary_axis`, and
    /// `secondary_body_id` from `config`.
    FixedPointingDynamics(Config& config);

    /// @brief Set the primary pointing axis and the body it should track.
    ///
    /// @param axis     Body-frame axis (X/Y/Z) to align with the direction to
    ///                 `body_id`.
    /// @param body_id  Target body whose direction defines the primary pointing
    ///                 axis.
    void SetPrimaryPointing(Axis axis, BodyId body_id) {
      primary_axis_ = axis;
      primary_body_id_ = body_id;
    }

    /// @brief Set the secondary pointing axis and the body it should track.
    ///
    /// @param axis     Body-frame axis (X/Y/Z) constrained toward the direction to
    ///                 `body_id`.
    /// @param body_id  Target body whose direction defines the secondary pointing
    ///                 axis.
    void SetSecondaryPointing(Axis axis, BodyId body_id) {
      secondary_axis_ = axis;
      secondary_body_id_ = body_id;
    }

    using AttitudeDynamics::Propagate;

    /// @brief Compute the two-vector pointing attitude at time `tf` from the
    /// object's position.
    ///
    /// Builds an orthonormal body-to-world rotation matrix `R_b2w(t)` whose columns
    /// align the primary/secondary body axes with the directions from the object's
    /// position `r = rv.head(3)` to `primary_body_id_`/`secondary_body_id_` at
    /// `GetLupntEpoch() + t`. The angular velocity is obtained from the analytic
    /// time-derivative of `R_b2w` (via autodiff) using `RotToAngularVelocity`.
    ///
    /// @param qw  Initial attitude state at `t0` (only its frame is used).
    /// @param t0  Initial epoch [s, TDB since J2000] (unused; attitude is computed
    ///            directly at `tf`).
    /// @param tf  Final epoch [s, TDB since J2000] at which to evaluate the
    ///            pointing attitude.
    /// @param rv  Object's Cartesian position/velocity state (`Cart6`); must share
    ///            the same frame as `qw` and have position in `rv.head(3)`.
    /// @return    Attitude state (quaternion + body angular velocity) realizing the
    ///            two-vector pointing law at time `tf`.
    State Propagate(const State& qw, Real t0, Real tf, const State& rv) override;

    /// @brief Unimplemented: fixed-pointing attitude requires the object's
    /// position/velocity.
    ///
    /// Always raises a `LUPNT_CHECK` failure ("Position and velocity required for
    /// fixed pointing dynamics") and returns `qw` unchanged. Callers must use the
    /// `Propagate(qw, t0, tf, rv)` overload instead.
    State Propagate(const State& qw, Real t0, Real tf, const State* u) override;

    /// @brief Return `Attitude::TYPE`.
    StateType GetStateType() const override { return Attitude::TYPE; }
  };
}  // namespace lupnt

#pragma once

#include "lupnt/dynamics/dynamics.h"
#include "lupnt/states/state.h"

namespace lupnt {

  /// @brief Base class for closed-form (non-integrated) orbit propagation models.
  ///
  /// `AnalyticalDynamics` subclasses advance an orbital-element or relative-motion
  /// state from `t0` to `tf` using an analytic state-transition matrix evaluated at
  /// `dt = tf - t0`, rather than numerically integrating an ODE (contrast with
  /// NumericalDynamics). Used wherever a fast/closed-form propagator is sufficient,
  /// e.g. Keplerian two-body propagation, Clohessy-Wiltshire relative motion, or
  /// formation-flying geometric mappings.
  class AnalyticalDynamics : public Dynamics {
  public:
    virtual ~AnalyticalDynamics() = default;

    using Dynamics::Propagate;

    /// @brief Propagate the state from `t0` to `tf` without computing a state
    /// transition matrix.
    ///
    /// Forwards to `Propagate(x0, t0, tf, u, nullptr)`.
    ///
    /// @param x0  Initial state at time `t0`.
    /// @param t0  Initial epoch [s].
    /// @param tf  Final epoch [s].
    /// @param u   Optional control/forcing input; unused by analytic propagators.
    /// @return    Propagated state at time `tf`.
    State Propagate(const State &x0, Real t0, Real tf, const State *u = nullptr) override;

    /// @brief Propagate the state to a sequence of output epochs, all relative to
    /// the initial epoch `tfs(0)`.
    ///
    /// Unlike the generic `Dynamics::Propagate(x0, tfs, u)`, each output row is
    /// computed directly as `Propagate(x0, tfs(0), tfs(i), u)` (parallelized with
    /// OpenMP), exploiting the fact that analytic propagators can jump to any `tf`
    /// in one step without re-using the previous step's result.
    ///
    /// @param x0   Initial state at `tfs(0)`.
    /// @param tfs  Vector of output epochs [s], including the initial epoch as
    ///             `tfs(0)`.
    /// @param u    Optional control/forcing input; unused by analytic propagators.
    /// @return     Matrix whose i-th row is the propagated state at `tfs(i)`.
    MatX Propagate(const State &x0, const VecX &tfs, const State *u = nullptr) override;

    /// @brief Propagate the state from `t0` to `tf` and, if requested, compute the
    /// analytic state transition matrix d(xf)/d(x0).
    ///
    /// Implemented by each concrete analytical propagator (KeplerianDynamics,
    /// ClohessyWiltshireDynamics, YamanakaAnkersenDynamics,
    /// RoeGeometricMappingDynamics) using its closed-form STM.
    ///
    /// @param x0   Initial state at time `t0`.
    /// @param t0   Initial epoch [s].
    /// @param tf   Final epoch [s].
    /// @param u    Optional control/forcing input; unused by analytic propagators.
    /// @param stm  Output: analytic state transition matrix d(xf)/d(x0); left
    ///             unmodified if nullptr.
    /// @return     Propagated state at time `tf`.
    virtual State Propagate(const State &x0, Real t0, Real tf, const State *u, MatXd *stm) override
        = 0;
  };

  /// @brief Keplerian (unperturbed two-body) propagation of an orbital-element state.
  ///
  /// Advances the mean/true anomaly-like element (depending on `T`) by the mean
  /// motion `n = sqrt(GM_ / a^3) * dt`, holding all other elements fixed. Used as a
  /// lightweight ground-truth or reference propagator for `ClassicalOE`,
  /// `QuasiNonsingularOE`, or `EquinoctialOE` states (e.g. in tests and analytic
  /// orbit examples) where high-fidelity perturbations are not needed.
  template <typename T = ClassicalOE> class KeplerianDynamics : public AnalyticalDynamics {
  private:
    Real GM_;

  public:
    /// @brief Construct with the gravitational parameter of the central body.
    ///
    /// @param GM  Gravitational parameter [m^3/s^2 or consistent unit system].
    KeplerianDynamics(Real GM) : GM_(GM) {}
    using AnalyticalDynamics::Propagate;

    /// @brief Propagate the orbital-element state by Keplerian mean motion over
    /// `dt = tf - t0`.
    ///
    /// Only the anomaly-like element (mean anomaly for `ClassicalOE`, mean argument
    /// of latitude for `QuasiNonsingularOE`, mean longitude for `EquinoctialOE`) is
    /// advanced; `stm` (when requested) is the identity plus, for `ClassicalOE`,
    /// the analytic d(M)/d(a) sensitivity term. STM computation is not implemented
    /// for the `QuasiNonsingularOE`/`EquinoctialOE` specializations.
    ///
    /// @param state  Initial orbital-element state at time `t0`.
    /// @param t0     Initial epoch [s].
    /// @param tf     Final epoch [s].
    /// @param u      Unused (no control input for Keplerian motion).
    /// @param stm    Output: state transition matrix d(xf)/d(x0); left unmodified if
    ///               nullptr.
    /// @return       Propagated orbital-element state at time `tf`.
    State Propagate(const State &state, Real t0, Real tf, const State *u, MatXd *stm) override;

    /// @brief Return the orbital-element `StateType` (`T::TYPE`) propagated by this
    /// instantiation.
    StateType GetStateType() const override { return T::TYPE; }
  };

  /// @brief Clohessy-Wiltshire (linearized two-body relative motion) propagation of a
  /// chief-relative Cartesian state.
  ///
  /// Propagates a `RelCart6` (relative position/velocity, e.g. in the chief's RTN
  /// frame) using the analytic CW state transition matrix for a circular reference
  /// orbit of semi-major axis `a` and mean motion `n`. Used for fast
  /// formation-flying / rendezvous relative-motion propagation where the CW
  /// linearization (circular reference, small separation) is adequate.
  class ClohessyWiltshireDynamics : public AnalyticalDynamics {
  private:
    Real a_, n_;
    VecX K_;
    Real t0_;

  public:
    /// @brief Construct with the chief orbit's semi-major axis and mean motion.
    ///
    /// @param a  Chief orbit semi-major axis [m or consistent length unit].
    /// @param n  Chief orbit mean motion [rad/s].
    ClohessyWiltshireDynamics(Real a, Real n);

    /// @brief Compute the Clohessy-Wiltshire state transition matrix for elapsed
    /// time `tf`.
    ///
    /// Returns `Phi = A * B(n_*tf)` mapping the integration-constant vector `K_`
    /// (solved for at the first call to Propagate) to the relative state at
    /// `t0_ + tf`.
    ///
    /// @param tf  Elapsed time since `t0_` [s].
    /// @return    6x6 CW state transition matrix.
    Mat6 ComputeMat(Real tf);
    using AnalyticalDynamics::Propagate;

    /// @brief Propagate the relative Cartesian state using the Clohessy-Wiltshire
    /// solution.
    ///
    /// On the first call (or whenever `t0` changes), solves for the CW integration
    /// constants `K_` from `x0` via `ComputeMat(t0)`, then evaluates
    /// `xf = ComputeMat(tf - t0) * K_`.
    ///
    /// @param state  Initial relative Cartesian state (`RelCart6`) at time `t0`.
    /// @param t0     Initial epoch [s].
    /// @param tf     Final epoch [s].
    /// @param u      Unused (no control input).
    /// @param stm    Output: CW state transition matrix d(xf)/d(x0); left unmodified
    ///               if nullptr.
    /// @return       Propagated relative Cartesian state at time `tf`.
    State Propagate(const State &state, Real t0, Real tf, const State *u, MatXd *stm) override;

    /// @brief Return `RelCart6::TYPE`, the relative-position/velocity state type
    /// propagated by this model.
    StateType GetStateType() const override { return RelCart6::TYPE; }
  };

  /// @brief Yamanaka-Ankersen propagation of relative motion about an eccentric
  /// reference orbit.
  ///
  /// Extends Clohessy-Wiltshire to eccentric chief orbits using the
  /// Yamanaka-Ankersen state transition matrix, parameterized by the true anomaly
  /// `f` of the chief. Intended for formation-flying relative-motion propagation
  /// when the chief orbit is non-circular.
  ///
  /// @note Propagate() currently throws `std::runtime_error("Not implemented")`;
  /// only ComputeMat/ComputeInverseMat are functional.
  class YamanakaAnkersenDynamics : public AnalyticalDynamics {
  private:
    Real a_, n_, e_, M0_;
    VecX K_;
    Real t0_;
    VecX rv_rtn_;

  public:
    /// @brief Construct from the chief's classical orbital elements and the
    /// deputy's initial RTN relative state.
    ///
    /// @param coe_c  Chief's classical orbital elements (provides semi-major axis,
    ///               eccentricity, mean anomaly).
    /// @param rv_rtn Deputy's initial relative position/velocity in the chief's RTN
    ///               frame.
    /// @param GM_    Gravitational parameter of the central body [m^3/s^2 or
    ///               consistent unit system].
    YamanakaAnkersenDynamics(const ClassicalOE &coe_c, const Cart6 &rv_rtn, Real GM_);

    /// @brief Compute the Yamanaka-Ankersen state transition matrix at elapsed time
    /// `t`.
    ///
    /// Evaluates the true anomaly `f` at `t` via `MeanToTrueAnomaly`, then forms
    /// `Phi = A(f) * B(f)` mapping the (as-yet-unsolved) integration-constant vector
    /// to the relative state.
    ///
    /// @param t  Elapsed time since `t0_` [s].
    /// @return   6x6 Yamanaka-Ankersen state transition matrix.
    MatX ComputeMat(Real t);

    /// @brief Compute the inverse of the Yamanaka-Ankersen state transition matrix
    /// at elapsed time `t`.
    ///
    /// Used to map a relative state at time `t` back to the integration-constant
    /// vector (`Phi^{-1} * x`), e.g. for initializing `K_` from `rv_rtn_`.
    ///
    /// @param t  Elapsed time since `t0_` [s].
    /// @return   6x6 inverse Yamanaka-Ankersen state transition matrix.
    MatX ComputeInverseMat(Real t);
    using AnalyticalDynamics::Propagate;

    /// @brief Propagate the relative Cartesian state using the Yamanaka-Ankersen
    /// solution.
    ///
    /// @note Currently unimplemented (throws `std::runtime_error`) for `dt != 0`.
    ///
    /// @param state  Initial relative Cartesian state at time `t0`.
    /// @param t0     Initial epoch [s].
    /// @param tf     Final epoch [s].
    /// @param u      Unused (no control input).
    /// @param stm    Output: state transition matrix d(xf)/d(x0); left unmodified if
    ///               nullptr.
    /// @return       Propagated relative Cartesian state at time `tf`.
    State Propagate(const State &state, Real t0, Real tf, const State *u, MatXd *stm) override;

    /// @brief Return `ClassicalOE::TYPE`, the state type associated with this model.
    StateType GetStateType() const override { return ClassicalOE::TYPE; }
  };

  /// @brief Geometric-mapping propagation of quasi-nonsingular relative orbital
  /// elements (ROE).
  ///
  /// Maps a fixed quasi-nonsingular ROE vector `K_` to a relative Cartesian state at
  /// time `t` via the geometric ROE-to-RTN transition matrix, parameterized by the
  /// chief's argument of latitude `u = f + w`. Intended for formation-flying
  /// analysis expressed in ROE coordinates.
  ///
  /// @note Propagate() currently throws `std::runtime_error("Not implemented")`;
  /// only ComputeMat is functional.
  class RoeGeometricMappingDynamics : public AnalyticalDynamics {
  private:
    Real a_, e_, i_, w_, M0_;
    Real ex_, ey_, n_;
    VecX K_;
    Real t0_;

  public:
    /// @brief Construct from the chief's classical orbital elements, the deputy's
    /// quasi-nonsingular ROE, and the central-body gravitational parameter.
    ///
    /// @param coe_c  Chief's classical orbital elements at the reference epoch.
    /// @param roe    Deputy's quasi-nonsingular relative orbital elements (ROE),
    ///               stored as the fixed integration-constant vector `K_`.
    /// @param GM     Gravitational parameter of the central body [m^3/s^2 or
    ///               consistent unit system].
    RoeGeometricMappingDynamics(const ClassicalOE coe_c, const QuasiNonsingROE &roe, Real GM);

    /// @brief Compute the ROE-to-relative-Cartesian geometric mapping matrix at
    /// elapsed time `t`.
    ///
    /// Evaluates the chief's true anomaly `f` at `t` via `MeanToTrueAnomaly`, forms
    /// the argument of latitude `u = f + w_`, and returns `Phi = A(u) * B(u)`
    /// mapping the fixed ROE vector `K_` to the relative Cartesian state.
    ///
    /// @param t  Elapsed time since `t0_` [s].
    /// @return   6x6 ROE geometric mapping matrix.
    MatX ComputeMat(Real t);
    using AnalyticalDynamics::Propagate;

    /// @brief Propagate the relative Cartesian state using the ROE geometric
    /// mapping.
    ///
    /// @note Currently unimplemented (throws `std::runtime_error`) for `dt != 0`.
    ///
    /// @param state  Initial relative Cartesian state at time `t0`.
    /// @param t0     Initial epoch [s].
    /// @param tf     Final epoch [s].
    /// @param u      Unused (no control input).
    /// @param stm    Output: state transition matrix d(xf)/d(x0); left unmodified if
    ///               nullptr.
    /// @return       Propagated relative Cartesian state at time `tf`.
    State Propagate(const State &state, Real t0, Real tf, const State *u, MatXd *stm) override;

    /// @brief Return `QuasiNonsingROE::TYPE`, the relative-orbital-element state
    /// type propagated by this model.
    StateType GetStateType() const override { return QuasiNonsingROE::TYPE; }
  };

}  // namespace lupnt

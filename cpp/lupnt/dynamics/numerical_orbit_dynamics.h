
/**
 * @file numerical_orbit_dynamics.h
 * @brief Header file for numerical orbit dynamics
 *
 */

#pragma once

#include "lupnt/dynamics/dynamics.h"
#include "lupnt/environment/body.h"
#include "lupnt/numerics/integrator.h"
#include "lupnt/states/params.h"
#include "lupnt/states/state.h"

namespace lupnt {

  struct TerminationInfo {
    bool terminated = false;  // true if stopped before tf due to user/max-iter
    TerminationReason reason = TerminationReason::ReachedTf;
    Real t_stop = 0.0;  // time actually reached
    int steps = 0;      // steps taken
  };

  /**
   * @class NumericalOrbitDynamics
   * @brief Class for numerical orbit dynamics.
   *
   * This class provides methods for propagating the state vector of an orbit
   * numerically using an ODE function and a numerical propagator.
   */
  class NumericalDynamics : public Dynamics {
  private:
    ODE odefunc_;                 ///< ODE function for the dynamics.
    Ptr<Integrator> integrator_;  ///< Numerical propagator.
    Real dt_ = 10.0;              ///< Time step for numerical integration.

  public:
    NumericalDynamics();
    NumericalDynamics(Config& dynamics_config);

    void SetTimeStep(Real dt);
    Real GetTimeStep() const;

    void SetODE(ODE odefunc);
    void SetIntegrator(IntegratorType integ);
    void SetIntegratorParams(IntegratorParams params);

    virtual VecX ComputeRates(Real t, const State& x) const;

    Integrator* GetIntegrator() { return integrator_.get(); }

    using Dynamics::Propagate;
    State Propagate(const State& x0, Real t0, Real tf, const State* u = nullptr) override;
    // State Propagate(const State &x0, Real t0, Real tf, const State *u, MatXd *stm) override;
    StateType GetStateType() const override { return Cart6::TYPE; }

    State PropagateEx(const State& x0, Real t0, Real tf, TerminationInfo* info);
    MatX PropagateEx(const State& x0, Real t0, const VecX& tf, TerminationInfo* info);
    State PropagateExStm(const State& x0, Real t0, Real tf, MatXd* stm, TerminationInfo* info);
  };

  // Cartesian Two-Body Dynamics
  class CartesianTwoBodyDynamics : public NumericalDynamics {
  private:
    Real GM_;

  public:
    CartesianTwoBodyDynamics(Real GM);
    VecX ComputeRates(Real t, const State& x) const override;
    using NumericalDynamics::Propagate;
    StateType GetStateType() const override { return Cart6::TYPE; }
  };

  // J2 Cartesian Two-Body Dynamics
  class JToCartTwoBodyDynamics : public NumericalDynamics {
  private:
    Real GM_, J2_, R_body_;
    Frame frame_ = Frame::GCRF;
    Frame body_fixed_frame_ = Frame::ITRF;

  public:
    JToCartTwoBodyDynamics(Real GM, Real J2, Real R_body, Frame frame = Frame::GCRF,
                           Frame body_fixed_frame = Frame::ITRF);
    VecX ComputeRates(Real t, const State& x) const override;
    using NumericalDynamics::Propagate;
    StateType GetStateType() const override { return Cart6::TYPE; }
    void SetFrame(Frame frame) { frame_ = frame; }
    Frame GetFrame() const { return frame_; }
    void SetBodyFixedFrame(Frame frame) { body_fixed_frame_ = frame; }
    Frame GetBodyFixedFrame() const { return body_fixed_frame_; }
  };

  class J2KeplerianDynamics : public NumericalDynamics {
  private:
    Real GM_, J2_, R_body_;

  public:
    J2KeplerianDynamics(Real GM, Real J2, Real R_body);
    VecX ComputeRates(Real t, const State& x) const override;
    using NumericalDynamics::Propagate;
    StateType GetStateType() const override { return Cart6::TYPE; }
  };

  class MoonMeanDynamics : public NumericalDynamics {
  private:
    const double n3_ = 2.66e-6;
    const double J2_ = 2.03e-4;
    const double k_ = 0.98785;

  public:
    MoonMeanDynamics();
    VecX ComputeRates(Real t, const State& x) const override;
    using NumericalDynamics::Propagate;
    StateType GetStateType() const override { return Cart6::TYPE; }
  };

  /**
   * @brief Numerical point-mass, gravity-field, SRP, drag, and relativity orbit dynamics.
   *
   * The model propagates Cartesian position and velocity in the selected integration
   * frame. Body ephemerides are obtained in SI units internally and scaled to the
   * configured UnitSystem at the dynamics boundary, so the propagated state, body
   * constants, and force-model parameters remain unit-consistent.
   */
  class NBodyDynamics : public NumericalDynamics {
  private:
    std::vector<Body> bodies_;
    bool use_srp_ = false;
    bool use_drag_ = false;
    bool use_relativity_ = true;
    bool use_ad_ = false;
    Frame frame_ = Frame::MOON_CI;
    UnitSystem units_ = SI_UNITS;

  public:
    NBodyDynamics();
    /**
     * @brief Construct an n-body dynamics model from a YAML configuration.
     *
     * Supported configuration fields include the integration frame, force-model
     * switches, SRP and drag ballistic coefficients, and the unit system.
     */
    NBodyDynamics(Config& dynamics_config);

    // Propagators
    using NumericalDynamics::Propagate;
    /**
     * @brief Propagate a Cartesian state and optionally compute the STM.
     *
     * @param x0 Initial Cartesian state in the configured frame and unit system.
     * @param t0 Initial TDB epoch in seconds.
     * @param tf Final TDB epoch in seconds.
     * @param u Optional control input.
     * @param stm Optional state transition matrix output.
     * @return Final Cartesian state in the configured frame and unit system.
     */
    State Propagate(const State& x0, Real t0, Real tf, const State* u, MatXd* stm) override;
    /**
     * @brief Compute Cartesian position and velocity rates.
     *
     * @param t TDB epoch in seconds.
     * @param x Cartesian state in the configured frame and unit system.
     * @return Time derivative of the state in the configured unit system.
     */
    VecX ComputeRates(Real t, const State& x) const override;

    // Modifiers
    /**
     * @brief Add a perturbing or central body to the dynamics model.
     *
     * Bodies should use the same UnitSystem as the dynamics model. SetUnits()
     * updates already-added built-in bodies to the requested units.
     */
    void AddBody(const Body& body);
    /**
     * @brief Remove a body from the dynamics model by matching its BodyId.
     */
    void RemoveBody(const Body& body);

    // Getters
    /// @brief Return true when solar radiation pressure is enabled.
    bool GetUseSrp() { return use_srp_; }
    /// @brief Return true when atmospheric drag is enabled.
    bool GetUseDrag() { return use_drag_; }
    /// @brief Return true when post-Newtonian relativistic corrections are enabled.
    bool GetUseRelativity() const { return use_relativity_; }
    /// @brief Return the SRP ballistic coefficient used by the force model.
    Real GetSrpCoeff() const { return GetParam("bcoeff_srp"); }
    /// @brief Return the drag ballistic coefficient used by the force model.
    Real GetDragCoeff() const { return GetParam("bcoeff_drag"); }
    StateType GetStateType() const override { return Cart6::TYPE; }
    /// @brief Return the integration frame.
    Frame GetFrame() { return frame_; }
    /// @brief Return the bodies currently included in the dynamics model.
    std::vector<Body> GetBodies() { return bodies_; }
    /// @brief Return the active distance, time, and mass unit system.
    UnitSystem GetUnits() const { return units_; }

    // Setters
    /// @brief Enable or disable automatic differentiation rates.
    void SetAutodiff(bool use_ad) { use_ad_ = use_ad; }
    /// @brief Set the SRP ballistic coefficient directly.
    void SetSrpCoeff(Real bcoeff);
    /// @brief Set the drag ballistic coefficient directly.
    void SetDragCoeff(Real bcoeff);
    /// @brief Set SRP coefficient from reflectivity coefficient, area, and mass.
    void SetSrpCoefficient(Real CR, Real area, Real mass);
    /// @brief Set drag coefficient from drag coefficient, area, and mass.
    void SetDragCoefficient(Real CD, Real area, Real mass);
    /// @brief Enable or disable solar radiation pressure.
    void SetUseSrp(bool use_srp) { use_srp_ = use_srp; }
    /// @brief Enable or disable atmospheric drag.
    void SetUseDrag(bool use_drag) { use_drag_ = use_drag; }
    /// @brief Enable or disable post-Newtonian relativistic corrections.
    void SetUseRelativity(bool use_relativity) { use_relativity_ = use_relativity; }
    /// @brief Set the frame used for numerical integration.
    void SetFrame(Frame frame) { frame_ = frame; }
    /**
     * @brief Set the active unit system for states, constants, and acceleration output.
     *
     * Built-in bodies already registered with the model are recreated in the new
     * unit system. Time epochs remain TDB seconds.
     */
    void SetUnits(const UnitSystem& units);
  };

}  // namespace lupnt

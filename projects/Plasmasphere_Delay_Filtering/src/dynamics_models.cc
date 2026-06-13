#include "src/dynamics_models.h"

#include <lupnt/lupnt.h>

namespace filtering_sim {
  using namespace lupnt;

  /****************************************************************************
   * Joint Orbit-Clock State
   ****************************************************************************/
  JointOrbitClockState::JointOrbitClockState(const State& orbit_state, const State& clock_state) {
    clock_state_size_ = clock_state.size();
    resize(orbit_state_size_ + clock_state_size_);
    head(orbit_state_size_) = orbit_state;
    tail(clock_state_size_) = clock_state;
    Init(clock_state_size_);
    SetFrame(orbit_state.GetFrame());  // Set frame to orbit state's frame
  }

  /****************************************************************************
   * Joint Orbit-Clock Dynamics
   ****************************************************************************/
  void JointOrbitClockDynamics::SetOrbitDynamics(Ptr<NBodyDynamics> orbit_dynamics) {
    orbit_dynamics_ = orbit_dynamics;
    // Update params
    SetFrame(orbit_dynamics_->GetFrame());
    BodyId center = GetFrameCenter(frame_);
    BodyData center_body_data = GetBodyData(center);
    GM_ = center_body_data.GM;

    SetParams(orbit_dynamics_->GetParams());
  }

  void JointOrbitClockDynamics::SetClockDynamics(Ptr<ClockDynamics> clock_dynamics) {
    clock_dynamics_ = clock_dynamics;
  }

  // Constructor ----------------------------------------------------------------
  JointOrbitClockDynamics::JointOrbitClockDynamics() : NumericalDynamics() {
    SetODE([this](Real t, const VecX& x) { return ComputeRates(t, x); });
    ParamState default_params = ParamState(VecX::Zero(2), {"bcoeff_srp", "bcoeff_drag"});
    this->SetParams(default_params);
  };

  JointOrbitClockDynamics::JointOrbitClockDynamics(Config& config) : NumericalDynamics(config) {
    SetODE([this](Real t, const VecX& x) { return ComputeRates(t, x); });
    Logger::Debug("Creating", "JointOrbitClockDynamics");
    ParamState default_params = ParamState(VecX::Zero(2), {"bcoeff_srp", "bcoeff_drag"});
    this->SetParams(default_params);
  }

  REGISTER_FACTORY_CLASS(Dynamics, JointOrbitClockDynamics)

  // Overwrite Propagaters ----------------------------------------------------------------
  State JointOrbitClockDynamics::PropagateWithParams(const State& x0, Real t0, Real tf,
                                                     const ParamState& params, const State* u) {
    // Set parameters to orbit dynamics
    orbit_dynamics_->SetParams(params);

    // Call base class Propagate
    return JointOrbitClockDynamics::Propagate(x0, t0, tf, u);
  }

  State JointOrbitClockDynamics::Propagate(const State& x0, Real t0, Real tf, const State* u) {
    // Call base class Propagate
    State x_prop = NumericalDynamics::Propagate(x0, t0, tf, u);

    if (add_noise_) {
      // Add process noise to clock states
      int clock_state_size_ = x0.size() - orbit_state_size_;
      double dt = tf - t0;
      VecX clock_process_noise = clock_dynamics_->GetProcessNoise(dt, clock_state_size_);
      x_prop.tail(clock_state_size_) += C * clock_process_noise;
    }

    return x_prop;
  }

  // Compute Rates ----------------------------------------------------------------
  VecX JointOrbitClockDynamics::ComputeRates(Real t, const State& x) const {
    // Split state into orbit and clock states
    int clock_state_size_ = x.size() - orbit_state_size_;
    State orbit_state = x.head(orbit_state_size_);
    State clock_state = x.tail(clock_state_size_);

    // Obtain the frame of the stat
    // Compute relativistic delays
    Real r = orbit_state.head(3).norm();
    Real v = orbit_state.segment(3, 3).norm();

    // Relativisitc clock drift
    // The clock state is in meters -> devide by C instead of C^2
    Real rel_drift = -(GM_ / r + v * v / 2) / (C);

    // Todo: add sagnac if the frame is rotating

    // Orbit rates
    VecX orbit_rates = orbit_dynamics_->ComputeRates(t, orbit_state);

    // Clock rates
    VecX clock_rates(clock_state_size_);
    if (clock_state_size_ == 2) {
      clock_rates(0) = clock_state(1);
      clock_rates(1) = rel_drift;
    } else {
      clock_rates(0) = clock_state(1);
      clock_rates(1) = clock_state(2) + rel_drift;
      clock_rates(2) = 0.0;
    }

    // Combine rates
    VecX rates(x.size());
    rates << orbit_rates, clock_rates;
    return rates;
  }

  /****************************************************************************
   * Clock Noise Matrices
   ****************************************************************************/
  Mat3 ThreeStateClockNoise(ClockModel clk_model, Real dt, double inflate_q2) {
    Mat3 Q(3, 3);
    double q1, q2, q3;
    std::tie(q1, q2, q3) = ClockDynamics::GetClockValues(clk_model);
    q2 *= inflate_q2;

    Q(0, 0) = q1 * dt + q2 * pow(dt, 3) / 3 + q3 * pow(dt, 5) / 20;
    Q(0, 1) = q2 * pow(dt, 2) / 2 + q3 * pow(dt, 4) / 8;
    Q(0, 2) = q3 * pow(dt, 3) / 6;
    Q(1, 0) = Q(0, 1);
    Q(1, 1) = q2 * dt + q3 * pow(dt, 3) / 3;
    Q(1, 2) = q3 * pow(dt, 2) / 2;
    Q(2, 0) = Q(0, 2);
    Q(2, 1) = Q(1, 2);
    Q(2, 2) = q3 * dt;
    return Q;
  }

  Mat2 TwoStateClockNoise(ClockModel clk_model, Real dt, double inflate_q2) {
    Mat2 Q(2, 2);
    double q1, q2, q3;
    std::tie(q1, q2, q3) = ClockDynamics::GetClockValues(clk_model);
    q2 *= inflate_q2;

    Q(0, 0) = q1 * dt + q2 * pow(dt, 3) / 3;
    Q(0, 1) = q2 * pow(dt, 2) / 2;
    Q(1, 0) = Q(0, 1);
    Q(1, 1) = q2 * dt;
    return Q;
  }

}  // namespace filtering_sim

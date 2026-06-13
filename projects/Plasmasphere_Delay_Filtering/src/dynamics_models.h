#pragma once

#include <lupnt/lupnt.h>

namespace filtering_sim {
  using namespace lupnt;

  class JointOrbitClockState : public State {
  private:
    int orbit_state_size_ = 6;
    int clock_state_size_ = 3;

    void Init(int clock_size) {
      SetType(JointOrbitClockState::TYPE);
      if (clock_size == 2) {
        SetNames({"x", "y", "z", "vx", "vy", "vz", "b", "d"});
        SetUnits({"m", "m", "m", "m/s", "m/s", "m/s", "s", "s/s"});
      } else if (clock_size == 3) {
        SetNames({"x", "y", "z", "vx", "vy", "vz", "b", "d", "dr"});
        SetUnits({"m", "m", "m", "m/s", "m/s", "m/s", "s", "s/s", "s/s^2"});
      } else {
        throw std::runtime_error("Unsupported clock state size");
      }
      clock_state_size_ = clock_size;
    }

  public:
    static constexpr StateType TYPE = "JointOrbitClock";
    JointOrbitClockState() : State() { Init(3); }
    JointOrbitClockState(int clock_size) : State() { Init(clock_size); }
    JointOrbitClockState(const State& x) : State() {
      if (x.size() < orbit_state_size_ + 2 || x.size() > orbit_state_size_ + 3) {
        throw std::runtime_error("Invalid state size for JointOrbitClockState");
      }
      orbit_state_size_ = 6;
      clock_state_size_ = x.size() - orbit_state_size_;
      SetType(JointOrbitClockState::TYPE);
    }
    JointOrbitClockState(const State& orbit_state, const State& clock_state);

    JointOrbitClockState& operator=(const State& x) {
      if (this != &x) {
        if (x.size() < orbit_state_size_ + 2 || x.size() > orbit_state_size_ + 3) {
          throw std::runtime_error("Invalid state size for JointOrbitClockState");
        }
        head(orbit_state_size_ + clock_state_size_) = x.head(orbit_state_size_ + clock_state_size_);
      }
      return *this;
    }

    ~JointOrbitClockState() = default;

    // Getters
    State GetOrbitState() const;
    State GetClockState() const;

    // Setters
    void SetOrbitState(const State& orbit_state);
    void SetClockState(const State& clock_state);
  };

  class JointOrbitClockDynamics : public NumericalDynamics {
  private:
    Ptr<NBodyDynamics> orbit_dynamics_;
    Ptr<ClockDynamics> clock_dynamics_;
    int orbit_state_size_ = 6;
    int clock_state_size_ = 3;
    State orbit_state_;
    State clock_state_;
    Frame frame_ = Frame::UNDEFINED;
    double GM_ = GM_MOON;  // Default to Moon GM
    bool add_noise_ = false;

  public:
    JointOrbitClockDynamics();
    JointOrbitClockDynamics(Config& config);

    ~JointOrbitClockDynamics() = default;

    // Setters
    void SetOrbitDynamics(Ptr<NBodyDynamics> orbit_dynamics);
    void SetClockDynamics(Ptr<ClockDynamics> clock_dynamics);
    void SetFrame(Frame frame) { frame_ = frame; }
    void SetAddNoise(bool add_noise) { add_noise_ = add_noise; }

    // Getters
    StateType GetStateType() const override { return JointOrbitClockState::TYPE; }
    Frame GetFrame() const { return frame_; }

    using NumericalDynamics::PropagateWithParams;
    State PropagateWithParams(const State& x0, Real t0, Real tf, const ParamState& params,
                              const State* u = nullptr) override;

    using NumericalDynamics::Propagate;
    State Propagate(const State& x0, Real t0, Real tf, const State* u = nullptr) override;

    VecX ComputeRates(Real t, const State& x) const override;
  };

  Mat3 ThreeStateClockNoise(ClockModel clk_model, Real dt, double inflate_q2);
  Mat2 TwoStateClockNoise(ClockModel clk_model, Real dt, double inflate_q2);

}  // namespace filtering_sim

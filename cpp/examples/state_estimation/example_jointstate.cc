// Demonstrates JointState construction for a combined lunar orbit, SRP
// parameter, and clock state, then propagates the joint STM for one minute.
#include <lupnt/lupnt.h>

using namespace lupnt;

int main() {
  Real t0 = GregorianToTime(2025, 7, 15, 1, 0, 0);  // [s] Start time (TAI)

  // Orbital elements
  Real a = 6541.4e3;    // [m] Semi-major axis
  Real e = 0.60;        // [-] Eccentricity
  Real i = 56.2 * RAD;  // [deg] Inclination
  Real O = 0.00 * RAD;  // [deg] Right ascension of the ascending node
  Real w = 90.0 * RAD;  // [deg] Argument of perigee
  Real M = 0.00 * RAD;  // [deg] Mean anomaly

  Vec6 coe0_op(a, e, i, O, w, M);
  ClassicalOE sat_coe = ClassicalOE(coe0_op, Frame::MOON_OP);
  Cart6 rv0_op = ClassicalToCart(coe0_op, GM_MOON);
  Cart6 rv0_ci = ConvertFrame(t0, rv0_op, Frame::MOON_CI);

  // Clock state
  ClockModel clock_model = ClockModel::OCXO;
  ClockState3 clock_state3;
  clock_state3.b() = 0.1;     // [s]
  clock_state3.d() = 0.01;    // [s/s]
  clock_state3.dr() = 0.001;  // [s/s^2]
  Ptr<Clock> clock = std::make_shared<Clock>();
  clock->SetState(clock_state3);
  clock->SetModel(clock_model);
  Ptr<ClockDynamics> clock_dynamics = clock->GetDynamicsSharedPtr();

  // Dynamics
  Ptr<NBodyDynamics> dyn = std::make_shared<NBodyDynamics>();
  dyn->SetIntegrator(IntegratorType::RKF45);
  dyn->AddBody(Body::Moon(5, 5));
  dyn->SetTimeStep(1.0);  // [s]
  dyn->SetFrame(Frame::MOON_CI);
  dyn->SetAutodiff(true);

  // SRP
  double CR = 1.5;
  double area = 1.0;    // [m^2]
  double mass = 500.0;  // [kg]
  dyn->SetSrpCoefficient(CR, area, mass);

  // Parameter State
  ParamState param_state = dyn->GetParams();

  // Process Noise
  auto process_noise = [](const State& x, Real t0, Real tf) -> MatXd {
    int n = x.size();
    MatXd Q = MatXd::Zero(n, n);
    Real dt = tf - t0;
    // Simple constant process noise
    for (int i = 0; i < n; i++) {
      Q(i, i) = 1e-6 * dt;  // Adjust noise level as needed
    }
    return Q;
  };

  auto proc_noise_clock = [clock_model, clock_dynamics](const State& x, Real t0, Real tf) -> MatXd {
    return clock_dynamics->ThreeStateNoise(clock_model, tf - t0);
  };

  auto proc_noise_ptr = std::make_shared<ProcessNoiseFunction>(process_noise);
  auto proc_noise_clock_ptr = std::make_shared<ProcessNoiseFunction>(proc_noise_clock);

  // Define Joint State
  JointState joint_state;
  joint_state.Add(rv0_ci, dyn, std::move(proc_noise_ptr), param_state,
                  std::vector<EstType>{ESTIMATED, FIXED});
  joint_state.Add(clock->GetState(), clock_dynamics, std::move(proc_noise_clock_ptr), ParamState(0),
                  std::vector<EstType>{});
  int size = joint_state.GetStmSize();
  int state_size = joint_state.GetStateSize();
  int param_size = joint_state.GetParamSize();

  std::cout << "Joint State Summary:" << std::endl;
  std::cout << "Joint State STM   Size: " << size << std::endl;
  std::cout << "Joint State State Size: " << state_size << std::endl;
  std::cout << "Joint State Param Size: " << param_size << std::endl;
  std::cout << "  - Estimated Params  : " << joint_state.GetEstimatedParamSize() << std::endl;
  std::cout << "  - Considered Params : " << joint_state.GetConsideredParamSize() << std::endl;
  std::cout << "  - Fixed Params      : " << joint_state.GetFixedParamSize() << std::endl;
  std::cout << std::endl;

  State x0 = joint_state.GetStmState();

  std::cout << "Initial Joint STM State:" << std::endl;
  std::cout << x0 << std::endl;
  std::cout << std::endl;

  FilterDynamicsFunction dyn_func = joint_state.GetDynamicsFunction();
  MatXd stm = MatXd::Identity(size, size);
  State xf = dyn_func(x0, t0, t0 + 60.0, nullptr, &stm);

  std::cout << "Propagated Joint STM State after 60 seconds:" << std::endl;
  std::cout << xf << std::endl;
  std::cout << std::endl;

  std::cout << "State Transition Matrix (STM):" << std::endl;
  std::cout << stm << std::endl;
  return 0;
}

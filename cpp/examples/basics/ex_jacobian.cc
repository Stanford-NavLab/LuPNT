// Demonstrates state transition matrix generation for NBodyDynamics with
// autodiff, then compares it against a finite-difference STM.
#include <lupnt/lupnt.h>

#include <chrono>

using namespace lupnt;
using namespace std::chrono;

int main() {
  // Dynamics
  auto dyn = MakeUnique<NBodyDynamics>();
  dyn->SetIntegrator(IntegratorType::RKF45);
  dyn->SetIntegratorParams(IntegratorParams(20, 1e-12, 1e-12));
  dyn->SetFrame(Frame::MOON_CI);
  dyn->AddBody(Body::Moon(10, 10));
  dyn->AddBody(Body::Earth());
  dyn->SetTimeStep(10);
  dyn->SetAutodiff(true);

  // Initial state
  Real a = 6541.4e3;       // [m]
  Real e = 0.6;            // [-]
  Real i = 65.5 * RAD;     // [rad]
  Real Omega = 0.0 * RAD;  // [rad]
  Real w = 90.0 * RAD;     // [rad]
  Real M = 0.0 * RAD;      // [rad]

  // Time
  Real et0_utc = GregorianToTime(2025, 1, 1, 12, 0, 0).val();  // [s] UTC
  double et0 = UtcToTai(et0_utc).val();                        // [s] TAI
  double dt = 30 * 60;                                         // [s] Timestep

  // Initial state
  State coe = ClassicalOE({a, e, i, Omega, w, M}, Frame::MOON_CI);
  State rv0 = ClassicalToCart(coe, GM_MOON);

  MatXd Phi;

  // omp_set_num_threads(1); // Disable openmp threading
  std::cout << "Automatic differentiation" << std::endl;
  VecX rv = dyn->Propagate(rv0, et0, et0 + dt, nullptr, &Phi);

  // Compare with numerical Derivative
  MatXd Phi_num(rv0.size(), rv0.size());
  double eps = 1e-6;
  std::cout << "Numerical differentiation" << std::endl;
#pragma omp parallel for
  for (int i = 0; i < rv0.size(); i++) {
    Vec6 x_plus = rv0;
    Vec6 x_min = rv0;
    x_plus(i) += eps;
    x_min(i) -= eps;
    Vec6 dx_pert = dyn->Propagate(x_plus, et0, et0 + dt) - dyn->Propagate(x_min, et0, et0 + dt);
#pragma omp critical
    {
      Phi_num.col(i) = dx_pert.cast<double>() / (2 * eps);
    }
  }

  std::cout << "<Dynamics Example>" << std::endl;
  std::cout << " Initial State: " << rv0.transpose() << std::endl << std::endl;
  std::cout << " Propagated State: " << rv.transpose() << std::endl << std::endl;

  std::cout << " Jacobian: " << std::endl << Phi << std::endl << std::endl;
  std::cout << " Numerical Jacobian: " << std::endl << Phi_num << std::endl << std::endl;

  MatXd diff = (Phi - Phi_num).array() / Phi_num.array();
  std::cout << " Difference in STM (ratio): " << std::endl << diff << std::endl << std::endl;

  return 0;
}

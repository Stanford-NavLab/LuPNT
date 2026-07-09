// Example 1: Propagating an ELFO Orbit
// ------------------------------------------------------------------------------
// C++ counterpart of python/examples/ex1_propagate_orbit.ipynb.
//
// Initialise a lunar Elliptical Frozen Orbit (ELFO) from classical orbital
// elements (COE) defined in the Moon-centred Orbital Plane (OP) frame, convert
// the state to the Moon-Centred Inertial (MCI) frame, propagate it for 30 days
// with full N-body dynamics + Solar Radiation Pressure (SRP), and inspect the
// long-term libration of the osculating orbital elements.
//
// The OP frame tracks the instantaneous orbital plane of the Earth about the
// Moon, which is the natural frame in which the canonical frozen-orbit
// parameters (a, e, i, omega) of Ely (2005) are defined.

#include <iomanip>
#include <iostream>

#include "lupnt/lupnt.h"

using namespace lupnt;

int main() {
  // --- Epoch (LuPNT propagates in TAI seconds from J2000) --------------------
  const Real t0 = GregorianToTime(2009, 7, 15, 1, 0, 0);
  std::cout << "Epoch: " << TimeToGregorianString(t0) << "\n";

  // --- Initial state: ELFO in classical orbital elements (OP frame) ----------
  // COE layout: [a, e, i, RAAN, omega, M]; angles in radians, a in metres.
  const Real a = 6541.4e3;      // [m]   semi-major axis
  const Real e = 0.6;           // [-]   eccentricity
  const Real i = 56.2 * RAD;    // [rad] inclination
  const Real raan = 0.0 * RAD;  // [rad] right ascension of the ascending node
  const Real w = 90.0 * RAD;    // [rad] argument of perigee
  const Real M0 = 0.0 * RAD;    // [rad] mean anomaly
  Vec6 coe0_op(a, e, i, raan, w, M0);

  const Real period = GetOrbitalPeriod(a, GM_MOON);
  std::cout << "Orbital period: " << period / SECS_HOUR << " hours\n";

  // --- Convert COE -> Cartesian (OP), then OP -> MCI at t0 -------------------
  const Vec6 rv0_op = ClassicalToCart(coe0_op, GM_MOON);
  const Vec6 rv0_mci = ConvertFrame(t0, rv0_op, Frame::MOON_OP, Frame::MOON_CI);
  std::cout << "r0 (MCI): " << (rv0_mci.head(3) / 1e3).transpose() << " km\n";

  // --- Force model -----------------------------------------------------------
  //   Moon    : 2x2 spherical harmonics (J2 + C22)
  //   Earth   : point-mass third body (dominant long-period perturbation)
  //   Sun     : point-mass third body + SRP source
  //   SRP     : cannon-ball, C_R = 1.8, A = 0.02 m^2, m = 10 kg
  NBodyDynamics dyn;
  dyn.SetFrame(Frame::MOON_CI);
  dyn.SetUnits(SI_UNITS);
  dyn.AddBody(Body::Moon(2, 2));
  dyn.AddBody(Body::Earth());
  dyn.AddBody(Body::Sun());
  dyn.SetSrpCoefficient(1.8, 0.02, 10.0);

  dyn.SetIntegrator(IntegratorType::RKF45);
  IntegratorParams params;
  params.reltol = 1e-10;
  params.abstol = 1e-6;  // [m] 1 um position floor
  dyn.SetIntegratorParams(params);

  // --- Propagate for 30 days -------------------------------------------------
  const Real sim_duration = 30.0 * SECS_DAY;
  const Real dt_out = 30.0 * SECS_MINUTE;
  const VecX tspan = Arange(Real(0.0), sim_duration + dt_out, dt_out);
  const VecX tfs = t0 + tspan.array();
  std::cout << "Propagating " << tspan.size() << " epochs over " << sim_duration / SECS_DAY
            << " days ...\n";
  const MatX6 rv_mci = dyn.Propagate(rv0_mci, tfs);

  // --- Convert the propagated states to OP and extract osculating COE --------
  const MatX6 rv_op = ConvertFrame(tfs, rv_mci, Frame::MOON_CI, Frame::MOON_OP);
  const MatX6 coe = CartToClassical(rv_op, GM_MOON);

  const int last = static_cast<int>(coe.rows()) - 1;
  std::cout << std::fixed << std::setprecision(4);
  std::cout << "\nOsculating COE libration over the propagation (OP frame):\n";
  std::cout << "  a    [km]  : " << coe(0, 0) / 1e3 << " -> " << coe(last, 0) / 1e3 << "\n";
  std::cout << "  e    [-]   : " << coe(0, 1) << " -> " << coe(last, 1) << "\n";
  std::cout << "  i    [deg] : " << coe(0, 2) * DEG << " -> " << coe(last, 2) * DEG << "\n";
  std::cout << "  omega[deg] : " << coe(0, 4) * DEG << " -> " << coe(last, 4) * DEG << "\n";

  return 0;
}

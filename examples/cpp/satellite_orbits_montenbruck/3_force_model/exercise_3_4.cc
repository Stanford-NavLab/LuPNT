#include <lupnt/lupnt.h>

#include <iostream>
#include <string>

using namespace lupnt;
using namespace std;

// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
// applications. Berlin : New York: Springer, 2000.
// doi: 10.1007/978-3-642-58351-3.
int main() {
  // Time
  Real mjd0_utc = Gregorian2MJD(1999, 03, 01, 00, 00, 0.0);  // [days]
  Real t_tai0 = MJD2Time(mjd0_utc);                          // [s]

  // Real t_tai0 = ConvertTime(mjd0_utc, Time::MJD_UTC, Time::TAI);  //
  // [s]

  // Propagation
  Real dt = 120.0;                           // [s] Time step
  Real tf = 1.0 * SECS_DAY;                  // [s] Final time
  int N_step = (int)(tf / dt).val();         // [-] Number of steps
  VecX ts = VecX::LinSpaced(0, tf, N_step);  // [s] Time vector
  VecX t_tais = ts.array() + t_tai0;         // [s] TAI time vector

  // Initial state
  Vec6 coe0(7178, 1e-3, 98.57 * RAD, 0, 0, 0);  // [km, -, rad]
  Vec6 rv0 = Classical2Cart(coe0, GM_EARTH);    // [km, km/s]

  // Spacecraft
  Real area = 5.0;     // [m^2] Cross-sectional area
  Real mass = 1000.0;  // [kg] Spacecraft mass
  Real CR = 1.3;       // [-] Radiation pressure coefficient
  Real CD = 2.3;       // [-] Drag coefficient

  // Gravity field
  std::string grav_file = "JGM3.cof";
  int n_max = 20, m_max = 20;
  bool normalized = true;
  GravityField grav = ReadHarmonicGravityField<double>(grav_file, n_max, m_max, normalized);

  Body earth = BodyT<double>::Earth();
  earth.use_gravity_field = true;
  earth.gravity_field = grav;

  Body moon = BodyT<double>::Moon();
  Body sun = BodyT<double>::Sun();

  // Dynamics
  NBodyDynamics dyn;
  dyn.SetArea(area);
  dyn.SetMass(mass);
  dyn.SetSrpCoeff(CR);
  dyn.SetDragCoeff(CD);
  dyn.SetFrame(Frame::ITRF);
  dyn.AddBody(sun);
  dyn.AddBody(moon);

  Vec6 rv_dot = dyn.ComputeRates(t_tai0, rv0);
  return 0;
}

#include "lupnt/physics/orbit_state/anomaly.h"

#include "lupnt/numerics/math_utils.h"

namespace lupnt {
  Real GetOrbitalPeriod(Real a, Real GM) { return 2 * PI * sqrt(pow(a, 3) / GM); }
  VEC_IMP_REAL_REAL(GetOrbitalPeriod);

  Real Ecc2TrueAnomaly(Real E, Real e) { return atan2(sqrt(1 - pow(e, 2)) * sin(E), cos(E) - e); }

  Real Ecc2MeanAnomaly(Real E, Real e) { return Wrap2Pi(E - e * sin(E)); }

  Real Mean2EccAnomaly(Real M, Real e) {
    Real MM = Wrap2Pi(M);

    // Initial estimate of E
    Real E = MM;
    Real Eest = E - (E - e * sin(E) - MM) / (1.0 - e * cos(E));

    int max_itr = 100;
    int itr = 0;
    while ((abs(Eest - E) >= EPS) && (itr <= max_itr)) {
      E = Eest;
      Eest = E - (E - e * sin(E) - M) / (1.0 - e * cos(E));
      itr++;
    }
    E = Eest;

    return Wrap2Pi(E);
  }

  Real Mean2TrueAnomaly(Real M, Real e) {
    Real E = Mean2EccAnomaly(M, e);
    Real nu = Ecc2TrueAnomaly(E, e);
    return Wrap2Pi(nu);
  }

  Real True2EccAnomaly(Real nu, Real e) {
    Real E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(nu / 2));
    return E;
  }

  Real True2MeanAnomaly(Real nu, Real e) {
    Real E = True2EccAnomaly(nu, e);
    Real M = Ecc2MeanAnomaly(E, e);
    return M;
  }

  VEC_IMP_REAL_REAL(Ecc2TrueAnomaly);
  VEC_IMP_REAL_REAL(Ecc2MeanAnomaly);
  VEC_IMP_REAL_REAL(Mean2EccAnomaly);
  VEC_IMP_REAL_REAL(True2EccAnomaly);
  VEC_IMP_REAL_REAL(Mean2TrueAnomaly);
  VEC_IMP_REAL_REAL(True2MeanAnomaly);
}  // namespace lupnt

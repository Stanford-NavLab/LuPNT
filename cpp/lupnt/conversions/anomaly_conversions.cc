#include "lupnt/conversions/anomaly_conversions.h"

#include "lupnt/numerics/math_utils.h"

namespace lupnt {
  Real GetOrbitalPeriod(Real a, Real GM) { return 2 * PI * sqrt(pow(a, 3) / GM); }
  VEC_IMP_REAL_REAL(GetOrbitalPeriod);

  Real EccToTrueAnomaly(Real E, Real e) { return atan2(sqrt(1 - pow(e, 2)) * sin(E), cos(E) - e); }

  Real EccToMeanAnomaly(Real E, Real e) { return WrapToPi(E - e * sin(E)); }

  Real MeanToEccAnomaly(Real M, Real e) {
    Real MM = WrapToPi(M);

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

    return WrapToPi(E);
  }

  Real MeanToTrueAnomaly(Real M, Real e) {
    Real E = MeanToEccAnomaly(M, e);
    Real nu = EccToTrueAnomaly(E, e);
    return WrapToPi(nu);
  }

  Real TrueToEccAnomaly(Real nu, Real e) {
    Real E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(nu / 2));
    return E;
  }

  Real TrueToMeanAnomaly(Real nu, Real e) {
    Real E = TrueToEccAnomaly(nu, e);
    Real M = EccToMeanAnomaly(E, e);
    return M;
  }

  VEC_IMP_REAL_REAL(EccToTrueAnomaly);
  VEC_IMP_REAL_REAL(EccToMeanAnomaly);
  VEC_IMP_REAL_REAL(MeanToEccAnomaly);
  VEC_IMP_REAL_REAL(TrueToEccAnomaly);
  VEC_IMP_REAL_REAL(MeanToTrueAnomaly);
  VEC_IMP_REAL_REAL(TrueToMeanAnomaly);
}  // namespace lupnt

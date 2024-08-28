#pragma once

#include "lupnt/numerics/vector_macros.h"
#include "orbit_states.h"

namespace lupnt {

  Real GetOrbitalPeriod(Real a, Real GM);
  VEC_DEF_REAL_REAL(GetOrbitalPeriod);

  // Anomaly
  Real Ecc2TrueAnomaly(Real E, Real e);
  Real Ecc2MeanAnomaly(Real E, Real e);
  Real Mean2EccAnomaly(Real M, Real e);
  Real True2EccAnomaly(Real nu, Real e);
  Real Mean2TrueAnomaly(Real M, Real e);
  Real True2MeanAnomaly(Real f, Real e);

  // Vector definitions
  VEC_DEF_REAL_REAL(GetPeriod);
  VEC_DEF_REAL_REAL(Ecc2TrueAnomaly);
  VEC_DEF_REAL_REAL(Ecc2MeanAnomaly);
  VEC_DEF_REAL_REAL(Mean2EccAnomaly);
  VEC_DEF_REAL_REAL(True2EccAnomaly);
  VEC_DEF_REAL_REAL(Mean2TrueAnomaly);
  VEC_DEF_REAL_REAL(True2MeanAnomaly);

}  // namespace lupnt

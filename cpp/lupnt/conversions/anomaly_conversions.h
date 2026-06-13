#pragma once

#include "lupnt/numerics/vector_macros.h"
#include "lupnt/states/state.h"

namespace lupnt {

  Real GetOrbitalPeriod(Real a, Real GM);
  VEC_DEF_REAL_REAL(GetOrbitalPeriod);

  // Anomaly
  Real EccToTrueAnomaly(Real E, Real e);
  Real EccToMeanAnomaly(Real E, Real e);
  Real MeanToEccAnomaly(Real M, Real e);
  Real TrueToEccAnomaly(Real nu, Real e);
  Real MeanToTrueAnomaly(Real M, Real e);
  Real TrueToMeanAnomaly(Real f, Real e);

  // Vector definitions
  VEC_DEF_REAL_REAL(GetOrbitalPeriod);
  VEC_DEF_REAL_REAL(EccToTrueAnomaly);
  VEC_DEF_REAL_REAL(EccToMeanAnomaly);
  VEC_DEF_REAL_REAL(MeanToEccAnomaly);
  VEC_DEF_REAL_REAL(TrueToEccAnomaly);
  VEC_DEF_REAL_REAL(MeanToTrueAnomaly);
  VEC_DEF_REAL_REAL(TrueToMeanAnomaly);

}  // namespace lupnt

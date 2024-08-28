#pragma once

#include <string>
#include <tuple>

#include "lupnt/core/constants.h"
#include "lupnt/numerics/graphs.h"
#include "lupnt/numerics/vector_macros.h"

namespace lupnt {

  Real ConvertTime(Real t, Time from, Time to);
  VecX ConvertTime(VecX t, Time from, Time to);

  Real UTC2UT1(Real t_utc);
  Real UT12UTC(Real t_ut1);

  Real TAI2UTC(Real t_tai);
  Real UTC2TAI(Real t_utc);

  Real TAI2TT(Real t_tai);
  Real TT2TAI(Real t_tt);

  Real TCG2TT(Real t_tcg);
  Real TT2TCG(Real t_tt);

  Real TT2TDB(Real t_tt);
  Real TDB2TT(Real t_tdb);

  Real TAI2GPS(Real t_tai);
  Real GPS2TAI(Real t_gps);

  Real TCB2TDB(Real t_tcb);
  Real TT2TCB(Real t_tdb);

  Real MJD2Time(Real mjd);
  Real Time2MJD(Real t);

  Real JD2Time(Real jd);
  Real Time2JD(Real t);

  Real EarthRotationAngle(Real t_ut1);
  Real Gregorian2MJD(int year, int month, int day, int hour = 0, int min = 0, Real sec = 0);
  Real Gregorian2Time(int year, int month, int day, int hour = 0, int min = 0, Real sec = 0);

  Real GreenwichMeanSiderealTime(Real mjd_ut1);
  Real GreenwichApparentSiderealTime(Real mjd_ut1);

  std::tuple<int, int, int, int, int, Real> MJD2Gregorian(Real mjd);

  std::string MJD2GregorianString(Real mjd, int precision = 3);
  std::string Time2GregorianString(Real t, int precision = 3);

  VEC_DEF_REAL(UTC2UT1)
  VEC_DEF_REAL(UT12UTC)
  VEC_DEF_REAL(TAI2UTC)
  VEC_DEF_REAL(UTC2TAI)
  VEC_DEF_REAL(TAI2TT)
  VEC_DEF_REAL(TT2TAI)
  VEC_DEF_REAL(TCG2TT)
  VEC_DEF_REAL(TT2TCG)
  VEC_DEF_REAL(TT2TDB)
  VEC_DEF_REAL(TDB2TT)
  VEC_DEF_REAL(TAI2GPS)
  VEC_DEF_REAL(GPS2TAI)
  VEC_DEF_REAL(TCB2TDB)
  VEC_DEF_REAL(TT2TCB)

  VEC_DEF_REAL(MJD2Time)
  VEC_DEF_REAL(Time2MJD)
  VEC_DEF_REAL(JD2Time)
  VEC_DEF_REAL(Time2JD)

  VEC_DEF_REAL(EarthRotationAngle)
  VEC_DEF_REAL(GreenwichMeanSiderealTime)
  VEC_DEF_REAL(GreenwichApparentSiderealTime)

}  // namespace lupnt

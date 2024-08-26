/**
 * @file FrameConverter.cpp
 * @author Stanford NAV LAB
 * @brief Coordinate conversion functions
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <functional>
#include <map>
#include <memory>
#include <ostream>

#include "cheby.h"
#include "lupnt/core/constants.h"

namespace lupnt {

  Vec6 ITRF2GCRF(Real t_tai, const Vec6 &rv_itrf);
  Vec6 GCRF2ITRF(Real t_tai, const Vec6 &rv_gcrf);

  Vec6 GCRF2EME(Real t_tai, const Vec6 &rv_gcrf);
  Vec6 EME2GCRF(Real t_tai, const Vec6 &rv_eme);

  Vec6 GCRF2ICRF(Real t_tai, const Vec6 &rv_gcrf);
  Vec6 ICRF2GCRF(Real t_tai, const Vec6 &rv_icrf);

  Vec6 GCRF2MoonCI(Real t_tai, const Vec6 &rv_gcrf);
  Vec6 MoonCI2GCRF(Real t_tai, const Vec6 &rv_mi);

  Vec6 MoonCI2MoonPA(Real t_tai, const Vec6 &rv_mi);
  Vec6 MoonPA2MoonCI(Real t_tai, const Vec6 &rv_pa);

  Vec6 MoonPA2MoonME(Real t_tai, const Vec6 &rv_pa);
  Vec6 MoonME2MoonPA(Real t_tai, const Vec6 &rv_me);

  Vec6 GCRF2EMR(Real t_tai, const Vec6 &rv_gcrf);
  Vec6 EMR2GCRF(Real t_tai, const Vec6 &rv_emr);

  Vec6 MoonCI2MoonOP(Real t_tai, const Vec6 &rv_ci);
  Vec6 MoonOP2MoonCI(Real t_tai, const Vec6 &rv_op);

  // Rotations
  Mat3 RotPrecessionNutation(Real t_tai);
  Mat3 RotSideralMotion(Real t_tai);
  Mat3 RotSideralMotionDot(Real t_tai);
  Mat3 RotPolarMotion(Real t_tai);

  Mat3d RotGCRF2EME();
  Mat3d RotGCRF2EMEFirstOrder();
  Mat3d RotGCRF2EMESecondOrder();

  std::pair<Mat3, Mat3> RotMoonCI2MoonPA(Real t_tai);
  Mat3d RotMoonPA2MoonME();
  Mat3 RotOP2CI(Real t_tai);

}  // namespace lupnt

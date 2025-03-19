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

#include "lupnt/physics/frame_conversions.h"

#include <filesystem>

#include "lupnt/core/constants.h"
#include "lupnt/data/eop.h"
#include "lupnt/data/iau_sofa.h"
#include "lupnt/data/kernels.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/cheby.h"
#include "lupnt/physics/frame_converter_spice.h"
#include "lupnt/physics/orbit_state.h"
#include "lupnt/physics/time_converter.h"

#define FRAME_CONVERSION(from, to, func) \
  {{Frame::from, Frame::to}, [](Real t, const Vec6& rv) -> Vec6 { return func(t, rv); }}

namespace lupnt {

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 34
  Mat3 RotPrecessionNutation(Real t_tai) {
    Real t_tt = ConvertTime(t_tai, Time::TAI, Time::TT);
    Real jd_tt = Time2JD(t_tt);

    IauSofaData iau_data = GetIauSofaData(jd_tt);
    Real X = iau_data.X * RAD_ARCSEC;
    Real Y = iau_data.Y * RAD_ARCSEC;
    Real s = iau_data.s * RAD_ARCSEC;

    Real a = 1. / (1. + sqrt(1. - X * X - Y * Y));
    Mat3 mat{{1. - a * X * X, -a * X * Y, -X},
             {-a * X * Y, 1. - a * Y * Y, -Y},
             {X, Y, 1. - a * (X * X + Y * Y)}};
    Mat3 R_pn = RotZ(-s) * mat;
    return R_pn;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 35
  Mat3 RotSideralMotion(Real t_tai) {
    Real t_ut1 = ConvertTime(t_tai, Time::TAI, Time::UT1);
    Real theta_era = EarthRotationAngle(t_ut1);
    Mat3 R_s = RotZ(theta_era);
    return R_s;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 36
  Mat3 RotSideralMotionDot(Real t_tai) {
    Real t_ut1 = ConvertTime(t_tai, Time::TAI, Time::UT1);
    Real theta_era = EarthRotationAngle(t_ut1);

    Real t_utc = ConvertTime(t_tai, Time::TAI, Time::UTC);
    Real mjd_utc = Time2MJD(t_utc);
    EopData eop = GetEopData(mjd_utc);
    Real lod = eop.lod;

    Real w_E = 7.292115146706979e-5 * (1. - lod / SECS_DAY);
    Mat3 R_s_dot{{-w_E * sin(theta_era), w_E * cos(theta_era), 0},
                 {-w_E * cos(theta_era), -w_E * sin(theta_era), 0},
                 {0, 0, 0}};
    return R_s_dot;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 36
  Mat3 RotPolarMotion(Real t_tai) {
    Real t_utc = ConvertTime(t_tai, Time::TAI, Time::UTC);
    Real t_tt = ConvertTime(t_tai, Time::TAI, Time::TT);
    Real mjd_utc = Time2MJD(t_utc);

    EopData eop = GetEopData(mjd_utc);
    Real xp = eop.x_pole;
    Real yp = eop.y_pole;
    Real sp = -47e-6 * RAD_ARCSEC * (t_tt / DAYS_CENTURY);

    Mat3 R_po = RotX(-yp) * RotY(-xp) * RotZ(sp);
    return R_po;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 37
  Vec6 GCRF2ITRF(Real t_tai, const Vec6& rv_gcrf) {
    Mat3 R_po = RotPolarMotion(t_tai);
    Mat3 R_pn = RotPrecessionNutation(t_tai);
    Mat3 R_s = RotSideralMotion(t_tai);
    Mat3 R_s_dot = RotSideralMotionDot(t_tai);

    Mat3 R_gcrf2itrf = R_po * R_s * R_pn;
    Mat3 R_gcrf2itrf_dot = R_po * R_s_dot * R_pn;

    Vec3 r_gcrf = rv_gcrf.head(3);
    Vec3 v_gcrf = rv_gcrf.tail(3);

    Vec3 r_itrf = R_gcrf2itrf * r_gcrf;
    Vec3 v_itrf = R_gcrf2itrf * v_gcrf + R_gcrf2itrf_dot * r_gcrf;

    Vec6 rv_itrf;
    rv_itrf << r_itrf, v_itrf;
    return rv_itrf;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 37
  Vec6 ITRF2GCRF(Real t_tai, const Vec6& rv_itrf) {
    Mat3 R_po = RotPolarMotion(t_tai);
    Mat3 R_pn = RotPrecessionNutation(t_tai);
    Mat3 R_s = RotSideralMotion(t_tai);
    Mat3 R_s_dot = RotSideralMotionDot(t_tai);

    Mat3 R_gcrf2itrf = R_po * R_s * R_pn;
    Mat3 R_gcrf2itrf_dot = R_po * R_s_dot * R_pn;

    Vec3 r_itrf = rv_itrf.head(3);
    Vec3 v_itrf = rv_itrf.tail(3);

    Vec3 r_gcrf = R_gcrf2itrf.transpose() * r_itrf;
    Vec3 v_gcrf = R_gcrf2itrf.transpose() * (v_itrf - R_gcrf2itrf_dot * r_gcrf);

    Vec6 rv_gcrf;
    rv_gcrf << r_gcrf, v_gcrf;
    return rv_gcrf;
  }

  /// @note
  /// Also known as Earth Frame Bias Matrix
  /// Astrodynamics Convention & Modeling Reference, Version 1.1, Page 38
  Mat3d RotGCRF2EME() {
    const double da = FRAME_BIAS_DALPHA0;
    const double xi = FRAME_BIAS_XI0;
    const double eta = FRAME_BIAS_ETA0;
    Mat3d B_e = RotX(-eta) * RotY(xi) * RotZ(da);
    return B_e;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 38
  Mat3d RotGCRF2EMEFirstOrder() {
    const double da = FRAME_BIAS_DALPHA0;
    const double xi = FRAME_BIAS_XI0;
    const double eta = FRAME_BIAS_ETA0;
    Mat3d B_e{{1, da, -xi}, {-da, 1, -eta}, {xi, eta, 1}};
    return B_e;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 38
  Mat3d RotGCRF2EMESecondOrder() {
    const double da = FRAME_BIAS_DALPHA0;
    const double xi = FRAME_BIAS_XI0;
    const double eta = FRAME_BIAS_ETA0;
    Mat3d B_e{{1 - 0.5 * (da * da + xi * xi), da, -xi},
              {-da - eta * xi, 1 - 0.5 * (da * da + eta * eta), -eta},
              {xi - eta * da, eta + xi * da, 1 - 0.5 * (eta * eta + xi * xi)}};
    return B_e;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 39
  Vec6 GCRF2EME(Real t_tai, const Vec6& rv_gcrf) {
    (void)t_tai;
    Mat3d B_e = RotGCRF2EME();
    Vec3 r_gcrf = rv_gcrf.head(3);
    Vec3 v_gcrf = rv_gcrf.tail(3);

    Vec3 r_eme = B_e * r_gcrf;
    Vec3 v_eme = B_e * v_gcrf;

    Vec6 rv_eme;
    rv_eme << r_eme, v_eme;
    return rv_eme;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 39
  Vec6 EME2GCRF(Real t_tai, const Vec6& rv_eme) {
    (void)t_tai;
    Mat3d B_e = RotGCRF2EME();
    Vec3 r_eme = rv_eme.head(3);
    Vec3 v_eme = rv_eme.tail(3);

    Mat3d B_e_inv = B_e.transpose();
    Vec3 r_gcrf = B_e_inv * r_eme;
    Vec3 v_gcrf = B_e_inv * v_eme;

    Vec6 rv_gcrf;
    rv_gcrf << r_gcrf, v_gcrf;
    return rv_gcrf;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 39
  Vec6 GCRF2ICRF(Real t_tai, const Vec6& rv_gcrf) {
    Vec6 rv_ssb2earth = GetBodyPosVel(t_tai, NaifId::SSB, NaifId::EARTH, Frame::GCRF);
    Vec6 rv_icrf = rv_gcrf + rv_ssb2earth;
    return rv_icrf;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 40
  Vec6 ICRF2GCRF(Real t_tai, const Vec6& rv_icrf) {
    Vec6 rv_ssb2earth = GetBodyPosVel(t_tai, NaifId::SSB, NaifId::EARTH, Frame::GCRF);
    Vec6 rv_gcrf = rv_icrf - rv_ssb2earth;
    return rv_gcrf;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 40
  Vec6 GCRF2MoonCI(Real t_tai, const Vec6& rv_gcrf) {
    Vec6 rv_earth2moon = GetBodyPosVel(t_tai, NaifId::EARTH, NaifId::MOON, Frame::GCRF);
    Vec6 rv_mi = rv_gcrf - rv_earth2moon;
    return rv_mi;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 40
  Vec6 MoonCI2GCRF(Real t_tai, const Vec6& rv_mi) {
    Vec6 rv_earth2moon = GetBodyPosVel(t_tai, NaifId::EARTH, NaifId::MOON, Frame::GCRF);
    Vec6 rv_gcrf = rv_mi + rv_earth2moon;
    return rv_gcrf;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 42
  std::pair<Mat3, Mat3> RotMoonCI2MoonPA(Real t_tai) {
    (void)t_tai;
    Vec6 lunar_mantle = GetLunarMantleData(t_tai);
    auto [phi, theta, psi, phi_dot, theta_dot, psi_dot] = unpack(lunar_mantle);

    Real spsi = sin(psi);
    Real cpsi = cos(psi);
    Mat3 mat{
        {-psi_dot * spsi, psi_dot * cpsi, 0}, {-psi_dot * cpsi, -psi_dot * spsi, 0}, {0, 0, 0}};
    Mat3 R_mi2pa = RotZ(psi) * RotX(theta) * RotZ(phi);
    Mat3 R_mi2pa_dot = mat * RotX(theta) * RotZ(phi);
    return {R_mi2pa, R_mi2pa_dot};
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 42
  Vec6 MoonCI2MoonPA(Real t_tai, const Vec6& rv_mi) {
    auto [R_mi2pa, R_mi2pa_dot] = RotMoonCI2MoonPA(t_tai);

    Vec3 r_mi = rv_mi.head(3);
    Vec3 v_mi = rv_mi.tail(3);

    Vec3 r_pa = R_mi2pa * r_mi;
    Vec3 v_pa = R_mi2pa * v_mi + R_mi2pa_dot * r_mi;

    Vec6 rv_pa;
    rv_pa << r_pa, v_pa;
    return rv_pa;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 42
  Vec6 MoonPA2MoonCI(Real t_tai, const Vec6& rv_pa) {
    auto [R_mi2pa, R_mi2pa_dot] = RotMoonCI2MoonPA(t_tai);

    Vec3 r_pa = rv_pa.head(3);
    Vec3 v_pa = rv_pa.tail(3);

    Vec3 r_mi = R_mi2pa.transpose() * r_pa;
    Vec3 v_mi = R_mi2pa.transpose() * (v_pa - R_mi2pa_dot * r_mi);
    Vec6 rv_mi;
    rv_mi << r_mi, v_mi;
    return rv_mi;
  }

  Mat3d RotMoonPA2MoonME() {
    Mat3d B_moon
        = RotX(-0.2785 * RAD_ARCSEC) * RotY(-78.6944 * RAD_ARCSEC) * RotZ(-67.8526 * RAD_ARCSEC);
    return B_moon;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 43
  Vec6 MoonPA2MoonME(Real t_tai, const Vec6& rv_pa) {
    (void)t_tai;
    Mat3d B_moon = RotMoonPA2MoonME();
    Vec3 r_pa = rv_pa.head(3);
    Vec3 v_pa = rv_pa.tail(3);

    Vec3 r_me = B_moon * r_pa;
    Vec3 v_me = B_moon * v_pa;

    Vec6 rv_me;
    rv_me << r_me, v_me;
    return rv_me;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 43
  Vec6 MoonME2MoonPA(Real t_tai, const Vec6& rv_me) {
    (void)t_tai;
    Mat3d B_moon = RotMoonPA2MoonME();
    Vec3 r_me = rv_me.head(3);
    Vec3 v_me = rv_me.tail(3);

    Mat3d B_moon_inv = B_moon.transpose();
    Vec3 r_pa = B_moon_inv * r_me;
    Vec3 v_pa = B_moon_inv * v_me;

    Vec6 rv_pa;
    rv_pa << r_pa, v_pa;
    return rv_pa;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 48
  Vec6 GCRF2EMR(Real t_tai, const Vec6& rv_gcrf) {
    Vec6 rv_earth2emb = GetBodyPosVel(t_tai, NaifId::EARTH, NaifId::EMB, Frame::GCRF);
    Vec6 rv_emr = Inertial2Synodic(rv_earth2emb, rv_gcrf);
    return rv_emr;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 48
  Vec6 EMR2GCRF(Real t_tai, const Vec6& rv_emr) {
    Vec6 rv_earth2emb = GetBodyPosVel(t_tai, NaifId::EARTH, NaifId::EMB, Frame::GCRF);
    Vec6 rv_gcrf = Synodic2Intertial(rv_earth2emb, rv_emr);
    return rv_gcrf;
  }

  /// @note
  /// T. A. Ely, ‘Stable Constellations of Frozen Elliptical Inclined Lunar Orbits’, J of
  /// Astronaut Sci, vol. 53, no. 3, pp. 301–316, Sep. 2005, doi: 10.1007/BF03546355.
  Mat3 RotOP2CI(Real t_tai) {
    Vec6 rv_m2e = GetBodyPosVel(t_tai, NaifId::MOON, NaifId::EARTH, Frame::GCRF);
    Vec3 r = rv_m2e.head(3);
    Vec3 v = rv_m2e.tail(3);

    Vec3 z_op = r.cross(v).normalized();
    Vec3 i_pole_me = Vec3::UnitZ();
    Vec3 i_pole = ConvertFrame(t_tai, i_pole_me, Frame::MOON_ME, Frame::MOON_CI);
    Vec3 x_op = i_pole.cross(z_op).normalized();
    Vec3 y_op = z_op.cross(x_op).normalized();

    // R = [x_OP, y_OP, z_OP]_CI (expressed in CI)
    Mat3 R_op2ci;
    R_op2ci << x_op, y_op, z_op;
    return R_op2ci;
  }

  /// @brief
  /// @param t_tai
  /// @param rv_ci
  /// @return Vec6
  /// @note
  /// T. A. Ely, ‘Stable Constellations of Frozen Elliptical Inclined Lunar Orbits’, J of
  /// Astronaut Sci, vol. 53, no. 3, pp. 301–316, Sep. 2005, doi: 10.1007/BF03546355.
  Vec6 MoonCI2MoonOP(Real t_tai, const Vec6& rv_ci) {
    Mat3 R_ci2op = RotOP2CI(t_tai).transpose();
    Vec3 r_op = R_ci2op * rv_ci.head(3);
    Vec3 v_op = R_ci2op * rv_ci.tail(3);
    Vec6 rv_op;
    rv_op << r_op, v_op;
    return rv_op;
  }

  /// @brief
  /// @param t_tai
  /// @param rv_op
  /// @return Vec6
  /// @note
  /// T. A. Ely, ‘Stable Constellations of Frozen Elliptical Inclined Lunar Orbits’, J of
  /// Astronaut Sci, vol. 53, no. 3, pp. 301–316, Sep. 2005, doi: 10.1007/BF03546355.
  Vec6 MoonOP2MoonCI(Real t_tai, const Vec6& rv_op) {
    Mat3 R_op2ci = RotOP2CI(t_tai);
    Vec3 r_ci = R_op2ci * rv_op.head(3);
    Vec3 v_ci = R_op2ci * rv_op.tail(3);
    Vec6 rv_ci;
    rv_ci << r_ci, v_ci;
    return rv_ci;
  }

}  // namespace lupnt

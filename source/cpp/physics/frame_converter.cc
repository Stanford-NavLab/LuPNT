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

#include "lupnt/physics/frame_converter.h"

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

  std::map<std::pair<Frame, Frame>, std::function<Vec6(Real, const Vec6& rv)>> frame_conversions = {
      FRAME_CONVERSION(GCRF, ITRF, GCRF2ITRF),
      FRAME_CONVERSION(ITRF, GCRF, ITRF2GCRF),
      FRAME_CONVERSION(GCRF, EME, GCRF2EME),
      FRAME_CONVERSION(EME, GCRF, EME2GCRF),
      FRAME_CONVERSION(GCRF, ICRF, GCRF2ICRF),
      FRAME_CONVERSION(ICRF, GCRF, ICRF2GCRF),
      FRAME_CONVERSION(GCRF, MOON_CI, GCRF2MoonCI),
      FRAME_CONVERSION(MOON_CI, GCRF, MoonMI2GCRF),
      FRAME_CONVERSION(MOON_CI, MOON_PA, MoonMI2MoonPA),
      FRAME_CONVERSION(MOON_PA, MOON_CI, MoonPA2MoonCI),
      FRAME_CONVERSION(MOON_PA, MOON_ME, MoonPA2MoonME),
      FRAME_CONVERSION(MOON_ME, MOON_PA, MoonME2MoonPA),
      FRAME_CONVERSION(GCRF, EMR, GCRF2EMR),
      FRAME_CONVERSION(EMR, GCRF, EMR2GCRF),
      FRAME_CONVERSION(MOON_ME, MOON_OP, MoonME2MoonOP),
      FRAME_CONVERSION(MOON_OP, MOON_ME, MoonOP2MoonME),
  };

  /// @brief Convert the state vector from one coordinate system to another
  /// (with integer ID input)
  /// @param t_tai Epoch of the state vector
  /// @param rv_in   State vector in the original coordinate system
  /// @param frame_in Coordinate system of the original state vector
  /// @param frame_out Coordinate system of the converted state vector
  /// @return Vec6  State vector in the converted coordinate system
  ///
  Vec6 ConvertFrame(Real t_tai, const Vec6& rv_in, Frame frame_in, Frame frame_out) {
    if (frame_in == frame_out) return rv_in;
    std::vector<Frame> path = FindShortestPath(frame_in, frame_out, frame_conversions);
    Vec6 rv_out = rv_in;
    for (size_t i = 0; i < path.size() - 1; i++) {
      std::function<Vec6(Real, const Vec6& rv)> f = frame_conversions[{path[i], path[i + 1]}];
      rv_out = f(t_tai, rv_out);
    }
    return rv_out;
  }

  CartesianOrbitState ConvertFrame(Real t_tai, const CartesianOrbitState& state_in,
                                   Frame frame_out) {
    Vec6 rv_in = state_in.GetVec();
    Vec6 rv_out = ConvertFrame(t_tai, rv_in, state_in.GetCoordSystem(), frame_out);
    return CartesianOrbitState(rv_out, frame_out);
  }

  Vec3 ConvertFrame(Real t_tai, const Vec3& r_in, Frame frame_in, Frame frame_out) {
    Vec6 rv_in;
    rv_in << r_in, Vec3::Zero();
    Vec6 rv_out_6 = ConvertFrame(t_tai, rv_in, frame_in, frame_out);
    return rv_out_6.head(3);
  }

  Mat<-1, 6> ConvertFrame(Real t_tai, const Mat<-1, 6>& rv_in, Frame frame_in, Frame frame_out) {
    Mat<-1, 6> rv_out(rv_in.rows(), 6);
    for (int i = 0; i < rv_in.rows(); i++) {
      rv_out.row(i) = ConvertFrame(t_tai, rv_in.row(i).transpose().eval(), frame_in, frame_out);
    }
    return rv_out;
  }

  Mat<-1, 3> ConvertFrame(Real t_tai, const Mat<-1, 3>& r_in, Frame frame_in, Frame frame_out) {
    Mat<-1, 3> r_out(r_in.rows(), 3);
    for (int i = 0; i < r_in.rows(); i++) {
      Vec6 rv_in;
      rv_in << r_in.row(i).transpose(), Vec3::Zero();
      r_out.row(i) = ConvertFrame(t_tai, rv_in, frame_in, frame_out).head(3);
    }
    return r_out;
  }

  Mat<-1, 6> ConvertFrame(VecX t_tai, const Vec6& rv_in, Frame frame_in, Frame frame_out) {
    Mat<-1, 6> rv_out(t_tai.size(), 6);
    for (int i = 0; i < t_tai.size(); i++) {
      rv_out.row(i) = ConvertFrame(t_tai(i), rv_in, frame_in, frame_out);
    }
    return rv_out;
  }

  Mat<-1, 3> ConvertFrame(VecX t_tai, const Vec3& r_in, Frame frame_in, Frame frame_out) {
    Mat<-1, 3> r_out(t_tai.size(), 3);
    Vec6 rv_in;
    rv_in << r_in, Vec3::Zero();
    for (int i = 0; i < t_tai.size(); i++)
      r_out.row(i) = ConvertFrame(t_tai(i), rv_in, frame_in, frame_out).head(3);
    return r_out;
  }

  Mat<-1, 6> ConvertFrame(VecX t_tai, const Mat<-1, 6>& rv_in, Frame frame_in, Frame frame_out) {
    assert(t_tai.size() == rv_in.rows() && "Epoch and rv_in must have same size");
    Mat<-1, 6> rv_out(t_tai.size(), 6);
    for (int i = 0; i < t_tai.size(); i++) {
      rv_out.row(i) = ConvertFrame(t_tai(i), rv_in.row(i).transpose().eval(), frame_in, frame_out);
    }
    return rv_out;
  }

  Mat<-1, 3> ConvertFrame(VecX t_tai, const Mat<-1, 3>& r_in, Frame frame_in, Frame frame_out) {
    assert(t_tai.size() == r_in.rows() && "Epoch and r_in must have same size");
    Mat<-1, 6> rv_in(t_tai.size(), 6);
    rv_in << r_in, Mat<-1, 3>::Zero(t_tai.size(), 3);
    Mat<-1, 6> rv_out = ConvertFrame(t_tai, rv_in, frame_in, frame_out);
    return rv_out.leftCols(3);
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 34
  Mat3 RotPrecessionNutation(Real t_tai) {
    Real t_tt = ConvertTime(t_tai, TimeSys::TAI, TimeSys::TT);
    Real jd_tt = TimeToJD(t_tt);

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

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 35
  Mat3 RotSideralMotion(Real t_tai) {
    Real t_ut1 = ConvertTime(t_tai, TimeSys::TAI, TimeSys::UT1);
    Real theta_era = EarthRotationAngle(t_ut1);
    Mat3 R_s = RotZ(theta_era);
    return R_s;
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 36
  Mat3 RotSideralMotionDot(Real t_tai) {
    Real t_ut1 = ConvertTime(t_tai, TimeSys::TAI, TimeSys::UT1);
    Real theta_era = EarthRotationAngle(t_ut1);

    Real t_utc = ConvertTime(t_tai, TimeSys::TAI, TimeSys::UTC);
    Real mjd_utc = TimeToMJD(t_utc);
    EopData eop = GetEopData(mjd_utc);
    Real lod = eop.lod;

    Real w_E = 7.292115146706979e-5 * (1. - lod / SECS_DAY);
    Mat3 R_s_dot{{-w_E * sin(theta_era), w_E * cos(theta_era), 0},
                 {-w_E * cos(theta_era), -w_E * sin(theta_era), 0},
                 {0, 0, 0}};
    return R_s_dot;
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 36
  Mat3 RotPolarMotion(Real t_tai) {
    Real t_utc = ConvertTime(t_tai, TimeSys::TAI, TimeSys::UTC);
    Real t_tt = ConvertTime(t_tai, TimeSys::TAI, TimeSys::TT);
    Real mjd_utc = TimeToMJD(t_utc);

    EopData eop = GetEopData(mjd_utc);
    Real xp = eop.x_pole;
    Real yp = eop.y_pole;
    Real sp = -47e-6 * RAD_ARCSEC * (t_tt / DAYS_CENTURY);

    Mat3 R_po = RotX(-yp) * RotY(-xp) * RotZ(sp);
    return R_po;
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 37
  Vec6 GCRF2ITRF(Real t_tai, const Vec6& rv_gcrf) {
    Mat3 R_po = RotPolarMotion(t_tai);
    Mat3 R_pn = RotPrecessionNutation(t_tai);
    Mat3 R_s = RotSideralMotion(t_tai);
    Mat3 R_s_dot = RotSideralMotionDot(t_tai);

    Mat3 R = R_po * R_s * R_pn;
    Mat3 R_dot = R_po * R_s_dot * R_pn;

    Vec3 r_gcrf = rv_gcrf.head(3);
    Vec3 v_gcrf = rv_gcrf.tail(3);

    Vec3 r_itrf = R * r_gcrf;
    Vec3 v_itrf = R * v_gcrf + R_dot * r_gcrf;

    Vec6 rv_itrf;
    rv_itrf << r_itrf, v_itrf;
    return rv_itrf;
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 37
  Vec6 ITRF2GCRF(Real t_tai, const Vec6& rv_itrf) {
    Mat3 R_po = RotPolarMotion(t_tai);
    Mat3 R_pn = RotPrecessionNutation(t_tai);
    Mat3 R_s = RotSideralMotion(t_tai);
    Mat3 R_s_dot = RotSideralMotionDot(t_tai);

    Mat3 R = R_po * R_s * R_pn;
    Mat3 R_dot = R_po * R_s_dot * R_pn;

    Mat3 R_inv = R.transpose();
    Mat3 R_dot_inv = R_dot.transpose();

    Vec3 r_itrf = rv_itrf.head(3);
    Vec3 v_itrf = rv_itrf.tail(3);

    Vec3 r_gcrf = R_inv * r_itrf;
    Vec3 v_gcrf = R_inv * v_itrf + R_dot_inv * r_itrf;

    Vec6 rv_gcrf;
    rv_gcrf << r_gcrf, v_gcrf;
    return rv_gcrf;
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 38
  Mat3d EarthFrameBiasMatrix() {
    const double da = FRAME_BIAS_DALPHA0;
    const double xi = FRAME_BIAS_XI0;
    const double eta = FRAME_BIAS_ETA0;
    Mat3d B_e = RotX(-eta) * RotY(xi) * RotZ(da);
    return B_e;
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 38
  Mat3d EarthFrameBiasMatrixFirstOrder() {
    const double da = FRAME_BIAS_DALPHA0;
    const double xi = FRAME_BIAS_XI0;
    const double eta = FRAME_BIAS_ETA0;
    Mat3d B_e{{1, da, -xi}, {-da, 1, -eta}, {xi, eta, 1}};
    return B_e;
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 38
  Mat3d EarthFrameBiasMatrixSecondOrder() {
    const double da = FRAME_BIAS_DALPHA0;
    const double xi = FRAME_BIAS_XI0;
    const double eta = FRAME_BIAS_ETA0;
    Mat3d B_e{{1 - 0.5 * (da * da + xi * xi), da, -xi},
              {-da - eta * xi, 1 - 0.5 * (da * da + eta * eta), -eta},
              {xi - eta * da, eta + xi * da, 1 - 0.5 * (eta * eta + xi * xi)}};
    return B_e;
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 39
  Vec6 GCRF2EME(Real t_tai, const Vec6& rv_gcrf) {
    Mat3d B_e = EarthFrameBiasMatrix();
    Vec3 r_gcrf = rv_gcrf.head(3);
    Vec3 v_gcrf = rv_gcrf.tail(3);

    Vec3 r_eme = B_e * r_gcrf;
    Vec3 v_eme = B_e * v_gcrf;

    Vec6 rv_eme;
    rv_eme << r_eme, v_eme;
    return rv_eme;
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 39
  Vec6 EME2GCRF(Real t_tai, const Vec6& rv_eme) {
    Mat3d B_e = EarthFrameBiasMatrix();
    Vec3 r_eme = rv_eme.head(3);
    Vec3 v_eme = rv_eme.tail(3);

    Mat3d B_e_inv = B_e.transpose();
    Vec3 r_gcrf = B_e_inv * r_eme;
    Vec3 v_gcrf = B_e_inv * v_eme;

    Vec6 rv_gcrf;
    rv_gcrf << r_gcrf, v_gcrf;
    return rv_gcrf;
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 39
  Vec6 GCRF2ICRF(Real t_tai, const Vec6& rv_gcrf) {
    Vec6 rv_earth2ssb = GetBodyPosVel(t_tai, NaifId::EARTH, NaifId::SSB);
    Vec6 rv_icrf = rv_gcrf + rv_earth2ssb;
    return rv_icrf;
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 40
  Vec6 ICRF2GCRF(Real t_tai, const Vec6& rv_icrf) {
    Vec6 rv_ssb2earth = GetBodyPosVel(t_tai, NaifId::SSB, NaifId::EARTH);
    Vec6 rv_gcrf = rv_icrf + rv_ssb2earth;
    return rv_gcrf;
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 40
  Vec6 GCRF2MoonCI(Real t_tai, const Vec6& rv_gcrf) {
    Vec6 rv_earth2moon = GetBodyPosVel(t_tai, NaifId::EARTH, NaifId::MOON);
    Vec6 rv_mi = rv_gcrf + rv_earth2moon;
    return rv_mi;
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 40
  Vec6 MoonMI2GCRF(Real t_tai, const Vec6& rv_mi) {
    Vec6 rv_earth2moon = GetBodyPosVel(t_tai, NaifId::EARTH, NaifId::MOON);
    Vec6 rv_gcrf = rv_mi - rv_earth2moon;
    return rv_gcrf;
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 42
  Mat3d RotMItoPA(Real t_tai) {
    double psi = 0;
    double theta = 0;
    double phi = 0;
    Mat3d Rot_mi2pa = RotZ(psi) * RotX(theta) * RotZ(phi);
    return Rot_mi2pa;
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 42
  Mat3d RotMItoPAdot(Real t_tai) {
    double phi = 0;
    double theta = 0;
    double psi = 0;
    double psi_dot = 0;
    double spsi = sin(psi);
    double cpsi = cos(psi);
    Mat3d mat{
        {-psi_dot * spsi, psi_dot * cpsi, 0}, {-psi_dot * cpsi, -psi_dot * spsi, 0}, {0, 0, 0}};
    Mat3d Rot_mi2pa_dot = mat * RotX(theta) * RotZ(phi);
    return Rot_mi2pa_dot;
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 42
  Vec6 MoonMI2MoonPA(Real t_tai, const Vec6& rv_mi) {
    Mat3d Rot_mi2pa = RotMItoPA(t_tai);
    Mat3d Rot_mi2pa_dot = RotMItoPAdot(t_tai);

    Vec3 r_mi = rv_mi.head(3);
    Vec3 v_mi = rv_mi.tail(3);

    Vec3 r_pa = Rot_mi2pa * r_mi;
    Vec3 v_pa = Rot_mi2pa * v_mi + Rot_mi2pa_dot * r_mi;

    Vec6 rv_pa;
    rv_pa << r_pa, v_pa;
    return rv_pa;
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 42
  Vec6 MoonPA2MoonCI(Real t_tai, const Vec6& rv_pa) {
    Mat3d Rot_mi2pa = RotMItoPA(t_tai);
    Mat3d Rot_mi2pa_dot = RotMItoPAdot(t_tai);

    Mat3d Rot_pa2mi = Rot_mi2pa.transpose();
    Mat3d Rot_pa2mi_dot = Rot_mi2pa_dot.transpose();

    Vec3 r_pa = rv_pa.head(3);
    Vec3 v_pa = rv_pa.tail(3);

    Vec3 r_mi = Rot_pa2mi * r_pa;
    Vec3 v_mi = Rot_pa2mi * v_pa + Rot_pa2mi_dot * r_pa;

    Vec6 rv_mi;
    rv_mi << r_mi, v_mi;
    return rv_mi;
  }

  Mat3d MoonFrameBiasMatrix() {
    Mat3d B_moon
        = RotX(-0.2785 * RAD_ARCSEC) * RotY(-78.6944 * RAD_ARCSEC) * RotZ(-67.8526 * RAD_ARCSEC);
    return B_moon;
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 43
  Vec6 MoonPA2MoonME(Real t_tai, const Vec6& rv_pa) {
    Mat3d B_moon = MoonFrameBiasMatrix();
    Vec3 r_pa = rv_pa.head(3);
    Vec3 v_pa = rv_pa.tail(3);

    Vec3 r_me = B_moon * r_pa;
    Vec3 v_me = B_moon * v_pa;

    Vec6 rv_me;
    rv_me << r_me, v_me;
    return rv_me;
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 43
  Vec6 MoonME2MoonPA(Real t_tai, const Vec6& rv_me) {
    Mat3d B_moon = MoonFrameBiasMatrix();
    Vec3 r_me = rv_me.head(3);
    Vec3 v_me = rv_me.tail(3);

    Mat3d B_moon_inv = B_moon.transpose();
    Vec3 r_pa = B_moon_inv * r_me;
    Vec3 v_pa = B_moon_inv * v_me;

    Vec6 rv_pa;
    rv_pa << r_pa, v_pa;
    return rv_pa;
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 48
  Vec6 GCRF2EMR(Real t_tai, const Vec6& rv_gcrf) {
    Vec6 rv_earth2emb = GetBodyPosVel(t_tai, NaifId::EARTH, NaifId::EMB);
    Vec6 rv_emr = Inertial2Synodic(rv_earth2emb, rv_gcrf);
    return rv_emr;
  }

  /// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 48
  Vec6 EMR2GCRF(Real t_tai, const Vec6& rv_emr) {
    Vec6 rv_earth2emb = GetBodyPosVel(t_tai, NaifId::EARTH, NaifId::EMB);
    Vec6 rv_gcrf = Synodic2Intertial(rv_earth2emb, rv_emr);
    return rv_gcrf;
  }

  /// @brief Transform from moon frame to earth frame
  /// @param t_tai
  /// @param rv_me
  /// @return Vec6
  /// @ref
  /// T. A. Ely, ‘Stable Constellations of Frozen Elliptical Inclined Lunar Orbits’, J of
  /// Astronaut Sci, vol. 53, no. 3, pp. 301–316, Sep. 2005, doi: 10.1007/BF03546355.
  Vec6 MoonME2MoonOP(Real t_tai, const Vec6& rv_me) {
    Vec6 rv_m2e = GetBodyPosVel(t_tai, NaifId::MOON, NaifId::EARTH);
    Vec3 r = rv_m2e.head(3);
    Vec3 v = rv_m2e.tail(3);

    Vec3 z_op = r.cross(v).normalized();
    Vec3 i_pole = Vec3::UnitZ();
    Vec3 x_op = i_pole.cross(z_op).normalized();
    Vec3 y_op = z_op.cross(x_op).normalized();

    Mat3 R_op2me;
    R_op2me << x_op, y_op, z_op;
    Mat3 R_me2op = R_op2me.transpose();

    Vec3 r_op = R_me2op * rv_me.head(3);
    Vec3 v_op = R_me2op * rv_me.tail(3);
    Vec6 rv_op;
    rv_op << r_op, v_op;
    return rv_op;
  }

  /// @brief Transform from moon frame to earth frame
  /// @param t_tai
  /// @param rv_op
  /// @return Vec6
  /// @ref
  /// T. A. Ely, ‘Stable Constellations of Frozen Elliptical Inclined Lunar Orbits’, J of
  /// Astronaut Sci, vol. 53, no. 3, pp. 301–316, Sep. 2005, doi: 10.1007/BF03546355.
  Vec6 MoonOP2MoonME(Real t_tai, const Vec6& rv_op) {
    Vec6 rv_m2e = GetBodyPosVel(t_tai, NaifId::MOON, NaifId::EARTH);
    Vec3 r = rv_m2e.head(3);
    Vec3 v = rv_m2e.tail(3);

    Vec3 z_op = r.cross(v).normalized();
    Vec3 i_pole = Vec3::UnitZ();
    Vec3 x_op = i_pole.cross(z_op).normalized();
    Vec3 y_op = z_op.cross(x_op).normalized();

    Mat3 R_op2me;
    R_op2me << x_op, y_op, z_op;
    Mat3 R_me2op = R_op2me.transpose();

    Vec3 r_me = R_op2me * rv_op.head(3);
    Vec3 v_me = R_op2me * rv_op.tail(3);
    Vec6 rv_me;
    rv_me << r_me, v_me;
    return rv_me;
  }

}  // namespace lupnt

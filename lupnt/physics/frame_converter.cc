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

#include "frame_converter.h"

#include <filesystem>

#include "cheby.h"
#include "frame_converter_spice.h"
#include "lupnt/core/constants.h"
#include "lupnt/data/eop.h"
#include "lupnt/data/iau_sofa.h"
#include "lupnt/data/kernels.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/orbit_state.h"
#include "lupnt/physics/time_converter.h"

#define FRAME_CONVERSION(from, to, func) \
  {{Frame::from, Frame::to},             \
   [](Real t, const Vec6& rv) -> Vec6 { return func(t, rv); }}

namespace lupnt {

std::map<std::pair<Frame, Frame>, std::function<Vec6(Real, const Vec6& rv)>>
    frame_conversions = {FRAME_CONVERSION(GCRF, ITRF, GCRFtoITRF),
                         FRAME_CONVERSION(ITRF, GCRF, ITRFtoGCRF),
                         FRAME_CONVERSION(GCRF, EME, GCRFtoEME),
                         FRAME_CONVERSION(EME, GCRF, EMEtoGCRF),
                         FRAME_CONVERSION(GCRF, ICRF, GCRFtoICRF),
                         FRAME_CONVERSION(ICRF, GCRF, ICRFtoGCRF),
                         FRAME_CONVERSION(GCRF, MOON_CI, GCRFtoMI),
                         FRAME_CONVERSION(MOON_CI, GCRF, MItoGCRF),
                         FRAME_CONVERSION(MOON_CI, MOON_PA, MItoPA),
                         FRAME_CONVERSION(MOON_PA, MOON_CI, PAtoMI),
                         FRAME_CONVERSION(MOON_PA, MOON_ME, PAtoME),
                         FRAME_CONVERSION(MOON_ME, MOON_PA, MEtoPA),
                         FRAME_CONVERSION(GCRF, EMR, GCRFtoEMR),
                         FRAME_CONVERSION(EMR, GCRF, EMRtoGCRF)};

/// @brief Convert the state vector from one coordinate system to another
/// (with integer ID input)
/// @param t_tai Epoch of the state vector
/// @param rv_in   State vector in the original coordinate system
/// @param frame_in Coordinate system of the original state vector
/// @param frame_out Coordinate system of the converted state vector
/// @return Vec6  State vector in the converted coordinate system
///
Vec6 ConvertFrame(Real t_tai, const Vec6& rv_in, Frame frame_in,
                  Frame frame_out) {
  if (frame_in == frame_out) return rv_in;
  std::vector<Frame> path =
      FindShortestPath(frame_in, frame_out, frame_conversions);
  Vec6 rv_out = rv_in;
  for (size_t i = 0; i < path.size() - 1; i++) {
    rv_out = frame_conversions[{path[i], path[i + 1]}](t_tai, rv_out);
  }
  return rv_out;
}

CartesianOrbitState ConvertFrame(Real t_tai,
                                 const CartesianOrbitState& state_in,
                                 Frame frame_out) {
  Vec6 rv_in = state_in.GetVec();
  Vec6 rv_out =
      ConvertFrame(t_tai, rv_in, state_in.GetCoordSystem(), frame_out);
  return CartesianOrbitState(rv_out, frame_out);
}

Vec3 ConvertFrame(Real t_tai, const Vec3& r_in, Frame frame_in,
                  Frame frame_out) {
  Vec6 rv_in;
  rv_in << r_in, Vec3::Zero();
  Vec6 rv_out_6 = ConvertFrame(t_tai, rv_in, frame_in, frame_out);
  return rv_out_6.head(3);
}

Mat<-1, 6> ConvertFrame(Real t_tai, const Mat<-1, 6>& rv_in, Frame frame_in,
                        Frame frame_out) {
  Mat<-1, 6> rv_out(rv_in.rows(), 6);
  for (int i = 0; i < rv_in.rows(); i++) {
    rv_out.row(i) = ConvertFrame(t_tai, rv_in.row(i).transpose().eval(),
                                 frame_in, frame_out);
  }
  return rv_out;
}

Mat<-1, 3> ConvertFrame(Real t_tai, const Mat<-1, 3>& r_in, Frame frame_in,
                        Frame frame_out) {
  Mat<-1, 6> rv_in(r_in.rows(), 6);
  rv_in << r_in, Mat<-1, 3>::Zero(r_in.rows(), 3);
  Mat<-1, 6> rv_out = ConvertFrame(t_tai, rv_in, frame_in, frame_out);
  return rv_out.leftCols(3);
}

Mat<-1, 6> ConvertFrame(VecX t_tai, const Vec6& rv_in, Frame frame_in,
                        Frame frame_out) {
  Mat<-1, 6> rv_out(t_tai.size(), 6);
  for (int i = 0; i < t_tai.size(); i++) {
    rv_out.row(i) = ConvertFrame(t_tai(i), rv_in, frame_in, frame_out);
  }
  return rv_out;
}

Mat<-1, 3> ConvertFrame(VecX t_tai, const Vec3& r_in, Frame frame_in,
                        Frame frame_out) {
  Mat<-1, 6> rv_in(t_tai.size(), 6);
  rv_in << r_in, Mat<-1, 3>::Zero(t_tai.size(), 3);
  Mat<-1, 6> rv_out = ConvertFrame(t_tai, rv_in, frame_in, frame_out);
  return rv_out.leftCols(3);
}

Mat<-1, 6> ConvertFrame(VecX t_tai, const Mat<-1, 6>& rv_in, Frame frame_in,
                        Frame frame_out) {
  assert(t_tai.size() == rv_in.rows() && "Epoch and rv_in must have same size");
  Mat<-1, 6> rv_out(t_tai.size(), 6);
  for (int i = 0; i < t_tai.size(); i++) {
    rv_out.row(i) = ConvertFrame(t_tai(i), rv_in.row(i).transpose().eval(),
                                 frame_in, frame_out);
  }
  return rv_out;
}

Mat<-1, 3> ConvertFrame(VecX t_tai, const Mat<-1, 3>& r_in, Frame frame_in,
                        Frame frame_out) {
  assert(t_tai.size() == r_in.rows() && "Epoch and r_in must have same size");
  Mat<-1, 6> rv_in(t_tai.size(), 6);
  rv_in << r_in, Mat<-1, 3>::Zero(t_tai.size(), 3);
  Mat<-1, 6> rv_out = ConvertFrame(t_tai, rv_in, frame_in, frame_out);
  return rv_out.leftCols(3);
}

Mat6 RotOp2Mi(Real t_tai) {
  Vec6 rv_earth_icrf = GetBodyPosVel(t_tai, NaifId::SUN, NaifId::EARTH);
  Vec6 rv_moon_icrf = GetBodyPosVel(t_tai, NaifId::SUN, NaifId::MOON);

  // IAU pole
  Vec3 iau_pole = Vec3::UnitZ();
  Vec3 iau_pole_icrf =
      ConvertFrame(t_tai, iau_pole, Frame::MOON_PA, Frame::GCRF);

  // MOON_OP unit vectors
  Vec3 dr = rv_earth_icrf.head(3) - rv_moon_icrf.head(3);
  Vec3 dv = rv_earth_icrf.tail(3) - rv_moon_icrf.tail(3);
  Vec3 z_op = dr.cross(dv).normalized();
  Vec3 x_op = iau_pole_icrf.cross(z_op).normalized();
  Vec3 y_op = z_op.cross(x_op).normalized();

  // create rotation matrix from MOON_OP to MOON_CI
  Mat3 R_op2mi;
  R_op2mi << x_op, y_op, z_op;

  Mat6 R_op2mi_tot = Mat6::Zero();
  R_op2mi_tot.topLeftCorner(3, 3) = R_op2mi;
  R_op2mi_tot.bottomRightCorner(3, 3) = R_op2mi;

  return R_op2mi_tot.cast<double>();
}

/// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 34
Mat3 RotPrecessionNutation(Real t_tai) {
  Real t_tt = ConvertTime(t_tai, TimeSys::TAI, TimeSys::TT);
  Real jd_tt = TimeToJD(t_tt);

  IauSofaData iau_data = GetIauSofaData(jd_tt);
  Real X = iau_data.X;
  Real Y = iau_data.Y;
  Real s = iau_data.s;

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
Vec6 GCRFtoITRF(Real t_tai, const Vec6& rv_gcrf) {
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
Vec6 ITRFtoGCRF(Real t_tai, const Vec6& rv_itrf) {
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
Vec6 GCRFtoEME(Real t_tai, const Vec6& rv_gcrf) {
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
Vec6 EMEtoGCRF(Real t_tai, const Vec6& rv_eme) {
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
Vec6 GCRFtoICRF(Real t_tai, const Vec6& rv_gcrf) {
  Vec6 rv_earth2ssb = GetBodyPosVel(t_tai, NaifId::EARTH, NaifId::SSB);
  Vec6 rv_icrf = rv_gcrf + rv_earth2ssb;
  return rv_icrf;
}

/// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 40
Vec6 ICRFtoGCRF(Real t_tai, const Vec6& rv_icrf) {
  Vec6 rv_ssb2earth = GetBodyPosVel(t_tai, NaifId::SSB, NaifId::EARTH);
  Vec6 rv_gcrf = rv_icrf + rv_ssb2earth;
  return rv_gcrf;
}

/// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 40
Vec6 GCRFtoMI(Real t_tai, const Vec6& rv_gcrf) {
  Vec6 rv_earth2moon = GetBodyPosVel(t_tai, NaifId::EARTH, NaifId::MOON);
  Vec6 rv_mi = rv_gcrf + rv_earth2moon;
  return rv_mi;
}

/// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 40
Vec6 MItoGCRF(Real t_tai, const Vec6& rv_mi) {
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
  Mat3d mat{{-psi_dot * spsi, psi_dot * cpsi, 0},
            {-psi_dot * cpsi, -psi_dot * spsi, 0},
            {0, 0, 0}};
  Mat3d Rot_mi2pa_dot = mat * RotX(theta) * RotZ(phi);
  return Rot_mi2pa_dot;
}

/// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 42
Vec6 MItoPA(Real t_tai, const Vec6& rv_mi) {
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
Vec6 PAtoMI(Real t_tai, const Vec6& rv_pa) {
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
  Mat3d B_moon = RotX(-0.2785 * RAD_ARCSEC) * RotY(-78.6944 * RAD_ARCSEC) *
                 RotZ(-67.8526 * RAD_ARCSEC);
  return B_moon;
}

/// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 43
Vec6 PAtoME(Real t_tai, const Vec6& rv_pa) {
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
Vec6 MEtoPA(Real t_tai, const Vec6& rv_me) {
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
Vec6 GCRFtoEMR(Real t_tai, const Vec6& rv_gcrf) {
  Vec6 rv_earth2emb = GetBodyPosVel(t_tai, NaifId::EARTH, NaifId::EMB);
  Vec6 rv_emr = Inertial2Synodic(rv_earth2emb, rv_gcrf);
  return rv_emr;
}

/// @ref Astrodynamics Convention & Modeling Reference, Version 1.1, Page 48
Vec6 EMRtoGCRF(Real t_tai, const Vec6& rv_emr) {
  Vec6 rv_earth2emb = GetBodyPosVel(t_tai, NaifId::EARTH, NaifId::EMB);
  Vec6 rv_gcrf = Synodic2Intertial(rv_earth2emb, rv_emr);
  return rv_gcrf;
}

}  // namespace lupnt
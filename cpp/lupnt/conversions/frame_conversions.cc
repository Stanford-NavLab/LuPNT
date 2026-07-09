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

#include "lupnt/conversions/frame_conversions.h"

#include <cspice/SpiceUsr.h>

#include <cmath>
#include <functional>
#include <optional>
#include <vector>

#include "lupnt/conversions/frame_converter.h"  // Frame enum + IsPlanet*Frame predicates
#include "lupnt/conversions/state_conversions.h"
#include "lupnt/conversions/time_conversions.h"
#include "lupnt/core/constants.h"
#include "lupnt/core/logger.h"
#include "lupnt/interfaces/eop.h"
#include "lupnt/interfaces/iau_sofa.h"
#include "lupnt/interfaces/kernels.h"
#include "lupnt/interfaces/spice.h"
#include "lupnt/numerics/cheby_fit.h"
#include "lupnt/numerics/math_utils.h"

#define FRAME_CONVERSION(from, to, func) \
  {{Frame::from, Frame::to}, [](Real t, const Vec6& rv) -> Vec6 { return func(t, rv); }}

namespace lupnt {

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 34
  Mat3 RotPrecessionNutation(Real t_tdb) {
    Real t_tt = ConvertTime(t_tdb, Time::TDB, Time::TT);
    Real jd_tt = TimeToJd(t_tt);

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

  namespace {
    // Forward declaration; defined below alongside the SPICE-fitted EOP
    // infrastructure. Lets RotSideralMotion/RotSideralMotionDot/
    // RotPolarMotion transparently consult a SPICE-derived EOP fit (see
    // InitFrameConversionFromSpice / ComputeEopFromSpice) before falling
    // back to the loaded IERS EOP table (GetEopData).
    bool TryGetFittedEop(Real t_tdb, Real* xp, Real* yp, Real* ut1_utc, Real* ut1_utc_dot);
  }  // namespace

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 35
  ///
  /// If InitFrameConversionFromSpice has fitted a SPICE-derived EOP model
  /// (see ComputeEopFromSpice) covering `t_tdb`, the fitted UT1-UTC is used
  /// in place of the loaded IERS EOP table.
  Mat3 RotSideralMotion(Real t_tdb) {
    Real t_ut1;
    Real ut1_utc;
    if (TryGetFittedEop(t_tdb, nullptr, nullptr, &ut1_utc, nullptr)) {
      Real t_utc = ConvertTime(t_tdb, Time::TDB, Time::UTC);
      t_ut1 = t_utc + ut1_utc;
    } else {
      t_ut1 = ConvertTime(t_tdb, Time::TDB, Time::UT1);
    }
    Real theta_era = EarthRotationAngle(t_ut1);
    Mat3 R_s = RotZ(theta_era);
    return R_s;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 36
  ///
  /// If InitFrameConversionFromSpice has fitted a SPICE-derived EOP model
  /// (see ComputeEopFromSpice) covering `t_tdb`, the fitted UT1-UTC (and its
  /// *exact* analytic time derivative) is used to compute theta_era and the
  /// Earth rotation rate w_E -- consistently with RotSideralMotion -- in
  /// place of the loaded IERS EOP table's UT1-UTC/LOD.
  Mat3 RotSideralMotionDot(Real t_tdb) {
    Real t_ut1;
    Real w_E;

    Real ut1_utc, ut1_utc_dot;
    if (TryGetFittedEop(t_tdb, nullptr, nullptr, &ut1_utc, &ut1_utc_dot)) {
      Real t_utc = ConvertTime(t_tdb, Time::TDB, Time::UTC);
      t_ut1 = t_utc + ut1_utc;
      w_E = (TWO_PI * 1.00273781191135448 / SECS_DAY) * (1. + ut1_utc_dot);
    } else {
      t_ut1 = ConvertTime(t_tdb, Time::TDB, Time::UT1);

      Real t_utc = ConvertTime(t_tdb, Time::TDB, Time::UTC);
      Real mjd_utc = TimeToMjd(t_utc);
      EopData eop = GetEopData(mjd_utc);
      Real lod = eop.lod;
      w_E = 7.292115146706979e-5 * (1. - lod / SECS_DAY);
    }

    Real theta_era = EarthRotationAngle(t_ut1);
    return RotZdot(theta_era, w_E);
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 36
  ///
  /// If InitFrameConversionFromSpice has fitted a SPICE-derived EOP model
  /// (see ComputeEopFromSpice) covering `t_tdb`, the fitted x_pole/y_pole are
  /// used in place of the loaded IERS EOP table.
  Mat3 RotPolarMotion(Real t_tdb) {
    Real t_tt = ConvertTime(t_tdb, Time::TDB, Time::TT);

    Real xp, yp;
    if (!TryGetFittedEop(t_tdb, &xp, &yp, nullptr, nullptr)) {
      Real t_utc = ConvertTime(t_tdb, Time::TDB, Time::UTC);
      Real mjd_utc = TimeToMjd(t_utc);
      EopData eop = GetEopData(mjd_utc);
      xp = eop.x_pole;
      yp = eop.y_pole;
    }
    Real sp = -47e-6 * RAD_ARCSEC * (t_tt / DAYS_CENTURY);

    Mat3 R_po = RotX(-yp) * RotY(-xp) * RotZ(sp);
    return R_po;
  }

  // =====================================================================
  // SPICE-derived Earth orientation parameters & SPICE-fitted lunar
  // orientation
  // =====================================================================
  //
  // GcrfToItrf/ItrfToGcrf and MoonCiToPa/MoonPaToCi normally compute the
  // Earth/lunar orientation rotations from purely-analytic models (IAU
  // 2006/2000A precession-nutation + EOP-based polar motion/sidereal
  // rotation for Earth -- see RotPrecessionNutation/RotPolarMotion/
  // RotSideralMotion above, with EOP (x_pole, y_pole, UT1-UTC, LOD) coming
  // from the loaded IERS EOP table, see lupnt/interfaces/eop.h; static,
  // pre-extracted DE Chebyshev "lunar mantle libration" angles for the
  // Moon -- see GetLunarMantleData in kernels.cc). These can disagree with
  // SPICE's high-accuracy binary PCK-based orientation
  // (lupnt::spice::ConvertFrameSpice, see frame_converter_spice.cc -- now
  // backed by the auto-downloaded earth_latest_high_prec.bpc /
  // moon_pa_de*.bpc + moon_de*.tf kernels, see spice::LoadSpiceKernel) by a
  // non-negligible amount -- in particular when the loaded EOP table does
  // not cover the requested epoch.
  //
  // ComputeEopFromSpice(t_tdb) (declared in frame_conversions.h) derives the
  // Earth orientation parameters (x_pole, y_pole, UT1-UTC) at a single epoch
  // directly from SPICE's "J2000"<->"ITRF93" rotation, by combining it with
  // LuPNT's analytic precession-nutation matrix (RotPrecessionNutation) and
  // inverting the polar-motion/sidereal-rotation decomposition
  // R_spice = R_po(xp,yp) * R_s(theta_era) * R_pn (closed form, see the
  // implementation below).
  //
  // InitFrameConversionFromSpice() (declared in frame_conversions.h) uses
  // this to close the Earth-orientation gap for a known simulation time
  // window: it samples ComputeEopFromSpice at a dense set of Chebyshev-Gauss
  // epochs spanning the window and fits piecewise Chebyshev polynomial
  // coefficients to the 3 EOP parameters -- the very same Chebyshev
  // representation (and evaluator, cheby_eval_ad) used internally for the
  // planetary/lunar position ephemeris data, see spice_cheby.h / kernels.cc.
  // Once initialized, RotPolarMotion/RotSideralMotion/RotSideralMotionDot
  // (and hence GcrfToItrf/ItrfToGcrf) transparently use these fitted EOP
  // values -- together with the *exact* analytic time derivative of
  // UT1-UTC, which a Chebyshev polynomial yields for free -- whenever the
  // requested epoch falls inside the fitted window, reproducing SPICE
  // closely without requiring SPICE at evaluation time. Outside of the
  // fitted window, the loaded IERS EOP table is used unchanged, so this
  // feature is purely additive/opt-in (nothing changes unless
  // InitFrameConversionFromSpice is explicitly called).
  //
  // The Moon CI<->PA fit (ZXZ Euler angles fitted directly to the SPICE
  // "J2000"<->"IAU_MOON" rotation) is unrelated and unchanged -- see
  // MoonCiToPa/MoonPaToCi and FitEulerAngleModel below.

  namespace {

    // ChebyshevFitSegment, ChebyshevFitModel, and FitChebyshevModel are now
    // defined in lupnt/numerics/cheby_fit.h (shared with time_conversions.cc).

    // A Chebyshev-fitted 3x3 rotation matrix R(t) together with its exact
    // Extracts ZXZ Euler angles (phi, theta, psi) from a rotation matrix R
    // satisfying R = Rz(psi) * Rx(theta) * Rz(phi) (passive convention,
    // matching RotZ / RotX in math_utils.cc).  Returns phi, psi in (-pi,pi]
    // and theta in [0, pi].
    //
    // Derived from the ZXZ product (cp=cos(psi), sp=sin(psi), etc.):
    //   R(2,2) = cos(theta)
    //   R(2,0) =  sin(phi)*sin(theta)   R(2,1) = -cos(phi)*sin(theta)
    //   R(0,2) =  sin(psi)*sin(theta)   R(1,2) =  cos(psi)*sin(theta)
    // => theta = acos(R(2,2)), phi = atan2(R(2,0), -R(2,1)),
    //    psi = atan2(R(0,2), R(1,2)).
    // Near gimbal lock (theta~0 or ~pi, where Earth's GCRF<->ITRF can be
    // transiently for very short windows), phi and psi become individually
    // sensitive, but the reconstructed R = Rz(psi)*Rx(theta)*Rz(phi) stays
    // exactly orthonormal regardless; in practice theta is always nonzero.
    inline void ExtractZXZAngles(const Mat3d& R, double& phi, double& theta, double& psi) {
      theta = std::acos(std::clamp(R(2, 2), -1.0, 1.0));
      phi = std::atan2(R(2, 0), -R(2, 1));
      psi = std::atan2(R(0, 2), R(1, 2));
    }

    // Adjusts `raw` by the nearest integer multiple of 2*pi so the result
    // is within pi of `ref` (standard single-step phase unwrap).
    inline double UnwrapAngle(double raw, double ref) {
      const double TWO_PI = 2.0 * PI;
      double diff = raw - ref;
      return raw - std::round(diff / TWO_PI) * TWO_PI;
    }

    // A Chebyshev-fitted rotation model that stores R(t) as ZXZ Euler angles
    // (phi, theta, psi) and reconstructs R = Rz(psi)*Rx(theta)*Rz(phi) at
    // evaluation time.
    //
    // KEY GUARANTEE: the reconstructed R_fit is ALWAYS exactly orthonormal
    // (to floating-point precision) regardless of the Chebyshev fit residual
    // in the angles -- because RotZ/RotX each satisfy R^T*R = I exactly via
    // sin^2+cos^2=1, and products of orthonormal matrices are orthonormal.
    // This eliminates the ~4e-8 orthonormality violation that occurred when
    // fitting the 9 matrix elements independently.
    struct FittedEulerAngleModel {
      ChebyshevFitModel fit;  // dim 0=phi, 1=theta, 2=psi

      bool Covers(Real t_tdb) const {
        double t = t_tdb.val();
        return t >= fit.t_start && t <= fit.t_end;
      }

      // Evaluates R(t) = Rz(psi)*Rx(theta)*Rz(phi) -- always exactly
      // orthonormal -- and, if R_dot != nullptr, its exact analytic dR/dt
      // via the product rule (same pattern as RotMoonCiToPa in kernels.cc).
      bool Eval(Real t_tdb, Mat3* R, Mat3* R_dot) const {
        VecX angles, angles_dot;
        if (!fit.Eval(t_tdb, &angles, R_dot != nullptr ? &angles_dot : nullptr)) return false;

        Real phi = angles(0), theta = angles(1), psi = angles(2);
        Mat3 R_z_phi = RotZ(phi);
        Mat3 R_x_theta = RotX(theta);
        Mat3 R_z_psi = RotZ(psi);
        *R = R_z_psi * R_x_theta * R_z_phi;  // always exactly orthonormal

        if (R_dot != nullptr) {
          Real phi_dot = angles_dot(0), theta_dot = angles_dot(1), psi_dot = angles_dot(2);
          *R_dot = RotZdot(psi, psi_dot) * R_x_theta * R_z_phi
                   + R_z_psi * RotXdot(theta, theta_dot) * R_z_phi
                   + R_z_psi * R_x_theta * RotZdot(phi, phi_dot);
        }
        return true;
      }
    };

    // Samples rotation matrices from `rotation_sampler` at Chebyshev-Gauss
    // nodes across [t_start, t_end] (split into ~segment_length segments),
    // extracts ZXZ Euler angles (phi, theta, psi) for R = Rz(psi)*Rx(theta)*
    // Rz(phi), UNWRAPS phi and psi continuously across the entire window so
    // the fitted functions are smooth, and returns a 3-dim piecewise
    // Chebyshev model ([dim 0]=phi, [dim 1]=theta, [dim 2]=psi).
    //
    // Unwrapping is essential for phi and psi: Earth's dominant spin angle
    // (psi~ERA) increases by ~2*pi per sidereal day, so without continuous
    // unwrapping it would wrap at +/-pi and break the polynomial fit.  Since
    // RotZ/RotX are 2*pi-periodic, R_fit is identical for wrapped or
    // unwrapped values; unwrapping only makes the polynomial smooth.
    ChebyshevFitModel FitEulerAngleModel(const std::function<Mat3d(double)>& rotation_sampler,
                                         double t_start, double t_end, double segment_length,
                                         int num_coeffs) {
      if (!(t_end > t_start))
        throw std::runtime_error("FitEulerAngleModel: t_end must be greater than t_start");
      if (num_coeffs < 2)
        throw std::runtime_error("FitEulerAngleModel: num_coeffs must be at least 2");
      if (!(segment_length > 0.0))
        throw std::runtime_error("FitEulerAngleModel: segment_length must be positive");

      double span = t_end - t_start;
      int num_segments = std::max(1, static_cast<int>(std::ceil(span / segment_length)));
      double seg_len = span / num_segments;

      ChebyshevFitModel model;
      model.t_start = t_start;
      model.t_end = t_end;
      model.num_dims = 3;  // [phi, theta, psi]
      model.num_coeffs = num_coeffs;
      model.segments.reserve(num_segments);

      // Running reference for inter-segment continuity: unwrapped phi, psi
      // of the most-recently-processed sample in temporal order.
      double ref_phi = 0.0, ref_psi = 0.0;
      bool first_sample = true;

      for (int s = 0; s < num_segments; s++) {
        double seg_start = t_start + s * seg_len;
        double seg_end = (s == num_segments - 1) ? t_end : seg_start + seg_len;
        double mid = 0.5 * (seg_start + seg_end);
        double radius = 0.5 * (seg_end - seg_start);

        // Chebyshev-Gauss nodes: t_k = mid + radius*cos((2k+1)*pi/(2N)),
        // k=0..N-1.  cos is strictly decreasing in k, so k=0 is the LATEST
        // node and k=N-1 is the EARLIEST.  Process in temporal order
        // (time-sorted index i=0..N-1 -> DCT index k=N-1-i) so ref_phi/psi
        // advances monotonically in time for correct unwrapping.
        std::vector<double> phi_dct(num_coeffs), theta_dct(num_coeffs), psi_dct(num_coeffs);

        for (int i = 0; i < num_coeffs; i++) {
          int k = num_coeffs - 1 - i;  // DCT index for the i-th earliest node
          double cheby_arg = (2.0 * k + 1.0) * PI / (2.0 * num_coeffs);
          double t = mid + radius * std::cos(cheby_arg);

          Mat3d Rmat = rotation_sampler(t);
          double phi_raw, theta_val, psi_raw;
          ExtractZXZAngles(Rmat, phi_raw, theta_val, psi_raw);

          double phi_uw, psi_uw;
          if (first_sample) {
            phi_uw = phi_raw;
            psi_uw = psi_raw;
            first_sample = false;
          } else {
            phi_uw = UnwrapAngle(phi_raw, ref_phi);
            psi_uw = UnwrapAngle(psi_raw, ref_psi);
          }
          ref_phi = phi_uw;
          ref_psi = psi_uw;

          phi_dct[k] = phi_uw;
          theta_dct[k] = theta_val;  // theta in [0,pi], no wrap possible
          psi_dct[k] = psi_uw;
        }

        // Discrete Chebyshev transform (DCT-II) -- identical formula to
        // FitChebyshevModel above.
        ChebyshevFitSegment seg;
        seg.t_mid = mid;
        seg.t_radius = radius;
        seg.coeffs.assign(3, std::vector<double>(num_coeffs, 0.0));
        for (int j = 0; j < num_coeffs; j++) {
          double w = (j == 0) ? 1.0 : 2.0;
          for (int k = 0; k < num_coeffs; k++) {
            double c = std::cos((2.0 * k + 1.0) * PI * j / (2.0 * num_coeffs));
            seg.coeffs[0][j] += phi_dct[k] * c;
            seg.coeffs[1][j] += theta_dct[k] * c;
            seg.coeffs[2][j] += psi_dct[k] * c;
          }
          seg.coeffs[0][j] *= w / num_coeffs;
          seg.coeffs[1][j] *= w / num_coeffs;
          seg.coeffs[2][j] *= w / num_coeffs;
        }
        model.segments.push_back(std::move(seg));
      }
      return model;
    }

    // Returns the 3x3 rotation block of the SPICE state-transformation matrix
    // `sxform_c(from_frame, to_frame, et)` at epoch t_tdb (seconds TDB since
    // J2000).
    Mat3d SampleSpiceRotationMatrix(double t_tdb, const std::string& from_frame,
                                    const std::string& to_frame) {
      SpiceDouble et = t_tdb;
      double xform[6][6];
#pragma omp critical
      { sxform_c(from_frame.c_str(), to_frame.c_str(), et, xform); }
      Mat3d R;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) R(i, j) = xform[i][j];
      return R;
    }

    std::optional<ChebyshevFitModel> spice_fit_eop;  // dims: 0=x_pole, 1=y_pole, 2=ut1_utc
    std::optional<FittedEulerAngleModel> spice_fit_moonci_to_moonpa;

    bool TryGetFittedEop(Real t_tdb, Real* xp, Real* yp, Real* ut1_utc, Real* ut1_utc_dot) {
      VecX f, df;
      bool ok = false;
#pragma omp critical
      {
        ok = spice_fit_eop.has_value()
             && spice_fit_eop->Eval(t_tdb, &f, ut1_utc_dot != nullptr ? &df : nullptr);
      }
      if (!ok) return false;
      if (xp != nullptr) *xp = f(0);
      if (yp != nullptr) *yp = f(1);
      if (ut1_utc != nullptr) *ut1_utc = f(2);
      if (ut1_utc_dot != nullptr) *ut1_utc_dot = df(2);
      return true;
    }

    bool TryGetFittedMoonCiToPa(Real t_tdb, Mat3* R, Mat3* R_dot) {
      bool ok = false;
#pragma omp critical
      {
        ok = spice_fit_moonci_to_moonpa.has_value()
             && spice_fit_moonci_to_moonpa->Eval(t_tdb, R, R_dot);
      }
      return ok;
    }

  }  // namespace

  SpiceEopParams ComputeEopFromSpice(Real t_tdb) {
    spice::LoadSpiceKernel();

    Mat3d R_spice = SampleSpiceRotationMatrix(t_tdb.val(), "J2000", "ITRF93");

    Mat3 R_pn = RotPrecessionNutation(t_tdb);
    Mat3d R_pn_d;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) R_pn_d(i, j) = R_pn(i, j).val();

    // R_spice ~= R_po(xp,yp) * R_s(theta_era) * R_pn
    //          = RotX(-yp) * RotY(-xp) * RotZ(sp + theta_era) * R_pn
    Mat3d M = R_spice * R_pn_d.transpose();

    double xp = std::asin(std::clamp(M(0, 2), -1.0, 1.0));
    double yp = -std::atan2(M(1, 2), M(2, 2));
    double psi = std::atan2(M(0, 1), M(0, 0));  // = sp + theta_era

    double t_tt = ConvertTime(t_tdb, Time::TDB, Time::TT).val();
    double sp = -47e-6 * RAD_ARCSEC * (t_tt / DAYS_CENTURY);
    double theta_era = std::atan2(std::sin(psi - sp), std::cos(psi - sp));

    // Invert EarthRotationAngle(t_ut1) for t_ut1 by linearizing about
    // t_ut1 ~= t_tdb (the true UT1-UTC offset is at most ~1 s, far smaller
    // than the ~1-day period of theta_era).
    constexpr double theta_0 = 0.7790572732640;
    constexpr double dtheta_dt = 1.00273781191135448;
    double era_at_tdb_arg = TWO_PI * (theta_0 + dtheta_dt * t_tdb.val() / SECS_DAY);
    double era_at_tdb = std::atan2(std::sin(era_at_tdb_arg), std::cos(era_at_tdb_arg));
    double dtheta = std::atan2(std::sin(theta_era - era_at_tdb), std::cos(theta_era - era_at_tdb));
    double dt_ut1 = dtheta / (TWO_PI * dtheta_dt / SECS_DAY);

    double t_utc = ConvertTime(t_tdb, Time::TDB, Time::UTC).val();
    double ut1_utc = (t_tdb.val() + dt_ut1) - t_utc;

    SpiceEopParams params;
    params.x_pole = xp;
    params.y_pole = yp;
    params.ut1_utc = ut1_utc;
    return params;
  }

  void InitFrameConversionFromSpice(Real t_start_tdb, Real t_end_tdb, Real segment_length_s,
                                    int num_coeffs) {
    double t0 = t_start_tdb.val();
    double t1 = t_end_tdb.val();
    if (!(t1 > t0))
      throw std::runtime_error(
          "InitFrameConversionFromSpice: t_end_tdb must be greater than t_start_tdb");

    spice::LoadSpiceKernel();

    // The native GcrfToItrf/MoonCiToPa functions ultimately need to match
    // lupnt::spice::ConvertFrameSpice -- the SPICE-based reference path
    // already wired into LuPNT (see frame_converter_spice.cc) -- so the
    // fits are computed against the *exact same* SPICE frame pairs that
    // ConvertFrameSpice itself uses for these conversions: "J2000" <->
    // "ITRF93" for GCRF<->ITRF (the high-accuracy Earth orientation frame
    // defined by the auto-downloaded earth_latest_high_prec.bpc), and
    // "J2000" <-> "IAU_MOON" for MOON_CI<->MOON_PA.
    double seg_len = segment_length_s.val();
    ChebyshevFitModel earth_fit = FitChebyshevModel(
        [](double t) -> VecXd {
          SpiceEopParams p = ComputeEopFromSpice(Real(t));
          VecXd f(3);
          f << p.x_pole.val(), p.y_pole.val(), p.ut1_utc.val();
          return f;
        },
        t0, t1, /*num_dims=*/3, seg_len, num_coeffs);
    ChebyshevFitModel moon_fit = FitEulerAngleModel(
        [](double t) { return SampleSpiceRotationMatrix(t, "J2000", "IAU_MOON"); }, t0, t1, seg_len,
        num_coeffs);

    Logger::Info(
        fmt::format("Fitted SPICE-derived Earth orientation parameters (x_pole, y_pole, "
                    "UT1-UTC) and MOON_CI<->MOON_PA orientation over epoch window [{:.3f}, "
                    "{:.3f}] s TDB ({} segment(s) of ~{:.0f} s, {} Chebyshev coefficients per "
                    "fitted quantity)",
                    t0, t1, earth_fit.segments.size(), seg_len, num_coeffs),
        "FrameConverter");

#pragma omp critical
    {
      spice_fit_eop = std::move(earth_fit);
      spice_fit_moonci_to_moonpa = FittedEulerAngleModel{std::move(moon_fit)};
    }
  }

  void ClearFrameConversionFit() {
#pragma omp critical
    {
      spice_fit_eop.reset();
      spice_fit_moonci_to_moonpa.reset();
    }
  }

  bool HasFittedEarthOrientation(Real t_tdb) {
    bool covers = false;
    double t = t_tdb.val();
#pragma omp critical
    {
      covers
          = spice_fit_eop.has_value() && t >= spice_fit_eop->t_start && t <= spice_fit_eop->t_end;
    }
    return covers;
  }

  bool HasFittedLunarOrientation(Real t_tdb) {
    bool covers = false;
#pragma omp critical
    {
      covers = spice_fit_moonci_to_moonpa.has_value() && spice_fit_moonci_to_moonpa->Covers(t_tdb);
    }
    return covers;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 37
  ///
  /// If InitFrameConversionFromSpice has fitted a SPICE-derived EOP model
  /// (see ComputeEopFromSpice) covering `t_tdb`, RotPolarMotion/
  /// RotSideralMotion[Dot] transparently use the fitted EOP, so that this
  /// matches lupnt::spice::ConvertFrameSpice closely.
  Vec6 GcrfToItrf(Real t_tdb, const Vec6& rv_gcrf) {
    Mat3 R_po = RotPolarMotion(t_tdb);
    Mat3 R_pn = RotPrecessionNutation(t_tdb);
    Mat3 R_s = RotSideralMotion(t_tdb);
    Mat3 R_s_dot = RotSideralMotionDot(t_tdb);

    Mat3 R_GcrfToItrf = R_po * R_s * R_pn;
    Mat3 R_GcrfToItrf_dot = R_po * R_s_dot * R_pn;

    Vec3 r_gcrf = rv_gcrf.head(3);
    Vec3 v_gcrf = rv_gcrf.tail(3);

    Vec3 r_itrf = R_GcrfToItrf * r_gcrf;
    Vec3 v_itrf = R_GcrfToItrf * v_gcrf + R_GcrfToItrf_dot * r_gcrf;

    Vec6 rv_itrf;
    rv_itrf << r_itrf, v_itrf;
    return rv_itrf;
  }

  Vec3 GcrfToItrf(Real t_tdb, const Vec3& r_gcrf) {
    Mat3 R_po = RotPolarMotion(t_tdb);
    Mat3 R_pn = RotPrecessionNutation(t_tdb);
    Mat3 R_s = RotSideralMotion(t_tdb);
    Mat3 R_GcrfToItrf = R_po * R_s * R_pn;

    Vec3 r_itrf = R_GcrfToItrf * r_gcrf;
    return r_itrf;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 37
  ///
  /// See the note on GcrfToItrf regarding InitFrameConversionFromSpice.
  Vec6 ItrfToGcrf(Real t_tdb, const Vec6& rv_itrf) {
    Mat3 R_po = RotPolarMotion(t_tdb);
    Mat3 R_pn = RotPrecessionNutation(t_tdb);
    Mat3 R_s = RotSideralMotion(t_tdb);
    Mat3 R_s_dot = RotSideralMotionDot(t_tdb);

    Mat3 R_GcrfToItrf = R_po * R_s * R_pn;
    Mat3 R_GcrfToItrf_dot = R_po * R_s_dot * R_pn;

    Vec3 r_itrf = rv_itrf.head(3);
    Vec3 v_itrf = rv_itrf.tail(3);

    Vec3 r_gcrf = R_GcrfToItrf.transpose() * r_itrf;
    Vec3 v_gcrf = R_GcrfToItrf.transpose() * (v_itrf - R_GcrfToItrf_dot * r_gcrf);

    Vec6 rv_gcrf;
    rv_gcrf << r_gcrf, v_gcrf;
    return rv_gcrf;
  }

  Vec3 ItrfToGcrf(Real t_tdb, const Vec3& r_itrf) {
    Mat3 R_po = RotPolarMotion(t_tdb);
    Mat3 R_pn = RotPrecessionNutation(t_tdb);
    Mat3 R_s = RotSideralMotion(t_tdb);
    Mat3 R_GcrfToItrf = R_po * R_s * R_pn;

    Vec3 r_gcrf = R_GcrfToItrf.transpose() * r_itrf;
    return r_gcrf;
  }

  /// @note
  /// Also known as Earth Frame Bias Matrix
  /// Astrodynamics Convention & Modeling Reference, Version 1.1, Page 38
  Mat3d RotGcrfToEme() {
    const double da = FRAME_BIAS_DALPHA0;
    const double xi = FRAME_BIAS_XI0;
    const double eta = FRAME_BIAS_ETA0;
    Mat3d B_e = RotX(-eta) * RotY(xi) * RotZ(da);
    return B_e;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 38
  Mat3d RotGcrfToEmeFirstOrder() {
    const double da = FRAME_BIAS_DALPHA0;
    const double xi = FRAME_BIAS_XI0;
    const double eta = FRAME_BIAS_ETA0;
    Mat3d B_e{{1, da, -xi}, {-da, 1, -eta}, {xi, eta, 1}};
    return B_e;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 38
  Mat3d RotGcrfToEmeSecondOrder() {
    const double da = FRAME_BIAS_DALPHA0;
    const double xi = FRAME_BIAS_XI0;
    const double eta = FRAME_BIAS_ETA0;
    Mat3d B_e{{1. - 0.5 * (da * da + xi * xi), da, -xi},
              {-da - eta * xi, 1. - 0.5 * (da * da + eta * eta), -eta},
              {xi - eta * da, eta + xi * da, 1. - 0.5 * (eta * eta + xi * xi)}};
    return B_e;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 39
  Vec6 GcrfToEme(const Vec6& rv_gcrf) {
    Mat3d B_e = RotGcrfToEme();
    Vec3 r_gcrf = rv_gcrf.head(3);
    Vec3 v_gcrf = rv_gcrf.tail(3);

    Vec3 r_eme = B_e * r_gcrf;
    Vec3 v_eme = B_e * v_gcrf;

    Vec6 rv_eme;
    rv_eme << r_eme, v_eme;
    return rv_eme;
  }
  Vec3 GcrfToEme(const Vec3& r_gcrf) {
    Mat3d B_e = RotGcrfToEme();
    Vec3 r_eme = B_e * r_gcrf;
    return r_eme;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 39
  Vec6 EmeToGcrf(const Vec6& rv_eme) {
    Mat3d B_e = RotGcrfToEme();
    Vec3 r_eme = rv_eme.head(3);
    Vec3 v_eme = rv_eme.tail(3);

    Mat3d B_e_inv = B_e.transpose();
    Vec3 r_gcrf = B_e_inv * r_eme;
    Vec3 v_gcrf = B_e_inv * v_eme;

    Vec6 rv_gcrf;
    rv_gcrf << r_gcrf, v_gcrf;
    return rv_gcrf;
  }
  Vec3 EmeToGcrf(const Vec3& r_eme) {
    Mat3d B_e = RotGcrfToEme();
    Mat3d B_e_inv = B_e.transpose();
    Vec3 r_gcrf = B_e_inv * r_eme;
    return r_gcrf;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 39
  Vec6 GcrfToIcrf(Real t_tdb, const Vec6& rv_gcrf) {
    Vec6 rv_ssb2earth = GetBodyPosVel(t_tdb, BodyId::SSB, BodyId::EARTH, Frame::GCRF);
    Vec6 rv_icrf = rv_gcrf + rv_ssb2earth;
    return rv_icrf;
  }
  Vec3 GcrfToIcrf(Real t_tdb, const Vec3& r_gcrf) {
    Vec3 r_ssb2earth = GetBodyPos(t_tdb, BodyId::SSB, BodyId::EARTH, Frame::GCRF);
    Vec3 r_icrf = r_gcrf + r_ssb2earth;
    return r_icrf;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 40
  Vec6 IcrfToGcrf(Real t_tdb, const Vec6& rv_icrf) {
    Vec6 rv_ssb2earth = GetBodyPosVel(t_tdb, BodyId::SSB, BodyId::EARTH, Frame::GCRF);
    Vec6 rv_gcrf = rv_icrf - rv_ssb2earth;
    return rv_gcrf;
  }
  Vec3 IcrfToGcrf(Real t_tdb, const Vec3& r_icrf) {
    Vec3 r_ssb2earth = GetBodyPos(t_tdb, BodyId::SSB, BodyId::EARTH, Frame::GCRF);
    Vec3 r_gcrf = r_icrf - r_ssb2earth;
    return r_gcrf;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 40
  Vec6 GcrfToMoonCi(Real t_tdb, const Vec6& rv_gcrf) {
    Vec6 rv_earth2moon = GetBodyPosVel(t_tdb, BodyId::EARTH, BodyId::MOON, Frame::GCRF);
    Vec6 rv_mi = rv_gcrf - rv_earth2moon;
    return rv_mi;
  }
  Vec3 GcrfToMoonCi(Real t_tdb, const Vec3& r_gcrf) {
    Vec3 r_earth2moon = GetBodyPos(t_tdb, BodyId::EARTH, BodyId::MOON, Frame::GCRF);
    Vec3 r_mi = r_gcrf - r_earth2moon;
    return r_mi;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 40
  Vec6 MoonCiToGcrf(Real t_tdb, const Vec6& rv_mi) {
    Vec6 rv_earth2moon = GetBodyPosVel(t_tdb, BodyId::EARTH, BodyId::MOON, Frame::GCRF);
    Vec6 rv_gcrf = rv_mi + rv_earth2moon;
    return rv_gcrf;
  }
  Vec3 MoonCiToGcrf(Real t_tdb, const Vec3& r_mi) {
    Vec3 r_earth2moon = GetBodyPos(t_tdb, BodyId::EARTH, BodyId::MOON, Frame::GCRF);
    Vec3 r_gcrf = r_mi + r_earth2moon;
    return r_gcrf;
  }

  struct CacheRotMoonCiToPa {
    Real t_tdb = NAN;
    Mat3 R_mi2pa;
    Mat3 R_mi2pa_dot;
    bool compute_vel;
  };
  static CacheRotMoonCiToPa cache_rot_moon_ci2pa;

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 42
  ///
  /// If InitFrameConversionFromSpice has fitted a lunar-orientation model
  /// covering `t_tdb`, the SPICE-fitted MOON_CI->MOON_PA rotation (and its
  /// exact analytic derivative) is used directly in place of the
  /// DE-Chebyshev "lunar mantle libration" angle reconstruction (see
  /// GetLunarMantleData), so that this matches
  /// lupnt::spice::ConvertFrameSpice closely.
  Mat3 RotMoonCiToPa(Real t_tdb, Mat3* R_mi2pa_dot) {
    bool compute_vel = (R_mi2pa_dot != nullptr);

    {
      Mat3 R_mi2pa_fit, R_mi2pa_dot_fit;
      if (TryGetFittedMoonCiToPa(t_tdb, &R_mi2pa_fit, compute_vel ? &R_mi2pa_dot_fit : nullptr)) {
        if (compute_vel) *R_mi2pa_dot = R_mi2pa_dot_fit;
        return R_mi2pa_fit;
      }
    }

    // Check cache
    if (abs(t_tdb - cache_rot_moon_ci2pa.t_tdb) < EPS
        && cache_rot_moon_ci2pa.compute_vel == compute_vel) {
      Mat3 R_mi2pa;
      bool cached = false;
#pragma omp critical
      {
        if (abs(t_tdb - cache_rot_moon_ci2pa.t_tdb) < EPS
            && cache_rot_moon_ci2pa.compute_vel == compute_vel) {
          if (compute_vel) *R_mi2pa_dot = cache_rot_moon_ci2pa.R_mi2pa_dot;
          R_mi2pa = cache_rot_moon_ci2pa.R_mi2pa;
          cached = true;
        }
      }
      if (cached) return R_mi2pa;
    }

    Vec6 lunar_mantle = GetLunarMantleData(t_tdb, compute_vel);
    auto [phi, theta, psi, phi_dot, theta_dot, psi_dot] = Unpack(lunar_mantle);

    Mat3 R_z_phi = RotZ(phi), R_x_theta = RotX(theta), R_z_psi = RotZ(psi);
    Mat3 R_mi2pa = R_z_psi * R_x_theta * R_z_phi;

    if (compute_vel) {
      Real spsi = sin(psi);
      Real cpsi = cos(psi);
      Mat3 mat = RotZdot(psi, psi_dot) * R_x_theta * R_z_phi
                 + R_z_psi * RotXdot(theta, theta_dot) * R_z_phi
                 + R_z_psi * R_x_theta * RotZdot(phi, phi_dot);
      *R_mi2pa_dot = mat;
    }

    // Update cache
#pragma omp critical
    {
      if (t_tdb > cache_rot_moon_ci2pa.t_tdb + EPS) {
        cache_rot_moon_ci2pa.t_tdb = t_tdb;
        cache_rot_moon_ci2pa.R_mi2pa = R_mi2pa;
        cache_rot_moon_ci2pa.compute_vel = compute_vel;
        if (compute_vel) cache_rot_moon_ci2pa.R_mi2pa_dot = *R_mi2pa_dot;
      }
    }
    return R_mi2pa;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 42
  Vec6 MoonCiToPa(Real t_tdb, const Vec6& rv_mi) {
    Mat3 R_mi2pa_dot;
    Mat3 R_mi2pa = RotMoonCiToPa(t_tdb, &R_mi2pa_dot);

    Vec3 r_mi = rv_mi.head(3);
    Vec3 v_mi = rv_mi.tail(3);

    Vec3 r_pa = R_mi2pa * r_mi;
    Vec3 v_pa = R_mi2pa * v_mi + R_mi2pa_dot * r_mi;

    Vec6 rv_pa;
    rv_pa << r_pa, v_pa;
    return rv_pa;
  }
  Vec3 MoonCiToPa(Real t_tdb, const Vec3& r_mi) {
    Mat3 R_mi2pa = RotMoonCiToPa(t_tdb);
    Vec3 r_pa = R_mi2pa * r_mi;
    return r_pa;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 42
  Vec6 MoonPaToCi(Real t_tdb, const Vec6& rv_pa) {
    Mat3 R_mi2pa_dot;
    Mat3 R_mi2pa = RotMoonCiToPa(t_tdb, &R_mi2pa_dot);

    Vec3 r_pa = rv_pa.head(3);
    Vec3 v_pa = rv_pa.tail(3);

    Vec3 r_mi = R_mi2pa.transpose() * r_pa;
    Vec3 v_mi = R_mi2pa.transpose() * (v_pa - R_mi2pa_dot * r_mi);
    Vec6 rv_mi;
    rv_mi << r_mi, v_mi;
    return rv_mi;
  }
  Vec3 MoonPaToCi(Real t_tdb, const Vec3& r_pa) {
    Mat3 R_mi2pa = RotMoonCiToPa(t_tdb);
    Vec3 r_mi = R_mi2pa.transpose() * r_pa;
    return r_mi;
  }

  Mat3d RotMoonPaToMe() {
    Mat3d B_moon
        = RotX(-0.2785 * RAD_ARCSEC) * RotY(-78.6944 * RAD_ARCSEC) * RotZ(-67.8526 * RAD_ARCSEC);
    return B_moon;
  }

  // ===========================================================================
  // Solar-system planet frames (generic IAU body orientation)
  // ===========================================================================
  namespace {
    // IAU linear orientation model, per body. Pole right ascension/declination
    // vary linearly with T (Julian centuries past J2000, TT); the prime meridian
    // W varies linearly with d (days past J2000, TT). Values from the IAU Working
    // Group on Cartographic Coordinates and Rotational Elements (2009/2015);
    // small periodic terms (e.g. for Mars/Neptune) are omitted -- the linear
    // model is accurate to well under a degree over multi-decade spans.
    struct IauOrientation {
      double ra0, ra0_T;    // pole RA:  ra0  + ra0_T  * T   [deg]
      double dec0, dec0_T;  // pole Dec: dec0 + dec0_T * T   [deg]
      double w0, w_dot;     // prime meridian: w0 + w_dot * d [deg]
    };

    const std::map<BodyId, IauOrientation>& IauOrientationTable() {
      static const std::map<BodyId, IauOrientation> kTable = {
          {BodyId::MERCURY, {281.0103, -0.0328, 61.4155, -0.0049, 329.5988, 6.1385108}},
          {BodyId::VENUS, {272.76, 0.0, 67.16, 0.0, 160.20, -1.4813688}},
          {BodyId::MARS, {317.68143, -0.1061, 52.88650, -0.0609, 176.630, 350.89198226}},
          {BodyId::JUPITER, {268.056595, -0.006499, 64.495303, 0.002413, 284.95, 870.5360000}},
          {BodyId::SATURN, {40.589, -0.036, 83.537, -0.004, 38.90, 810.7939024}},
          {BodyId::URANUS, {257.311, 0.0, -15.175, 0.0, 203.81, -501.1600928}},
          {BodyId::NEPTUNE, {299.36, 0.0, 43.46, 0.0, 249.978, 541.1397757}},
      };
      return kTable;
    }
  }  // namespace

  bool HasIauOrientation(BodyId body) { return IauOrientationTable().count(body) > 0; }

  bool IsPlanetFixedFrame(Frame frame) {
    switch (frame) {
      case Frame::MERCURY_FIXED:
      case Frame::VENUS_FIXED:
      case Frame::MARS_FIXED:
      case Frame::JUPITER_FIXED:
      case Frame::SATURN_FIXED:
      case Frame::URANUS_FIXED:
      case Frame::NEPTUNE_FIXED: return true;
      default: return false;
    }
  }

  bool IsPlanetCiFrame(Frame frame) {
    switch (frame) {
      case Frame::MERCURY_CI:
      case Frame::VENUS_CI:
      case Frame::MARS_CI:
      case Frame::JUPITER_CI:
      case Frame::SATURN_CI:
      case Frame::URANUS_CI:
      case Frame::NEPTUNE_CI: return true;
      default: return false;
    }
  }

  bool IsPlanetFrame(Frame frame) { return IsPlanetFixedFrame(frame) || IsPlanetCiFrame(frame); }

  Mat3 RotBodyCiToFixed(Real t_tdb, BodyId body, Mat3* R_dot) {
    auto it = IauOrientationTable().find(body);
    LUPNT_CHECK(it != IauOrientationTable().end(),
                fmt::format("No IAU orientation model for body {}", static_cast<int>(body)),
                "FrameConverter");
    const IauOrientation& o = it->second;

    Real t_tt = ConvertTime(t_tdb, Time::TDB, Time::TT);
    Real d = TimeToJd(t_tt) - JD_J2000_TT;  // days past J2000 (TT)
    Real T = d / 36525.0;                   // Julian centuries past J2000 (TT)

    Real ra = (o.ra0 + o.ra0_T * T) * RAD;
    Real dec = (o.dec0 + o.dec0_T * T) * RAD;
    Real W = (o.w0 + o.w_dot * d) * RAD;

    // ICRF -> body-fixed: R = Rz(W) * Rx(pi/2 - dec) * Rz(pi/2 + ra)
    Mat3 Rx_pole = RotX(M_PI / 2 - dec);
    Mat3 Rz_node = RotZ(M_PI / 2 + ra);
    Mat3 R = RotZ(W) * Rx_pole * Rz_node;
    if (R_dot != nullptr) {
      // Only the prime-meridian spin contributes appreciably (pole rates are
      // ~1e-2 deg/century -> negligible for velocity rotation).
      Real w_dot = o.w_dot * RAD / SECS_DAY;  // [rad/s]
      *R_dot = RotZdot(W, w_dot) * Rx_pole * Rz_node;
    }
    return R;
  }

  Vec6 PlanetCiToIcrf(Real t_tdb, const Vec6& rv_ci, BodyId body) {
    return rv_ci + GetBodyPosVel(t_tdb, body, Frame::ICRF);
  }
  Vec3 PlanetCiToIcrf(Real t_tdb, const Vec3& r_ci, BodyId body) {
    return r_ci + GetBodyPos(t_tdb, body, Frame::ICRF);
  }
  Vec6 IcrfToPlanetCi(Real t_tdb, const Vec6& rv_icrf, BodyId body) {
    return rv_icrf - GetBodyPosVel(t_tdb, body, Frame::ICRF);
  }
  Vec3 IcrfToPlanetCi(Real t_tdb, const Vec3& r_icrf, BodyId body) {
    return r_icrf - GetBodyPos(t_tdb, body, Frame::ICRF);
  }

  Vec6 BodyCiToFixed(Real t_tdb, const Vec6& rv_ci, BodyId body) {
    Mat3 R_dot;
    Mat3 R = RotBodyCiToFixed(t_tdb, body, &R_dot);
    Vec3 r_ci = rv_ci.head(3);
    Vec3 v_ci = rv_ci.tail(3);
    Vec6 rv_fixed;
    rv_fixed << R * r_ci, R * v_ci + R_dot * r_ci;
    return rv_fixed;
  }
  Vec3 BodyCiToFixed(Real t_tdb, const Vec3& r_ci, BodyId body) {
    return RotBodyCiToFixed(t_tdb, body) * r_ci;
  }
  Vec6 BodyFixedToCi(Real t_tdb, const Vec6& rv_fixed, BodyId body) {
    Mat3 R_dot;
    Mat3 R = RotBodyCiToFixed(t_tdb, body, &R_dot);
    Vec3 r_fixed = rv_fixed.head(3);
    Vec3 v_fixed = rv_fixed.tail(3);
    Vec3 r_ci = R.transpose() * r_fixed;
    Vec3 v_ci = R.transpose() * (v_fixed - R_dot * r_ci);
    Vec6 rv_ci;
    rv_ci << r_ci, v_ci;
    return rv_ci;
  }
  Vec3 BodyFixedToCi(Real t_tdb, const Vec3& r_fixed, BodyId body) {
    return RotBodyCiToFixed(t_tdb, body).transpose() * r_fixed;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 43
  Vec6 MoonPaToMe(const Vec6& rv_pa) {
    Mat3d B_moon = RotMoonPaToMe();
    Vec3 r_pa = rv_pa.head(3);
    Vec3 v_pa = rv_pa.tail(3);

    Vec3 r_me = B_moon * r_pa;
    Vec3 v_me = B_moon * v_pa;

    Vec6 rv_me;
    rv_me << r_me, v_me;
    return rv_me;
  }
  Vec3 MoonPaToMe(const Vec3& r_pa) {
    Mat3d B_moon = RotMoonPaToMe();
    Vec3 r_me = B_moon * r_pa;
    return r_me;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 43
  Vec6 MoonMeToPa(const Vec6& rv_me) {
    Mat3d B_moon = RotMoonPaToMe();
    Vec3 r_me = rv_me.head(3);
    Vec3 v_me = rv_me.tail(3);

    Mat3d B_moon_inv = B_moon.transpose();
    Vec3 r_pa = B_moon_inv * r_me;
    Vec3 v_pa = B_moon_inv * v_me;

    Vec6 rv_pa;
    rv_pa << r_pa, v_pa;
    return rv_pa;
  }
  Vec3 MoonMeToPa(const Vec3& r_me) {
    Mat3d B_moon = RotMoonPaToMe();
    Mat3d B_moon_inv = B_moon.transpose();
    Vec3 r_pa = B_moon_inv * r_me;
    return r_pa;
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 48
  Vec6 GcrfToEmr(Real t_tdb, const Vec6& rv_gcrf) {
    Vec6 rv_earth2emb = GetBodyPosVel(t_tdb, BodyId::EARTH, BodyId::EMB, Frame::GCRF);
    Vec6 rv_emr = InertialToSynodic(rv_earth2emb, rv_gcrf);
    return rv_emr;
  }
  Vec3 GcrfToEmr(Real t_tdb, const Vec3& r_gcrf) {
    Vec6 rv_earth2emb = GetBodyPosVel(t_tdb, BodyId::EARTH, BodyId::EMB, Frame::GCRF);
    Vec6 rv_gcrf;
    rv_gcrf << r_gcrf, Vec3::Zero();
    Vec6 rv_emr = InertialToSynodic(rv_earth2emb, rv_gcrf);
    return rv_emr.head(3);
  }

  /// @note Astrodynamics Convention & Modeling Reference, Version 1.1, Page 48
  Vec6 EmrToGcrf(Real t_tdb, const Vec6& rv_emr) {
    Vec6 rv_earth2emb = GetBodyPosVel(t_tdb, BodyId::EARTH, BodyId::EMB, Frame::GCRF);
    Vec6 rv_gcrf = SynodicToInertial(rv_earth2emb, rv_emr);
    return rv_gcrf;
  }
  Vec3 EmrToGcrf(Real t_tdb, const Vec3& r_emr) {
    Vec6 rv_earth2emb = GetBodyPosVel(t_tdb, BodyId::EARTH, BodyId::EMB, Frame::GCRF);
    Vec6 rv_emr;
    rv_emr << r_emr, Vec3::Zero();
    Vec6 rv_gcrf = SynodicToInertial(rv_earth2emb, rv_emr);
    return rv_gcrf.head(3);
  }

  /// @note
  /// T. A. Ely, ‘Stable Constellations of Frozen Elliptical Inclined Lunar
  /// Orbits’, J of Astronaut Sci, vol. 53, no. 3, pp. 301–316, Sep. 2005,
  /// doi: 10.1007/BF03546355.
  Mat3 RotMoonOpToCi(Real t_tdb) {
    Vec6 rv_m2e = GetBodyPosVel(t_tdb, BodyId::MOON, BodyId::EARTH, Frame::GCRF);
    Vec3 r = rv_m2e.head(3);
    Vec3 v = rv_m2e.tail(3);

    Vec3 z_op = r.cross(v).normalized();
    Vec3 i_pole_me = Vec3::UnitZ();
    Vec3 i_pole = ConvertFrame(t_tdb, i_pole_me, Frame::MOON_ME, Frame::MOON_CI);
    Vec3 x_op = i_pole.cross(z_op).normalized();
    Vec3 y_op = z_op.cross(x_op).normalized();

    // R = [x_OP, y_OP, z_OP]_CI (expressed in CI)
    Mat3 R_op2ci;
    R_op2ci << x_op, y_op, z_op;
    return R_op2ci;
  }

  /// @brief
  /// @param t_tdb
  /// @param rv_ci
  /// @return Vec6
  /// @note
  /// T. A. Ely, ‘Stable Constellations of Frozen Elliptical Inclined Lunar
  /// Orbits’, J of Astronaut Sci, vol. 53, no. 3, pp. 301–316, Sep. 2005,
  /// doi: 10.1007/BF03546355.
  Vec6 MoonCiToOp(Real t_tdb, const Vec6& rv_ci) {
    Mat3 R_ci2op = RotMoonOpToCi(t_tdb).transpose();
    Vec3 r_op = R_ci2op * rv_ci.head(3);
    Vec3 v_op = R_ci2op * rv_ci.tail(3);
    Vec6 rv_op;
    rv_op << r_op, v_op;
    return rv_op;
  }
  Vec3 MoonCiToOp(Real t_tdb, const Vec3& r_ci) {
    Mat3 R_ci2op = RotMoonOpToCi(t_tdb).transpose();
    Vec3 r_op = R_ci2op * r_ci;
    return r_op;
  }

  /// @brief
  /// @param rv_op
  /// @return Vec6
  /// @note
  /// T. A. Ely, ‘Stable Constellations of Frozen Elliptical Inclined Lunar
  /// Orbits’, J of Astronaut Sci, vol. 53, no. 3, pp. 301–316, Sep. 2005,
  /// doi: 10.1007/BF03546355.
  Vec6 MoonOpToCi(Real t_tdb, const Vec6& rv_op) {
    Mat3 R_op2ci = RotMoonOpToCi(t_tdb);
    Vec3 r_ci = R_op2ci * rv_op.head(3);
    Vec3 v_ci = R_op2ci * rv_op.tail(3);
    Vec6 rv_ci;
    rv_ci << r_ci, v_ci;
    return rv_ci;
  }
  Vec3 MoonOpToCi(Real t_tdb, const Vec3& r_op) {
    Mat3 R_op2ci = RotMoonOpToCi(t_tdb);
    Vec3 r_ci = R_op2ci * r_op;
    return r_ci;
  }

}  // namespace lupnt

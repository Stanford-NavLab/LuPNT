#pragma once

#include "lupnt/core/constants.h"
#include "lupnt/interfaces/spice_cheby.h"

namespace lupnt {

  /// @brief Converts a state vector from ITRF to GCRF (ECEF to ECI) at epoch `t_tdb`.
  ///
  /// Used throughout the simulator wherever Earth-fixed states (ground-station
  /// positions, ECEF-frame measurements, ITRF-referenced SPICE/EOP data) must be
  /// expressed in the inertial frame -- e.g. by `ConvertFrame` when propagating
  /// dynamics or computing line-of-sight vectors in GCRF.
  ///
  /// @param t_tdb    Epoch [s, TDB since J2000]
  /// @param rv_itrf  Position+velocity in ITRF (ECEF) [m, m/s]
  /// @return         Position+velocity in GCRF (ECI) [m, m/s]
  Vec6 ItrfToGcrf(Real t_tdb, const Vec6& rv_itrf);
  /// @brief Vec3 (position-only) overload of ItrfToGcrf().
  Vec3 ItrfToGcrf(Real t_tdb, const Vec3& rv_itrf);

  /// @brief Converts a state vector from GCRF to ITRF (ECI to ECEF) at epoch `t_tdb`.
  ///
  /// Used throughout the simulator wherever inertial states (propagated dynamics,
  /// SPICE ephemerides) must be expressed in the Earth-fixed frame -- e.g. ground
  /// station visibility, ECEF-frame measurements, and ITRF-referenced outputs.
  ///
  /// @param t_tdb    Epoch [s, TDB since J2000]
  /// @param rv_gcrf  Position+velocity in GCRF (ECI) [m, m/s]
  /// @return         Position+velocity in ITRF (ECEF) [m, m/s]
  Vec6 GcrfToItrf(Real t_tdb, const Vec6& rv_gcrf);
  /// @brief Vec3 (position-only) overload of GcrfToItrf().
  Vec3 GcrfToItrf(Real t_tdb, const Vec3& rv_gcrf);

  /// @brief Converts a state vector from GCRF to EME2000 (J2000 mean equator/equinox)
  /// by applying the constant IAU frame-bias rotation.
  ///
  /// EME2000 (a.k.a. ECI/J2000) is the conventional inertial frame used by many
  /// legacy ephemerides and dynamics models; `ConvertFrame` calls this when
  /// converting a GCRF state to/from `Frame::EME`/`Frame::ECI`. Unlike most other
  /// conversions in this header, the GCRF<->EME rotation is a fixed (time-independent)
  /// bias, so no epoch argument is needed.
  ///
  /// @param rv_gcrf  Position+velocity in GCRF [m, m/s]
  /// @return         Position+velocity in EME2000 (ECI) [m, m/s]
  Vec6 GcrfToEme(const Vec6& rv_gcrf);
  /// @brief Vec3 (position-only) overload of GcrfToEme().
  Vec3 GcrfToEme(const Vec3& rv_gcrf);
  /// @brief Converts a state vector from EME2000 (J2000 mean equator/equinox) to GCRF
  /// by applying the inverse of the constant IAU frame-bias rotation.
  ///
  /// Used by `ConvertFrame` when converting an EME/ECI-frame state into GCRF (the
  /// hub frame from which all other Earth/Moon frame conversions are reached).
  ///
  /// @param rv_eme  Position+velocity in EME2000 (ECI) [m, m/s]
  /// @return        Position+velocity in GCRF [m, m/s]
  Vec6 EmeToGcrf(const Vec6& rv_eme);
  /// @brief Vec3 (position-only) overload of EmeToGcrf().
  Vec3 EmeToGcrf(const Vec3& rv_eme);

  /// @brief Converts a state vector from GCRF (Earth-centered) to ICRF
  /// (solar-system-barycenter-centered), by adding the Earth's barycentric
  /// position/velocity at `t_tdb`.
  ///
  /// Used by `ConvertFrame` to bridge between the Earth-centered GCRF hub frame and
  /// the barycentric ICRF frame, e.g. when an agent's state must be expressed
  /// relative to the solar-system barycenter for high-fidelity N-body dynamics or
  /// SPICE-based ephemeris comparisons.
  ///
  /// @param t_tdb    Epoch [s, TDB since J2000]
  /// @param rv_gcrf  Position+velocity in GCRF (Earth-centered) [m, m/s]
  /// @return         Position+velocity in ICRF (solar-system-barycenter-centered) [m, m/s]
  Vec6 GcrfToIcrf(Real t_tdb, const Vec6& rv_gcrf);
  /// @brief Vec3 (position-only) overload of GcrfToIcrf().
  Vec3 GcrfToIcrf(Real t_tdb, const Vec3& rv_gcrf);
  /// @brief Converts a state vector from ICRF (solar-system-barycenter-centered) to
  /// GCRF (Earth-centered), by subtracting the Earth's barycentric position/velocity
  /// at `t_tdb`.
  ///
  /// Used by `ConvertFrame` when converting a barycentric ICRF state into the
  /// Earth-centered GCRF hub frame.
  ///
  /// @param t_tdb    Epoch [s, TDB since J2000]
  /// @param rv_icrf  Position+velocity in ICRF (solar-system-barycenter-centered) [m, m/s]
  /// @return         Position+velocity in GCRF (Earth-centered) [m, m/s]
  Vec6 IcrfToGcrf(Real t_tdb, const Vec6& rv_icrf);
  /// @brief Vec3 (position-only) overload of IcrfToGcrf().
  Vec3 IcrfToGcrf(Real t_tdb, const Vec3& rv_icrf);

  /// @brief Converts a state vector from GCRF (Earth-centered inertial) to MoonCi
  /// (Moon-centered inertial, axes aligned with ICRF/GCRF), by subtracting the
  /// Earth-to-Moon position/velocity at `t_tdb`.
  ///
  /// Used by `ConvertFrame` as the entry point into the Moon-centered frame family
  /// (MoonCi/MoonPa/MoonMe/MoonOp) -- e.g. when expressing a lunar orbiter's or
  /// surface asset's inertial state relative to the Moon instead of the Earth.
  ///
  /// @param t_tdb    Epoch [s, TDB since J2000]
  /// @param rv_gcrf  Position+velocity in GCRF (Earth-centered) [m, m/s]
  /// @return         Position+velocity in MoonCi (Moon-centered inertial) [m, m/s]
  Vec6 GcrfToMoonCi(Real t_tdb, const Vec6& rv_gcrf);
  /// @brief Vec3 (position-only) overload of GcrfToMoonCi().
  Vec3 GcrfToMoonCi(Real t_tdb, const Vec3& rv_gcrf);
  /// @brief Converts a state vector from MoonCi (Moon-centered inertial) to GCRF
  /// (Earth-centered inertial), by adding the Earth-to-Moon position/velocity at
  /// `t_tdb`.
  ///
  /// Used by `ConvertFrame` as the exit point from the Moon-centered frame family
  /// back to the GCRF hub frame.
  ///
  /// @param t_tdb  Epoch [s, TDB since J2000]
  /// @param rv_mi  Position+velocity in MoonCi (Moon-centered inertial) [m, m/s]
  /// @return       Position+velocity in GCRF (Earth-centered) [m, m/s]
  Vec6 MoonCiToGcrf(Real t_tdb, const Vec6& rv_mi);
  /// @brief Vec3 (position-only) overload of MoonCiToGcrf().
  Vec3 MoonCiToGcrf(Real t_tdb, const Vec3& rv_mi);

  /// @brief Converts a state vector from MoonCi (Moon-centered inertial, ICRF-aligned
  /// axes) to MoonPa (Moon-fixed, principal-axis frame) at epoch `t_tdb`.
  ///
  /// Used by `ConvertFrame` to express lunar-orbiter/lander states in the
  /// body-fixed principal-axis frame -- e.g. for lunar gravity-field evaluation,
  /// surface feature lookups, and as the bridge to MoonMe (via MoonPaToMe).
  /// Internally applies RotMoonCiToPa (and its time derivative) to rotate both
  /// position and velocity.
  ///
  /// @param t_tdb  Epoch [s, TDB since J2000]
  /// @param rv_mi  Position+velocity in MoonCi [m, m/s]
  /// @return       Position+velocity in MoonPa (Moon-fixed principal axes) [m, m/s]
  Vec6 MoonCiToPa(Real t_tdb, const Vec6& rv_mi);
  /// @brief Vec3 (position-only) overload of MoonCiToPa().
  Vec3 MoonCiToPa(Real t_tdb, const Vec3& rv_mi);
  /// @brief Converts a state vector from MoonPa (Moon-fixed, principal-axis frame) to
  /// MoonCi (Moon-centered inertial, ICRF-aligned axes) at epoch `t_tdb`.
  ///
  /// Used by `ConvertFrame` as the inverse of MoonCiToPa, e.g. to bring a
  /// gravity-model or surface-relative lunar state back into the inertial frame for
  /// dynamics propagation. Internally applies the transpose of RotMoonCiToPa (and
  /// its time derivative) to rotate both position and velocity.
  ///
  /// @param t_tdb  Epoch [s, TDB since J2000]
  /// @param rv_pa  Position+velocity in MoonPa (Moon-fixed principal axes) [m, m/s]
  /// @return       Position+velocity in MoonCi [m, m/s]
  Vec6 MoonPaToCi(Real t_tdb, const Vec6& rv_pa);
  /// @brief Vec3 (position-only) overload of MoonPaToCi().
  Vec3 MoonPaToCi(Real t_tdb, const Vec3& rv_pa);

  /// @brief Converts a state vector from MoonPa (Moon-fixed principal axes) to MoonMe
  /// (Moon-fixed mean-Earth/polar axes) by applying the constant IAU
  /// principal-axis-to-mean-Earth bias rotation.
  ///
  /// Used by `ConvertFrame` (and directly after MoonCiToPa for MoonCi->MoonMe) to
  /// express lunar surface/body-fixed states in the mean-Earth/polar frame commonly
  /// used for lunar cartography and surface-station coordinates. Like GcrfToEme,
  /// this is a fixed (time-independent) bias rotation, so no epoch is needed.
  ///
  /// @param rv_pa  Position+velocity in MoonPa (Moon-fixed principal axes) [m, m/s]
  /// @return       Position+velocity in MoonMe (Moon-fixed mean-Earth/polar axes) [m, m/s]
  Vec6 MoonPaToMe(const Vec6& rv_pa);
  /// @brief Vec3 (position-only) overload of MoonPaToMe().
  Vec3 MoonPaToMe(const Vec3& rv_pa);
  /// @brief Converts a state vector from MoonMe (Moon-fixed mean-Earth/polar axes) to
  /// MoonPa (Moon-fixed principal axes) by applying the inverse of the constant IAU
  /// principal-axis-to-mean-Earth bias rotation.
  ///
  /// Used by `ConvertFrame` as the inverse of MoonPaToMe, to bring a
  /// mean-Earth/polar-frame state (e.g. a surface-station position defined in MoonMe)
  /// into the principal-axis frame before further conversion (e.g. to MoonCi via
  /// MoonPaToCi).
  ///
  /// @param rv_me  Position+velocity in MoonMe (Moon-fixed mean-Earth/polar axes) [m, m/s]
  /// @return       Position+velocity in MoonPa (Moon-fixed principal axes) [m, m/s]
  Vec6 MoonMeToPa(const Vec6& rv_me);
  /// @brief Vec3 (position-only) overload of MoonMeToPa().
  Vec3 MoonMeToPa(const Vec3& rv_me);

  /// @brief Converts a state vector from GCRF (Earth-centered inertial) to EMR
  /// (Earth-Moon Rotating frame), via InertialToSynodic relative to the
  /// instantaneous Earth-Moon-barycenter direction.
  ///
  /// Used by `ConvertFrame` to express an Earth-centered inertial state in the
  /// rotating Earth-Moon frame -- e.g. for visualizing/analyzing libration-point
  /// trajectories or cislunar transfer geometry relative to the Earth-Moon line.
  ///
  /// @param t_tdb    Epoch [s, TDB since J2000]
  /// @param rv_gcrf  Position+velocity in GCRF (Earth-centered inertial) [m, m/s]
  /// @return         Position+velocity in EMR (Earth-Moon Rotating frame) [m, m/s]
  Vec6 GcrfToEmr(Real t_tdb, const Vec6& rv_gcrf);
  /// @brief Vec3 (position-only) overload of GcrfToEmr(). Velocity is set to zero
  /// internally before the synodic transform, so only the rotated position is
  /// meaningful.
  Vec3 GcrfToEmr(Real t_tdb, const Vec3& rv_gcrf);
  /// @brief Converts a state vector from EMR (Earth-Moon Rotating frame) to GCRF
  /// (Earth-centered inertial), via SynodicToInertial relative to the instantaneous
  /// Earth-Moon-barycenter direction.
  ///
  /// Used by `ConvertFrame` as the inverse of GcrfToEmr, to bring an
  /// Earth-Moon-rotating-frame state (e.g. a target point defined in EMR) back into
  /// the GCRF hub frame for dynamics propagation.
  ///
  /// @param t_tdb   Epoch [s, TDB since J2000]
  /// @param rv_emr  Position+velocity in EMR (Earth-Moon Rotating frame) [m, m/s]
  /// @return        Position+velocity in GCRF (Earth-centered inertial) [m, m/s]
  Vec6 EmrToGcrf(Real t_tdb, const Vec6& rv_emr);
  /// @brief Vec3 (position-only) overload of EmrToGcrf(). Velocity is set to zero
  /// internally before the inertial transform, so only the rotated position is
  /// meaningful.
  Vec3 EmrToGcrf(Real t_tdb, const Vec3& rv_emr);

  /// @brief Converts a state vector from MoonCi (Moon-centered inertial) to MoonOp
  /// (Moon orbit-plane frame: z-axis along the Earth-Moon orbit-normal, x-axis toward
  /// the lunar north pole projected into the orbit plane) at epoch `t_tdb`.
  ///
  /// Used by `ConvertFrame` to express lunar-orbiter states in the frozen-orbit
  /// "orbit plane" frame of Ely (2005), e.g. for designing/analyzing stable frozen
  /// elliptical inclined lunar orbits. Internally applies the transpose of
  /// RotMoonOpToCi.
  ///
  /// @param t_tdb  Epoch [s, TDB since J2000]
  /// @param rv_ci  Position+velocity in MoonCi (Moon-centered inertial) [m, m/s]
  /// @return       Position+velocity in MoonOp (Moon orbit-plane frame) [m, m/s]
  Vec6 MoonCiToOp(Real t_tdb, const Vec6& rv_ci);
  /// @brief Vec3 (position-only) overload of MoonCiToOp().
  Vec3 MoonCiToOp(Real t_tdb, const Vec3& rv_ci);
  /// @brief Converts a state vector from MoonOp (Moon orbit-plane frame) to MoonCi
  /// (Moon-centered inertial) at epoch `t_tdb`.
  ///
  /// Used by `ConvertFrame` as the inverse of MoonCiToOp, to bring a frozen-orbit
  /// "orbit plane" state back into the Moon-centered inertial frame for dynamics
  /// propagation. Internally applies RotMoonOpToCi.
  ///
  /// @param t_tdb  Epoch [s, TDB since J2000]
  /// @param rv_op  Position+velocity in MoonOp (Moon orbit-plane frame) [m, m/s]
  /// @return       Position+velocity in MoonCi (Moon-centered inertial) [m, m/s]
  Vec6 MoonOpToCi(Real t_tdb, const Vec6& rv_op);
  /// @brief Vec3 (position-only) overload of MoonOpToCi().
  Vec3 MoonOpToCi(Real t_tdb, const Vec3& rv_op);

  // Rotations

  /// @brief Returns the precession-nutation rotation matrix R_pn at epoch `t_tdb`,
  /// which rotates a vector from GCRF (ICRF-aligned, mean-equator-of-J2000) into the
  /// true-equator-and-equinox-of-date frame (CIRS-like intermediate frame), i.e.
  /// `v_intermediate = R_pn * v_GCRF` (IAU 2006/2000A precession-nutation, computed
  /// from the SOFA (X, Y, s) CIP coordinates via GetIauSofaData).
  ///
  /// Used internally by GcrfToItrf/ItrfToGcrf (and by ComputeEopFromSpice, which
  /// inverts the same decomposition) as the first stage of the full
  /// GCRF<->ITRF rotation `R_GcrfToItrf = R_po * R_s * R_pn`.
  ///
  /// @param t_tdb  Epoch [s, TDB since J2000]
  /// @return       3x3 precession-nutation rotation matrix R_pn (dimensionless)
  Mat3 RotPrecessionNutation(Real t_tdb);

  /// @brief Returns the sidereal-motion (Earth-rotation) rotation matrix R_s at epoch
  /// `t_tdb`, which rotates a vector from the true-equator-and-equinox-of-date
  /// (intermediate) frame into the Earth-fixed (pseudo-ITRF, pre-polar-motion) frame
  /// by the Earth Rotation Angle (ERA) about the z-axis, i.e.
  /// `v_pseudo_itrf = R_s * v_intermediate`.
  ///
  /// Used internally by GcrfToItrf/ItrfToGcrf as the second stage of
  /// `R_GcrfToItrf = R_po * R_s * R_pn`. If InitFrameConversionFromSpice has fitted a
  /// SPICE-derived EOP model covering `t_tdb`, the fitted UT1-UTC is used in place of
  /// the loaded IERS EOP table when computing UT1 (and hence the ERA).
  ///
  /// @param t_tdb  Epoch [s, TDB since J2000]
  /// @return       3x3 sidereal-motion rotation matrix R_s (dimensionless)
  Mat3 RotSideralMotion(Real t_tdb);

  /// @brief Returns the exact time derivative dR_s/dt of the sidereal-motion
  /// rotation matrix (see RotSideralMotion) at epoch `t_tdb`, evaluated using the
  /// Earth rotation rate w_E (from the EOP length-of-day, or the analytic derivative
  /// of the SPICE-fitted UT1-UTC if InitFrameConversionFromSpice has fitted a model
  /// covering `t_tdb`).
  ///
  /// Used internally by GcrfToItrf/ItrfToGcrf to compute the velocity-transformation
  /// term `R_GcrfToItrf_dot = R_po * R_s_dot * R_pn`, which accounts for the
  /// Earth-fixed frame's rotation rate when converting velocities between GCRF and
  /// ITRF.
  ///
  /// @param t_tdb  Epoch [s, TDB since J2000]
  /// @return       3x3 matrix dR_s/dt [1/s] (time derivative of the sidereal-motion
  ///               rotation matrix)
  Mat3 RotSideralMotionDot(Real t_tdb);

  /// @brief Returns the polar-motion rotation matrix R_po at epoch `t_tdb`, which
  /// rotates a vector from the Earth-fixed pseudo-ITRF (post-sidereal-rotation)
  /// frame into the true ITRF frame by the IERS polar-motion angles (x_pole, y_pole)
  /// and the TIO locator s', i.e. `v_ITRF = R_po * v_pseudo_itrf`.
  ///
  /// Used internally by GcrfToItrf/ItrfToGcrf as the final stage of
  /// `R_GcrfToItrf = R_po * R_s * R_pn`. If InitFrameConversionFromSpice has fitted a
  /// SPICE-derived EOP model covering `t_tdb`, the fitted x_pole/y_pole are used in
  /// place of the loaded IERS EOP table.
  ///
  /// @param t_tdb  Epoch [s, TDB since J2000]
  /// @return       3x3 polar-motion rotation matrix R_po (dimensionless)
  Mat3 RotPolarMotion(Real t_tdb);

  /// @brief Returns the constant (time-independent) GCRF-to-EME2000 frame-bias
  /// rotation matrix B_e (a.k.a. the Earth Frame Bias Matrix), exact to all orders in
  /// the IAU 2006 frame-bias angles, such that `v_EME = B_e * v_GCRF`.
  ///
  /// Used internally by GcrfToEme/EmeToGcrf to convert between the GCRF (ICRF-aligned)
  /// and EME2000 (J2000 mean equator/equinox) inertial frames.
  ///
  /// @return  3x3 GCRF-to-EME2000 frame-bias rotation matrix B_e (dimensionless)
  Mat3d RotGcrfToEme();
  /// @brief First-order (small-angle, linearized) approximation of RotGcrfToEme(),
  /// accurate to first order in the IAU 2006 frame-bias angles.
  Mat3d RotGcrfToEmeFirstOrder();
  /// @brief Second-order approximation of RotGcrfToEme(), accurate to second order in
  /// the IAU 2006 frame-bias angles.
  Mat3d RotGcrfToEmeSecondOrder();

  /// @brief Returns the MoonCi-to-MoonPa rotation matrix R_mi2pa at epoch `t_tdb`,
  /// which rotates a vector from MoonCi (Moon-centered inertial, ICRF-aligned axes)
  /// into MoonPa (Moon-fixed principal-axis frame) via the ZXZ Euler-angle libration
  /// sequence `R_mi2pa = Rz(psi) * Rx(theta) * Rz(phi)`, i.e. `v_PA = R_mi2pa * v_CI`.
  ///
  /// Used internally by MoonCiToPa/MoonPaToCi to rotate (and, via `R_mi2pa_dot`,
  /// differentiate) position/velocity between the inertial and body-fixed lunar
  /// frames. The libration angles (phi, theta, psi) and their rates normally come
  /// from the DE Chebyshev "lunar mantle libration" data (GetLunarMantleData), but if
  /// InitFrameConversionFromSpice has fitted a SPICE-derived lunar-orientation model
  /// covering `t_tdb`, that fitted rotation (and derivative) is used instead so that
  /// this matches `lupnt::spice::ConvertFrameSpice` closely. Results are cached per
  /// epoch to avoid recomputation across repeated calls at the same `t_tdb`.
  ///
  /// @param t_tdb        Epoch [s, TDB since J2000]
  /// @param R_mi2pa_dot  [out, optional] If non-null, receives the exact time
  ///                      derivative dR_mi2pa/dt [1/s]
  /// @return              3x3 MoonCi-to-MoonPa rotation matrix R_mi2pa (dimensionless)
  Mat3 RotMoonCiToPa(Real t_tdb, Mat3* R_mi2pa_dot = nullptr);

  /// @brief Returns the constant (time-independent) MoonPa-to-MoonMe frame-bias
  /// rotation matrix B_moon, built from the IAU mean-Earth/polar-axis offset angles,
  /// such that `v_ME = B_moon * v_PA`.
  ///
  /// Used internally by MoonPaToMe/MoonMeToPa to convert between the Moon-fixed
  /// principal-axis (PA) and mean-Earth/polar-axis (ME) frames.
  ///
  /// @return  3x3 MoonPa-to-MoonMe frame-bias rotation matrix B_moon (dimensionless)
  Mat3d RotMoonPaToMe();

  /// @brief Returns the MoonOp-to-MoonCi rotation matrix R_op2ci at epoch `t_tdb`,
  /// whose columns [x_OP, y_OP, z_OP] (expressed in MoonCi) define the "orbit plane"
  /// frame of Ely (2005): z_OP along the Earth-Moon orbital angular-momentum
  /// direction, x_OP the projection of the lunar mean-Earth pole onto the plane
  /// normal to z_OP, and y_OP completing the right-handed triad, such that
  /// `v_CI = R_op2ci * v_OP`.
  ///
  /// Used internally by MoonOpToCi/MoonCiToOp (the latter via its transpose) to
  /// rotate position/velocity between the Moon-centered inertial frame and the
  /// frozen-orbit "orbit plane" frame.
  ///
  /// @param t_tdb  Epoch [s, TDB since J2000]
  /// @return       3x3 MoonOp-to-MoonCi rotation matrix R_op2ci (dimensionless)
  Mat3 RotMoonOpToCi(Real t_tdb);

  // ---------------------------------------------------------------------
  // SPICE-derived Earth orientation parameters (EOP)
  // ---------------------------------------------------------------------

  /// @brief Earth orientation parameters (polar motion + UT1-UTC) derived
  /// from SPICE's high-accuracy J2000<->ITRF93 rotation at a single epoch.
  struct SpiceEopParams {
    Real x_pole;   ///< Polar motion x_p [rad]
    Real y_pole;   ///< Polar motion y_p [rad]
    Real ut1_utc;  ///< UT1 - UTC [s]
  };

  /// @brief Derives the Earth orientation parameters (x_pole, y_pole,
  /// UT1-UTC) at `t_tdb` directly from SPICE's "J2000"->"ITRF93" rotation,
  /// combined with LuPNT's IAU 2006/2000A precession-nutation matrix
  /// (RotPrecessionNutation). This gives a SPICE-consistent alternative to
  /// the IERS-EOP-table-based values returned by GetEopData -- useful for
  /// epochs not covered by the loaded EOP table, or for comparing/fitting
  /// against SPICE's high-accuracy `earth_latest_high_prec.bpc` orientation
  /// (see InitFrameConversionFromSpice).
  SpiceEopParams ComputeEopFromSpice(Real t_tdb);

  // ---------------------------------------------------------------------
  // SPICE-fitted high-accuracy Earth/lunar orientation (initialization)
  // ---------------------------------------------------------------------

  /// @brief Calibrates LuPNT's native (non-SPICE) high-accuracy frame
  /// conversions -- GcrfToItrf / ItrfToGcrf (Earth orientation) and
  /// MoonCiToPa / MoonPaToCi (lunar orientation) -- to closely reproduce
  /// `lupnt::spice::ConvertFrameSpice` over a given simulation time window.
  ///
  /// LuPNT normally computes these rotations from purely-analytic models:
  /// IAU 2006/2000A precession-nutation (RotPrecessionNutation) plus
  /// EOP-based polar motion and sidereal rotation (RotPolarMotion /
  /// RotSideralMotion / RotSideralMotionDot) for Earth, where the EOP
  /// (x_pole, y_pole, UT1-UTC, LOD) come from the loaded IERS EOP table (see
  /// lupnt/interfaces/eop.h); and static, pre-extracted DE Chebyshev "lunar mantle
  /// libration" coefficients for the Moon (see GetLunarMantleData). Those
  /// can disagree with SPICE's high-accuracy binary-PCK-based orientation
  /// (which `lupnt::spice::LoadSpiceKernel` now downloads automatically --
  /// see `earth_latest_high_prec.bpc` / `moon_pa_de*.bpc` + `moon_de*.tf`)
  /// by a non-negligible amount -- in particular when the loaded EOP table
  /// does not cover the requested epoch.
  ///
  /// This function closes that gap for a known time window:
  ///  - Earth: samples ComputeEopFromSpice(t) -- which derives (x_pole,
  ///    y_pole, UT1-UTC) from SPICE's "J2000"<->"ITRF93" rotation -- at the
  ///    Chebyshev-Gauss nodes spanning [t_start_tdb, t_end_tdb], and fits
  ///    piecewise Chebyshev polynomials to those 3 EOP parameters (the same
  ///    Chebyshev representation used internally for the planetary/lunar
  ///    ephemeris data -- see spice_cheby.h / kernels.cc). Once initialized,
  ///    RotPolarMotion / RotSideralMotion / RotSideralMotionDot (and hence
  ///    GcrfToItrf / ItrfToGcrf) transparently use these SPICE-fitted EOP
  ///    values -- together with the *exact* analytic time derivative of
  ///    UT1-UTC, via cheby_eval_ad -- whenever the requested epoch falls
  ///    within the fitted window, instead of the loaded IERS EOP table.
  ///  - Moon: samples the SPICE "J2000"<->"IAU_MOON" rotation matrix at the
  ///    same Chebyshev-Gauss nodes, extracts ZXZ Euler angles (phi, theta,
  ///    psi) for the decomposition R = Rz(psi)*Rx(theta)*Rz(phi),
  ///    continuously unwraps phi and psi (so the fitted functions are smooth
  ///    despite the 2*pi-periodic branch cuts), and fits piecewise Chebyshev
  ///    polynomials to those 3 angles. At evaluation time the angles (and
  ///    their *exact* analytic time derivatives, via cheby_eval_ad) are used
  ///    to reconstruct R = Rz(psi)*Rx(theta)*Rz(phi), which is GUARANTEED
  ///    exactly orthonormal (to floating-point precision) because products
  ///    of elementary rotation matrices are always orthonormal --
  ///    regardless of the Chebyshev fit residual in the angles. Once
  ///    initialized, MoonCiToPa / MoonPaToCi automatically use the fitted
  ///    model whenever the requested epoch falls within the fitted window.
  ///
  /// Outside the fitted window(s), the original analytic computation is used
  /// as before, so this feature is purely additive/opt-in (nothing changes
  /// unless InitFrameConversionFromSpice is explicitly called).
  ///
  /// @param t_start_tdb     Start of the time window to fit [s, TDB since J2000]
  /// @param t_end_tdb       End of the time window to fit [s, TDB since J2000]
  /// @param segment_length_s  Length of each Chebyshev fit segment [s]
  ///                          (default: 1 day; shorter segments improve
  ///                          accuracy at the cost of memory/init time)
  /// @param num_coeffs      Number of Chebyshev coefficients (polynomial
  ///                        degree + 1) per fitted quantity and segment
  ///                        (default: 13)
  void InitFrameConversionFromSpice(Real t_start_tdb, Real t_end_tdb,
                                    Real segment_length_s = 86400.0, int num_coeffs = 13);

  /// @brief Clears any SPICE-fitted coefficients previously installed by
  /// InitFrameConversionFromSpice, reverting GcrfToItrf / ItrfToGcrf /
  /// MoonCiToPa / MoonPaToCi to their original analytic computation.
  void ClearFrameConversionFit();

  /// @brief True if a SPICE-fitted Earth orientation parameter (EOP) model
  /// covering `t_tdb` is currently active (see InitFrameConversionFromSpice
  /// and ComputeEopFromSpice).
  bool HasFittedEarthOrientation(Real t_tdb);

  /// @brief True if a SPICE-fitted lunar orientation (MOON_CI<->MOON_PA)
  /// model covering `t_tdb` is currently active.
  bool HasFittedLunarOrientation(Real t_tdb);

  // ---------------------------------------------------------------------------
  // Solar-system planet body orientation (generic IAU model). Declared here as
  // BodyId-keyed helpers; the Frame-keyed predicates and the ConvertFrame wiring
  // live in frame_converter.h (which owns the Frame enum).
  // ---------------------------------------------------------------------------

  /// @brief True if an IAU orientation model is available for `body`
  /// (Mercury..Neptune).
  bool HasIauOrientation(BodyId body);

  /// @brief Rotation matrix from a planet-centered inertial frame (ICRF-aligned)
  /// to that planet's body-fixed frame at `t_tdb`, from the IAU linear orientation
  /// model. If `R_dot` is non-null it receives the exact time derivative (used to
  /// rotate velocities), dominated by the prime-meridian spin rate.
  /// @param t_tdb  Epoch [s, TDB since J2000]
  /// @param body   Planet body id (Mercury..Neptune)
  /// @param R_dot  Optional out: d/dt of the rotation matrix [1/s]
  Mat3 RotBodyCiToFixed(Real t_tdb, BodyId body, Mat3* R_dot = nullptr);

  /// @brief Translate a planet-centered inertial (ICRF-aligned) state to the
  /// SSB-centered ICRF hub, by adding the planet's ICRF ephemeris state.
  Vec6 PlanetCiToIcrf(Real t_tdb, const Vec6& rv_ci, BodyId body);
  /// @brief Vec3 (position-only) overload of PlanetCiToIcrf().
  Vec3 PlanetCiToIcrf(Real t_tdb, const Vec3& r_ci, BodyId body);
  /// @brief Inverse of PlanetCiToIcrf(): SSB-centered ICRF -> planet-centered
  /// inertial (ICRF-aligned).
  Vec6 IcrfToPlanetCi(Real t_tdb, const Vec6& rv_icrf, BodyId body);
  /// @brief Vec3 (position-only) overload of IcrfToPlanetCi().
  Vec3 IcrfToPlanetCi(Real t_tdb, const Vec3& r_icrf, BodyId body);
  /// @brief Rotate a planet-centered inertial state into the body-fixed frame
  /// (applies RotBodyCiToFixed and, for a Vec6, its derivative).
  Vec6 BodyCiToFixed(Real t_tdb, const Vec6& rv_ci, BodyId body);
  /// @brief Vec3 (position-only) overload of BodyCiToFixed().
  Vec3 BodyCiToFixed(Real t_tdb, const Vec3& r_ci, BodyId body);
  /// @brief Inverse of BodyCiToFixed(): body-fixed -> planet-centered inertial.
  Vec6 BodyFixedToCi(Real t_tdb, const Vec6& rv_fixed, BodyId body);
  /// @brief Vec3 (position-only) overload of BodyFixedToCi().
  Vec3 BodyFixedToCi(Real t_tdb, const Vec3& r_fixed, BodyId body);

}  // namespace lupnt

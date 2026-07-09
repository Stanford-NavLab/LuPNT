/*
Greg Miller (gmiller@gregmiller.net) 2022
Released as public domain
http://www.celestialprogramming.com/

Class to read binary versions of JPL's Development Ephemeris.  Files in
the propper format can be obtained from:
ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux

#    Properties       Units          Center Description
0    x,y,z            km             SSB    Mercury
1    x,y,z            km             SSB    Venus
2    x,y,z            km             SSB    Earth-Moon barycenter
3    x,y,z            km             SSB    Mars
4    x,y,z            km             SSB    Jupiter
5    x,y,z            km             SSB    Saturn
6    x,y,z            km             SSB    Uranus
7    x,y,z            km             SSB    Neptune
8    x,y,z            km             SSB    Pluto
9    x,y,z            km             Earth  Moon (geocentric)
10   x,y,z            km             SSB    Sun
11   dPsi,dEps        radians               Earth Nutations in lon and obliquity
12   phi,theta,psi    radians               Lunar mantle libration
13   Ox,Oy,Oz         radians/day           Lunar mantle angular velocity
14   t                seconds               TT-TDB (at geocenter)

Example: (prints x coordinate of venus using first JD available)

JPLDE.DE de = new JPLDE.DE(@"E:\Astronomy\_Ephemeris\JPLDEBinaries\jpleph.405");
Console.WriteLine(de.getPlanet(1, de.getHeader().jdStart)[0]);

24857048.3412405
*/
#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "lupnt/core/constants.h"
#include "lupnt/interfaces/spice.h"

#define MAXCOEFF 1020

namespace lupnt {

  /// @brief Get the TT-TDB time difference at a geocentric epoch from the DE ephemeris.
  ///
  /// Intended to provide the periodic TT-TDB term (DE block 14, "t [s], TT-TDB at geocenter")
  /// for high-precision time conversions, as an alternative to the analytic series used
  /// elsewhere in conversions/time_conversions.cc.
  ///
  /// @param t_tdb Epoch [s, TDB since J2000]
  /// @return      TT - TDB at the geocenter [s]
  /// @note Currently unimplemented (throws std::runtime_error).
  double GetTtTdbDifference(double t_tdb);

  /// @brief Evaluate the DE Chebyshev "lunar mantle libration" angles (and optionally their time
  /// derivatives) at a given epoch.
  ///
  /// Called by RotMoonCiToPa/RotMoonPaToCi (frame_conversions.cc) to reconstruct the Moon
  /// principal-axes (PA) <-> Moon-centered-inertial (CI) rotation as
  /// R = Rz(psi)*Rx(theta)*Rz(phi), unless a SPICE-fitted lunar orientation model is active (see
  /// InitFrameConversionFromSpice/HasFittedLunarOrientation in frame_conversions.h). Lazily loads
  /// the DE440 ephemeris (LoadEphemerisData) on first use and evaluates the Chebyshev
  /// polynomials for ephemeris block 12 (lunar mantle libration).
  ///
  /// @param t_tdb       Epoch [s, TDB since J2000]
  /// @param compute_vel If true, also compute the time derivatives phi_dot/theta_dot/psi_dot
  ///                    (default: true); if false, the last 3 elements of the return value are
  ///                    zero.
  /// @return            6-vector [phi, theta, psi, phi_dot, theta_dot, psi_dot]: ZXZ Euler angles
  ///                    [rad] of the lunar mantle libration and their time derivatives [rad/s].
  Vec6 GetLunarMantleData(Real t_tdb, bool compute_vel = true);

  /// @brief Vectorized (per-epoch) overload of GetLunarMantleData(Real, bool); each row of the
  /// returned matrix is the 6-vector for the corresponding epoch in `t_tdb`.
  MatX6 GetLunarMantleData(VecX t_tdb, bool compute_vel = true);

  /// @brief Get the position and velocity of `target` relative to `center` from the DE
  /// ephemeris, expressed in `frame`.
  ///
  /// Used throughout the simulator wherever a third-body state is needed: third-body
  /// gravitational perturbations in NumericalOrbitDynamics::CalcContrib
  /// (dynamics/numerical_orbit_dynamics.cc), Sun/body direction vectors for attitude pointing
  /// (dynamics/attitude_dynamics.cc), light-time/Shapiro-delay corrections in
  /// conversions/time_conversions.cc and the measurement models, and joint orbit/clock dynamics
  /// (dynamics/joint_orbit_clock_dynamics.cc). Lazily loads the DE440 ephemeris on first call.
  /// Internally evaluates the relevant DE Chebyshev blocks (combining Earth, Earth-Moon
  /// barycenter, and Moon data as needed for Earth/Moon/EMB queries) and, unless `frame` is GCRF,
  /// rotates both states into `frame` via ConvertFrame before differencing.
  ///
  /// @param t_tdb   Epoch [s, TDB since J2000]
  /// @param center  Body at the origin of the returned relative state.
  /// @param target  Body whose state is returned.
  /// @param frame   Reference frame in which the position/velocity are expressed (and in which
  ///                the subtraction `target - center` is performed).
  /// @return        6-vector [rx, ry, rz, vx, vy, vz]: position [m] and velocity [m/s] of
  ///                `target` relative to `center`, expressed in `frame`. Zero if
  ///                `center == target`.
  Vec6 GetBodyPosVel(Real t_tdb, BodyId center, BodyId target, Frame frame);

  /// @brief Overload of GetBodyPosVel() that uses the natural center body of `frame` (via
  /// `frame_centers`) as the reference body.
  Vec6 GetBodyPosVel(Real t_tdb, BodyId target, Frame frame);

  /// @brief Overload of GetBodyPosVel(Real, BodyId, BodyId, Frame) that additionally rescales the
  /// result from SI units to `units` and applies a relativistic coordinate-scale correction (see
  /// CoordinateScale) before unit conversion.
  Vec6 GetBodyPosVel(Real t_tdb, BodyId center, BodyId target, Frame frame, const UnitSystem& units,
                     CoordinateScale scale = CoordinateScale::TDB);

  /// @brief Overload of GetBodyPosVel(Real, BodyId, Frame, const UnitSystem&, CoordinateScale)
  /// that uses the natural center body of `frame` as the reference body.
  Vec6 GetBodyPosVel(Real t_tdb, BodyId target, Frame frame, const UnitSystem& units,
                     CoordinateScale scale = CoordinateScale::TDB);

  /// @brief Vectorized (per-epoch) overload of GetBodyPosVel(Real, BodyId, Frame); each row of
  /// the returned matrix is the 6-vector for the corresponding epoch in `t_tdb`.
  MatX6 GetBodyPosVel(const VecX& t_tdb, BodyId target, Frame frame);

  /// @brief Vectorized (per-epoch) overload of GetBodyPosVel(Real, BodyId, BodyId, Frame); each
  /// row of the returned matrix is the 6-vector for the corresponding epoch in `t_tdb`.
  MatX6 GetBodyPosVel(const VecX& t_tdb, BodyId center, BodyId target, Frame frame);

  /// @brief Vectorized (per-epoch) overload of
  /// GetBodyPosVel(Real, BodyId, Frame, const UnitSystem&, CoordinateScale); each row of the
  /// returned matrix is the 6-vector for the corresponding epoch in `t_tdb`.
  MatX6 GetBodyPosVel(const VecX& t_tdb, BodyId target, Frame frame, const UnitSystem& units,
                      CoordinateScale scale = CoordinateScale::TDB);

  /// @brief Vectorized (per-epoch) overload of
  /// GetBodyPosVel(Real, BodyId, BodyId, Frame, const UnitSystem&, CoordinateScale); each row of
  /// the returned matrix is the 6-vector for the corresponding epoch in `t_tdb`.
  MatX6 GetBodyPosVel(const VecX& t_tdb, BodyId center, BodyId target, Frame frame,
                      const UnitSystem& units, CoordinateScale scale = CoordinateScale::TDB);

  /// @brief Get the position of `target` relative to `center` from the DE ephemeris, expressed
  /// in `frame`.
  ///
  /// Used wherever only a third-body position (not velocity) is needed -- e.g. Sun/body
  /// direction vectors for eclipse and attitude-pointing computations (dynamics/
  /// attitude_dynamics.cc), and third-body gravity gradients in NumericalOrbitDynamics
  /// (dynamics/numerical_orbit_dynamics.cc). Equivalent to
  /// `GetBodyPosVel(t_tdb, center, target, frame).head(3)`.
  ///
  /// @param t_tdb   Epoch [s, TDB since J2000]
  /// @param center  Body at the origin of the returned relative position.
  /// @param target  Body whose position is returned.
  /// @param frame   Reference frame in which the position is expressed.
  /// @return        Position of `target` relative to `center` [m], expressed in `frame`. Zero if
  ///                `center == target`.
  Vec3 GetBodyPos(Real t_tdb, BodyId center, BodyId target, Frame frame);

  /// @brief Overload of GetBodyPos() that uses the natural center body of `frame` (via
  /// `frame_centers`) as the reference body.
  Vec3 GetBodyPos(Real t_tdb, BodyId target, Frame frame);

  /// @brief Overload of GetBodyPos(Real, BodyId, BodyId, Frame) that additionally rescales the
  /// result from SI units to `units` and applies a relativistic coordinate-scale correction (see
  /// CoordinateScale) before unit conversion.
  Vec3 GetBodyPos(Real t_tdb, BodyId center, BodyId target, Frame frame, const UnitSystem& units,
                  CoordinateScale scale = CoordinateScale::TDB);

  /// @brief Overload of GetBodyPos(Real, BodyId, Frame, const UnitSystem&, CoordinateScale) that
  /// uses the natural center body of `frame` as the reference body.
  Vec3 GetBodyPos(Real t_tdb, BodyId target, Frame frame, const UnitSystem& units,
                  CoordinateScale scale = CoordinateScale::TDB);

  /// @brief Vectorized (per-epoch) overload of GetBodyPos(Real, BodyId, Frame); each row of the
  /// returned matrix is the position for the corresponding epoch in `t_tdb`.
  MatX3 GetBodyPos(const VecX& t_tdb, BodyId target, Frame frame);

  /// @brief Vectorized (per-epoch) overload of GetBodyPos(Real, BodyId, BodyId, Frame); each row
  /// of the returned matrix is the position for the corresponding epoch in `t_tdb`.
  MatX3 GetBodyPos(const VecX& t_tdb, BodyId center, BodyId target, Frame frame);

  /// @brief Vectorized (per-epoch) overload of
  /// GetBodyPos(Real, BodyId, Frame, const UnitSystem&, CoordinateScale); each row of the
  /// returned matrix is the position for the corresponding epoch in `t_tdb`.
  MatX3 GetBodyPos(const VecX& t_tdb, BodyId target, Frame frame, const UnitSystem& units,
                   CoordinateScale scale = CoordinateScale::TDB);

  /// @brief Vectorized (per-epoch) overload of
  /// GetBodyPos(Real, BodyId, BodyId, Frame, const UnitSystem&, CoordinateScale); each row of the
  /// returned matrix is the position for the corresponding epoch in `t_tdb`.
  MatX3 GetBodyPos(const VecX& t_tdb, BodyId center, BodyId target, Frame frame,
                   const UnitSystem& units, CoordinateScale scale = CoordinateScale::TDB);

}  // namespace lupnt

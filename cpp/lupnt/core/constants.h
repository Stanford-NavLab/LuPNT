/**
 * @file constants.h
 * @author Stanford NAV LAB
 * @brief List of constants
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "lupnt/core/definitions.h"
#include "lupnt/core/file.h"

namespace lupnt {

  // Math constants
  static constexpr double PI = 3.14159265358979323846264338327950288419716939937511;
  static constexpr double TWO_PI = 2.0 * PI;
  static constexpr double PI_OVER_TWO = PI / 2.0;
  static constexpr double E = 2.71828182845904523536028747135266249775724709369996;
  static constexpr double EPS = 1.0e-16;

  // Unit system constants ******************************************************
  //
  // LuPNT's historical constants remain SI-valued for backward compatibility.
  // Use UnitSystem and GetPhysicalConstants(...) when a simulation state is
  // expressed in another coherent unit set, e.g. kilometers, seconds, kilograms.
  static constexpr double METER = 1.0;      // [m/m]
  static constexpr double KILOMETER = 1e3;  // [m/km]
  static constexpr double SECOND = 1.0;     // [s/s]
  static constexpr double KILOGRAM = 1.0;   // [kg/kg]

  /// @brief Raise `base` to an integer power (positive, negative, or zero), without
  /// relying on `std::pow`.
  ///
  /// Used by `UnitSystem::FromSI`/`UnitSystem::ToSI` to build dimensional scale
  /// factors (length^a * time^b * mass^c) at compile time.
  ///
  /// @param base     Base value
  /// @param exponent Integer exponent (may be negative)
  /// @return         `base` raised to `exponent`
  constexpr double PowInt(double base, int exponent) {
    double result = 1.0;
    int n = exponent < 0 ? -exponent : exponent;
    for (int i = 0; i < n; i++) result *= base;
    return exponent < 0 ? 1.0 / result : result;
  }

  /**
   * @brief Coherent simulation unit system used to scale SI constants.
   *
   * Each field stores the number of SI base units in one simulation unit.
   * For example, `KM_S_KG_UNITS.length == 1000.0`, so a quantity expressed
   * in meters is divided by 1000 when converted to simulation length units.
   */
  struct UnitSystem {
    double length = METER;   // SI meters per simulation length unit
    double time = SECOND;    // SI seconds per simulation time unit
    double mass = KILOGRAM;  // SI kilograms per simulation mass unit

    /// @brief Convert an SI value into this unit system using dimensional powers.
    ///
    /// Called by the per-quantity helpers below (`Length`, `Velocity`,
    /// `GravitationalParameter`, etc.) and by `GetPhysicalConstants(units)` to rescale
    /// every SI-valued physical constant into a simulation's chosen unit system (e.g.
    /// `KM_S_KG_UNITS`) before propagation.
    ///
    /// @param value        Value in SI units
    /// @param length_power Exponent of length in the quantity's dimension
    /// @param time_power   Exponent of time in the quantity's dimension
    /// @param mass_power   Exponent of mass in the quantity's dimension
    /// @return             Value expressed in this unit system
    constexpr double FromSI(double value, int length_power, int time_power = 0,
                            int mass_power = 0) const {
      return value
             / (PowInt(length, length_power) * PowInt(time, time_power) * PowInt(mass, mass_power));
    }

    /// @brief Convert a value expressed in this unit system back into SI units.
    ///
    /// Inverse of `FromSI`; used wherever a quantity stored in a simulation's working
    /// unit system needs to be reported or compared in SI units.
    ///
    /// @param value        Value in this unit system
    /// @param length_power Exponent of length in the quantity's dimension
    /// @param time_power   Exponent of time in the quantity's dimension
    /// @param mass_power   Exponent of mass in the quantity's dimension
    /// @return             Value in SI units
    constexpr double ToSI(double value, int length_power, int time_power = 0,
                          int mass_power = 0) const {
      return value * PowInt(length, length_power) * PowInt(time, time_power)
             * PowInt(mass, mass_power);
    }

    /// @brief Convert a length [m] into this unit system.
    constexpr double Length(double value_m) const { return FromSI(value_m, 1); }
    /// @brief Convert an area [m^2] into this unit system.
    constexpr double Area(double value_m2) const { return FromSI(value_m2, 2); }
    /// @brief Convert a volume [m^3] into this unit system.
    constexpr double Volume(double value_m3) const { return FromSI(value_m3, 3); }
    /// @brief Convert a duration [s] into this unit system.
    constexpr double Duration(double value_s) const { return FromSI(value_s, 0, 1); }
    /// @brief Convert a mass [kg] into this unit system.
    constexpr double Mass(double value_kg) const { return FromSI(value_kg, 0, 0, 1); }
    /// @brief Convert a velocity [m/s] into this unit system.
    constexpr double Velocity(double value_m_s) const { return FromSI(value_m_s, 1, -1); }
    /// @brief Convert an acceleration [m/s^2] into this unit system.
    constexpr double Acceleration(double value_m_s2) const { return FromSI(value_m_s2, 1, -2); }
    /// @brief Convert a gravitational parameter GM [m^3/s^2] into this unit system.
    constexpr double GravitationalParameter(double value_m3_s2) const {
      return FromSI(value_m3_s2, 3, -2);
    }
    /// @brief Convert an angular velocity [rad/s] into this unit system.
    constexpr double AngularVelocity(double value_rad_s) const {
      return FromSI(value_rad_s, 0, -1);
    }
    /// @brief Convert a frequency [Hz] into this unit system.
    constexpr double Frequency(double value_hz) const { return FromSI(value_hz, 0, -1); }
    /// @brief Convert a force [N] into this unit system.
    constexpr double Force(double value_n) const { return FromSI(value_n, 1, -2, 1); }
    /// @brief Convert a pressure [Pa] into this unit system.
    constexpr double Pressure(double value_pa) const { return FromSI(value_pa, -1, -2, 1); }
    /// @brief Convert an area-to-mass ratio [m^2/kg] into this unit system (e.g. for SRP
    /// area-to-mass terms).
    constexpr double AreaPerMass(double value_m2_kg) const { return FromSI(value_m2_kg, 2, 0, -1); }

    /// @brief True if both unit systems use the same length, time, and mass scales.
    constexpr bool operator==(const UnitSystem& other) const {
      return length == other.length && time == other.time && mass == other.mass;
    }
    /// @brief True if the unit systems differ in any of length, time, or mass scale.
    constexpr bool operator!=(const UnitSystem& other) const { return !(*this == other); }
  };

  static constexpr UnitSystem SI_UNITS{METER, SECOND, KILOGRAM};
  static constexpr UnitSystem M_S_KG_UNITS = SI_UNITS;
  static constexpr UnitSystem KM_S_KG_UNITS{KILOMETER, SECOND, KILOGRAM};

  // Angle conversion
  static constexpr double RAD = PI / 180.0;               // [rad/deg]
  static constexpr double DEG = 180.0 / PI;               // [deg/rad]
  static constexpr double ARCSEC_DEG = 3600.0;            // [arcsec/deg]
  static constexpr double DEG_ARCSEC = 1.0 / ARCSEC_DEG;  // [deg/arcsec]
  static constexpr double RAD_ARCSEC = DEG_ARCSEC * RAD;  // [rad/arcsec]
  static constexpr double ARCSEC_RAD = 1.0 / RAD_ARCSEC;  // [arcsec/rad]

  // Length
  static constexpr double INCH_M = 0.0254;    // [in/m]
  static constexpr double FOOT_M = 0.3048;    // [ft/m]
  static constexpr double MILE_M = 1609.344;  // [mile/m]
  static constexpr double KM_M = 1e-3;        // [km/m]
  static constexpr double M_KM = 1e3;         // [m/km]
  static constexpr double MM_KM = 1.0e6;      // [mm/km]
  static constexpr double KM_MM = 1.0e-6;     // [km/mm]
  static constexpr double MM_M = 1e3;         // [mm/m]
  static constexpr double M_MM = 1e-3;        // [m/mm]
  static constexpr double M_CM = 1e-2;        // [m/cm]
  static constexpr double CM_M = 1e2;         // [cm/m]

  // Time system constants ******************************************************
  static constexpr double SECS_DAY = 86400.0;         // [s/day]
  static constexpr double SECS_HOUR = 3600.0;         // [s/hour]
  static constexpr double SECS_MINUTE = 60.0;         // [s/minute]
  static constexpr double MINS_HOUR = 60.0;           // [min/hour]
  static constexpr double MINS_DAY = 1440.0;          // [min/day]
  static constexpr double HOURS_DAY = 24.0;           // [hour/day]
  static constexpr double DAYS_SEC = 1.0 / SECS_DAY;  // [day/s]
  static constexpr double DAYS_WEEK = 7.0;            // [days/week]
  static constexpr double DAYS_YEAR = 365.25;         // [days/year]
  static constexpr double DAYS_CENTURY = 36525.00;    // [days/century]

  static constexpr double JD_MJD_OFFSET = 2400000.5;  // [days]
  static constexpr double TT_TAI_OFFSET = 32.184;     // [s]
  static constexpr double A1_TAI_OFFSET = 0.0343817;  // [s]

  static constexpr double JD_CCSDS_TAI = 2436203.5;  // [days]
  static constexpr double JD_J2000_TT = 2451545.0;   // [days]
  static constexpr double MJD_CCSDS_TAI = 36203.0;   // [days]
  static constexpr double MJD_J2000_TT = 51544.5;    // [days]

  // Vallado page 194
  static constexpr double MJD_COORDINATE_TAI = 2443144.5 - JD_MJD_OFFSET;  // [days]
  static constexpr double MJD_COORDINATE_TT_TCG_TCB
      = MJD_COORDINATE_TAI + TT_TAI_OFFSET / SECS_DAY;  // [days]

  // IAU/Turyshev (2026, ApJ 997:97) Table 1 relativistic rate constants.
  // L_B, L_G, L_L define: TDB = TCB − L_B·(TCB−T_0) + TDB_0
  //                        TT  = TCG − L_G·(TCG−T_0)
  //                        TL  = TCL − L_L·(TCL−T_0)
  static constexpr double L_B = 1.550519768e-8;   // TCB→TDB rate  (IAU 2006 B3)
  static constexpr double L_G = 6.969290134e-10;  // TCG→TT  rate  (IAU 1997 B1.9, defining)
  // Turyshev et al. 2025 (ApJ 985:140) Table 2 / Eq.(35), R_MQ = 1738.0 km.
  // Keep all quoted digits: truncating to 3.13905e-11 costs 4e-17 s/s (1.3 ns/yr) in TL.
  static constexpr double L_L = 3.139054e-11;    // TCL→TL  rate  (Turyshev 2025, selenoid)
  static constexpr double L_H = 1.48253624e-8;   // TCL−TCB mean rate
  static constexpr double L_M = 1.485675290e-8;  // TL−TCB  mean rate
  // L_EM = L_H - L_C (Turyshev et al. 2025, Table 2 and Eq.(76)).
  // Currently unused: TdbToLtMinusTt() integrates this quantity from the
  // ephemeris rather than reading the constant. Kept correct for a future
  // closed-form fast path.
  static constexpr double L_EM = 1.709385e-11;  // TL−TT   mean rate [s/s TDB]
  // DE405 TDB offset: TDB_0 = −65.5 μs  (Turyshev 2026, Table 1)
  static constexpr double TDB_0 = -65.5e-6;  // [s]

  /**
   * @brief Relativistic coordinate scale attached to dimensional quantities.
   *
   * This is separate from the epoch time scale (`Time`) and from numerical
   * units (`UnitSystem`). For IAU-style scaled systems, spatial coordinates and
   * gravitational parameters are rescaled by a defining constant, while
   * velocities remain invariant because length and coordinate time are scaled by
   * the same factor.
   */
  enum class CoordinateScale {
    TCB,  ///< Barycentric coordinate quantities compatible with TCB.
    TDB,  ///< Barycentric coordinate quantities compatible with TDB.
    TCG,  ///< Geocentric coordinate quantities compatible with TCG.
    TT,   ///< Geocentric coordinate quantities compatible with TT.
    TCL,  ///< Lunicentric coordinate quantities compatible with TCL.
    TL,   ///< Lunar-time-compatible lunicentric quantities.
  };

  /// @brief Check whether two `CoordinateScale`s are related by a single constant IAU
  /// scale factor (i.e. belong to the same barycentric/geocentric/lunicentric pair, or
  /// are identical).
  ///
  /// Used as a guard by `CoordinateScaleRatio` (returns NaN if not convertible) and by
  /// `CheckedCoordinateScaleRatio` (throws if not convertible) before any of the
  /// `Scale*ForCoordinateScale*` helpers below rescale a position, state, or
  /// gravitational parameter between two `CoordinateScale`s.
  ///
  /// @param from Source coordinate scale
  /// @param to   Target coordinate scale
  /// @return     True if `from` and `to` are related by a constant scale factor
  constexpr bool AreCoordinateScalesConvertible(CoordinateScale from, CoordinateScale to) {
    if (from == to) return true;
    bool barycentric_pair = (from == CoordinateScale::TCB || from == CoordinateScale::TDB)
                            && (to == CoordinateScale::TCB || to == CoordinateScale::TDB);
    bool geocentric_pair = (from == CoordinateScale::TCG || from == CoordinateScale::TT)
                           && (to == CoordinateScale::TCG || to == CoordinateScale::TT);
    bool lunicentric_pair = (from == CoordinateScale::TCL || from == CoordinateScale::TL)
                            && (to == CoordinateScale::TCL || to == CoordinateScale::TL);
    return barycentric_pair || geocentric_pair || lunicentric_pair;
  }

  /// @brief Look up the IAU defining scale factor (1 - L_x) for a `CoordinateScale`.
  ///
  /// `TCB`/`TCG`/`TCL` (proper "coordinate time" scales) map to 1.0, while
  /// `TDB`/`TT`/`TL` map to `1 - L_B`/`1 - L_G`/`1 - L_L` respectively. Used by
  /// `CoordinateScaleRatio`/`CheckedCoordinateScaleRatio` to form the ratio between two
  /// scales.
  ///
  /// @param scale Coordinate scale
  /// @return      IAU defining scale factor for `scale`
  constexpr double CoordinateScaleFactor(CoordinateScale scale) {
    switch (scale) {
      case CoordinateScale::TCB:
      case CoordinateScale::TCG:
      case CoordinateScale::TCL: return 1.0;
      case CoordinateScale::TDB: return 1.0 - L_B;
      case CoordinateScale::TT: return 1.0 - L_G;
      case CoordinateScale::TL: return 1.0 - L_L;
    }
    return 1.0;
  }

  /// @brief Ratio `CoordinateScaleFactor(to) / CoordinateScaleFactor(from)` used to
  /// rescale dimensional quantities between two `CoordinateScale`s, or NaN if the two
  /// scales are not related by a constant factor (see `AreCoordinateScalesConvertible`).
  ///
  /// @param from Source coordinate scale
  /// @param to   Target coordinate scale
  /// @return     Scale ratio `to`/`from`, or NaN if not convertible
  constexpr double CoordinateScaleRatio(CoordinateScale from, CoordinateScale to) {
    return AreCoordinateScalesConvertible(from, to)
               ? CoordinateScaleFactor(to) / CoordinateScaleFactor(from)
               : std::numeric_limits<double>::quiet_NaN();
  }

  /// @brief Rescale a length-dimensioned value from one `CoordinateScale` to another,
  /// or return NaN if the scales are not related by a constant factor. Non-throwing
  /// counterpart of `ScaleLengthForCoordinateScaleChecked`.
  constexpr double ScaleLengthForCoordinateScale(double value, CoordinateScale from,
                                                 CoordinateScale to) {
    return value * CoordinateScaleRatio(from, to);
  }

  /// @brief Rescale a gravitational parameter GM from one `CoordinateScale` to another,
  /// or return NaN if the scales are not related by a constant factor. Non-throwing
  /// counterpart of `ScaleGravitationalParameterForCoordinateScaleChecked`.
  constexpr double ScaleGravitationalParameterForCoordinateScale(double value, CoordinateScale from,
                                                                 CoordinateScale to) {
    return value * CoordinateScaleRatio(from, to);
  }

  /// @brief Ratio `CoordinateScaleFactor(to) / CoordinateScaleFactor(from)`, throwing if
  /// the two scales are not related by a constant IAU scale factor.
  ///
  /// Used by `ScaleLengthForCoordinateScaleChecked`,
  /// `ScaleGravitationalParameterForCoordinateScaleChecked`,
  /// `ScalePositionForCoordinateScale`, and `ScaleStateForCoordinateScale` (and exposed
  /// to Python as `coordinate_scale_ratio`) to convert ephemeris positions/states
  /// produced in `CoordinateScale::TDB` (e.g. by `GetBodyPosVel`/`GetBodyPos` in
  /// `lupnt/interfaces/kernels.cc`) into another requested coordinate scale.
  ///
  /// @param from Source coordinate scale
  /// @param to   Target coordinate scale
  /// @return     Scale ratio `to`/`from`
  /// @throws std::invalid_argument if `from` and `to` are not related by a constant scale factor
  inline double CheckedCoordinateScaleRatio(CoordinateScale from, CoordinateScale to) {
    if (!AreCoordinateScalesConvertible(from, to)) {
      throw std::invalid_argument(
          "Coordinate scales are not related by a constant IAU scale factor");
    }
    return CoordinateScaleFactor(to) / CoordinateScaleFactor(from);
  }

  /// @brief Rescale a length-dimensioned value from one `CoordinateScale` to another.
  /// @throws std::invalid_argument if the scales are not related by a constant factor (see
  ///         `CheckedCoordinateScaleRatio`)
  inline double ScaleLengthForCoordinateScaleChecked(double value, CoordinateScale from,
                                                     CoordinateScale to) {
    return value * CheckedCoordinateScaleRatio(from, to);
  }

  /// @brief Rescale a gravitational parameter GM from one `CoordinateScale` to another.
  /// @throws std::invalid_argument if the scales are not related by a constant factor (see
  ///         `CheckedCoordinateScaleRatio`)
  inline double ScaleGravitationalParameterForCoordinateScaleChecked(double value,
                                                                     CoordinateScale from,
                                                                     CoordinateScale to) {
    return value * CheckedCoordinateScaleRatio(from, to);
  }

  /// @brief Rescale a Cartesian position vector from one `CoordinateScale` to another.
  ///
  /// Called by `GetBodyPos(..., units, scale)` in `lupnt/interfaces/kernels.cc` to convert a
  /// TDB-scale ephemeris position into the requested `CoordinateScale` before unit
  /// conversion.
  ///
  /// @param r    Position vector in `from` scale [m]
  /// @param from Source coordinate scale
  /// @param to   Target coordinate scale
  /// @return     Position vector rescaled to `to` [m]
  /// @throws std::invalid_argument if the scales are not related by a constant factor
  inline Vec3 ScalePositionForCoordinateScale(const Vec3& r, CoordinateScale from,
                                              CoordinateScale to) {
    return r * CheckedCoordinateScaleRatio(from, to);
  }

  /// @brief Rescale a Cartesian position+velocity state from one `CoordinateScale` to
  /// another, rescaling the position (via `ScalePositionForCoordinateScale`) and
  /// leaving the velocity unchanged (velocities are invariant under IAU coordinate
  /// scale changes because length and coordinate time scale by the same factor).
  ///
  /// Called by `GetBodyPosVel(..., units, scale)` in `lupnt/interfaces/kernels.cc` to convert
  /// a TDB-scale ephemeris state into the requested `CoordinateScale` before unit
  /// conversion.
  ///
  /// @param rv   Position+velocity state in `from` scale [m, m/s]
  /// @param from Source coordinate scale
  /// @param to   Target coordinate scale
  /// @return     State with position rescaled to `to` [m, m/s]
  /// @throws std::invalid_argument if the scales are not related by a constant factor
  inline Vec6 ScaleStateForCoordinateScale(const Vec6& rv, CoordinateScale from,
                                           CoordinateScale to) {
    Vec6 out = rv;
    out.head(3) = ScalePositionForCoordinateScale(rv.head(3), from, to);
    return out;
  }

  // Coordinate system constants DE440 *******************************************
  static constexpr double GM_SUN = 132712440041.279419e9;          // [m^3/s^2]
  static constexpr double GM_MERCURY = 22031.868551e9;             // [m^3/s^2]
  static constexpr double GM_VENUS = 324858.592000e9;              // [m^3/s^2]
  static constexpr double GM_EARTH = 398600.435507e9;              // [m^3/s^2]
  static constexpr double GM_MOON = 4902.800118e9;                 // [m^3/s^2]
  static constexpr double GM_MARS_SYSTEM = 42828.375816e9;         // [m^3/s^2]
  static constexpr double GM_JUPITER_SYSTEM = 126712764.100000e9;  // [m^3/s^2]
  static constexpr double GM_SATURN_SYSTEM = 37940584.841800e9;    // [m^3/s^2]
  static constexpr double GM_URANUS_SYSTEM = 5794556.400000e9;     // [m^3/s^2]
  static constexpr double GM_NEPTUNE_SYSTEM = 6836527.100580e9;    // [m^3/s^2]
  static constexpr double GM_PLUTO_SYSTEM = 977.000000e9;          // [m^3/s^2]
  // Planet-only mass parameters (JPL SSD planetary physical parameters).
  //
  // A planet-only GM is ~0.9996x its system GM -- the moons are only ~4e-4 of
  // the system mass -- so each value below sits just under the corresponding
  // *_SYSTEM constant above. A value one tenth of its *_SYSTEM counterpart is
  // wrong by a decimal place.
  //
  // Which to use: DE440 supplies only *system barycenters* for the outer
  // planets, so a potential or third-body term evaluated at the position
  // returned for BodyId::JUPITER etc. should be paired with the *_SYSTEM
  // value. Use these planet-only values only when the planet centre itself is
  // meant.
  static constexpr double GM_MARS = 0.4282837566395650e14;  // [m^3/s^2]
  static constexpr double GM_JUPITER = 126686531.900e9;     // [m^3/s^2]
  static constexpr double GM_SATURN = 37931206.234e9;       // [m^3/s^2]
  static constexpr double GM_URANUS = 5793951.256e9;        // [m^3/s^2]
  static constexpr double GM_NEPTUNE = 6835099.970e9;       // [m^3/s^2]

  static constexpr double GM_CERES = 62.62890e9;   // [m^3/s^2]
  static constexpr double GM_VESTA = 17.288245e9;  // [m^3/s^2]

  static constexpr double D_EARTH_MOON = 384400.0e3;  // [m]
  static constexpr double D_EARTH_EMB = 4671.0e3;     // [m]
  static constexpr double R_EARTH = 6378.137e3;       // [m]
  static constexpr double R_MOON = 1737.4e3;          // [m]
  static constexpr double R_SUN = 696342.0e3;         // [m]
  static constexpr double R_MERCURY = 2439.7e3;       // [m]
  static constexpr double R_VENUS = 6051.8e3;         // [m]
  static constexpr double R_MARS = 3396.0e3;          // [m]
  static constexpr double R_JUPITER = 71492.0e3;      // [m]
  static constexpr double R_SATURN = 60268.0e3;       // [m]
  static constexpr double R_URANUS = 25559.0e3;       // [m]
  static constexpr double R_NEPTUNE = 24764.0e3;      // [m]
  static constexpr double R_PLUTO = 1188.3e3;         // [m]

  static constexpr double OMEGA_EARTH_MOON = 2.6617e-6;             // [rad/s]
  static constexpr double D_MOON_EMB = D_EARTH_MOON - D_EARTH_EMB;  // [m]

  // Flattening
  static constexpr double WGS84_A = 6378.137e3;           // [m]
  static constexpr double WGS84_F = 1.0 / 298.257223563;  // [-]

  static constexpr double SUN_F = 0.0000;        // [-]
  static constexpr double EARTH_F = WGS84_F;     // [-]
  static constexpr double MOON_F = 0.0012;       // [-]
  static constexpr double MERCURY_F = 0.0009;    // [-]
  static constexpr double VENUS_F = 0.0000;      // [-]
  static constexpr double MARS_F = 1.0 / 169.8;  // [-]
  static constexpr double JUPITER_F = 0.06487;   // [-]
  static constexpr double SATURN_F = 0.09796;    // [-]
  static constexpr double URANUS_F = 0.02293;    // [-]
  static constexpr double NEPUTUNE_F = 0.01708;  // [-]

  // Sideral Rotation rate
  static constexpr double OMEGA_SUN = 2.903e-6;           // [rad/s]
  static constexpr double OMEGA_MERCURY = 1.244e-5;       // [rad/s]
  static constexpr double OMEGA_VENUS = -1.9521515e-7;    // [rad/s]
  static constexpr double OMEGA_EARTH = 7.2921151467e-5;  // [rad/s]
  static constexpr double OMEGA_MOON = 2.6617e-6;         // [rad/s]
  static constexpr double OMEGA_MARS = 7.0882185e-5;      // [rad/s]
  static constexpr double OMEGA_JUPITER = 1.758e-4;       // [rad/s]
  static constexpr double OMEGA_SATURN = 1.624e-4;        // [rad/s]
  static constexpr double OMEGA_URANUS = -1.036e-4;       // [rad/s]
  static constexpr double OMEGA_NEPTUNE = 1.083e-4;       // [rad/s]

  // Spherical harmonics
  static constexpr double J2_EARTH = 1.08262668e-3;
  // static constexpr double J2_MOON = 9.08901807506000e-5;
  static constexpr double J2_MOON
      = 9.094278450270e-5;  // Zonal value adjusted for permanent tide - Rigid J2
  static constexpr double C22_MOON
      = 3.470983013194e-5;  // Sectorial value adjusted for perm. tide - Rigid C22
  static constexpr double J2_MARS = 1.96045e-3;  // J2 value for Mars
  // Unnormalized second-degree zonal harmonic of the Sun, as estimated in
  // DE440 (Park et al. 2021, AJ 161:105). Used by the solar-oblateness
  // potential term w_LE in the DE440 Eq. (3) TDB-TT relation.
  static constexpr double J2_SUN = 2.246e-7;
  // Mean obliquity of the ecliptic at J2000 (Park et al. 2021, Eq. 24 at T=0):
  // 84381".448. Used to obtain the heliocentric ecliptic latitude of Earth.
  static constexpr double OBLIQUITY_J2000 = 84381.448 / 3600.0 * PI / 180.0;  // [rad]

  // Transformations Between GCRF and Mean Equator and Equinox at J2000
  static constexpr double FRAME_BIAS_XI0 = -16.6170e-3 * RAD_ARCSEC;   // [rad]
  static constexpr double FRAME_BIAS_ETA0 = -6.8192e-3 * RAD_ARCSEC;   // [rad]
  static constexpr double FRAME_BIAS_DALPHA0 = -14.6e-3 * RAD_ARCSEC;  // [rad]

  // Solar Radiation Pressure Constants
  /// Astronomical unit [m]. Exact by the IAU (2012) definition.
  static constexpr double AU = 149597870700.0;

  // --- Small-body populations that DE440 integrates but for which LuPNT carries
  // --- no ephemerides. Modelled as uniform circular rings (see RingPotential()).
  //
  // Values are taken from DE440's own integration header (header.440t), not from
  // independently published masses: DE440 fitted these during the integration, so
  // they are what reproduces its TT-TDB. LuPNT has no ephemerides for the bodies,
  // so each population is modelled as one uniform circular ring.
  //
  // Main asteroid belt: the 343 discrete asteroids DE440 integrates, summed from
  // its own header (header.440t, MA0001..MA1467 plus 14 extras). MA0001 is Ceres
  // at 6.2629e10 m^3/s^2. Total 1.7053e11 m^3/s^2 = 12.85e-10 Msun -- 1.05x the
  // independently published belt mass (Pitjeva EPM2014, 12.25e-10 Msun).
  static constexpr double GM_ASTEROID_BELT = 1.7053e11;  // [m^3/s^2]
  static constexpr double A_ASTEROID_BELT = 2.7 * AU;    // [m] mean belt radius
  // Kuiper belt: DE440's OWN constants, from the DE440/LE440 integration header
  // (header.440t, GROUP 1041), converted from AU^3/day^2 (factor 4.484859e23):
  //
  //   MA8201..MA8236  36 equal point masses forming the circular ring at 44 au,
  //                   0.552276997169882142e-12 each  ->  8.9168e12 m^3/s^2
  //   MA8001..MA8030  the 30 individual KBOs         ->  2.3154e12 m^3/s^2
  //                                           total  ->  1.1232e13 m^3/s^2
  //
  // These are the values DE440 actually integrated, so they are what reproduces
  // its TT-TDB. An independent published Kuiper-belt mass does not: Pitjeva &
  // Pitjev 2018 give 1.97e-2 M_Earth = 7.85e12, 1.43x smaller, because DE440
  // *fitted* its ring mass rather than adopting a published total.
  //
  // The 30 discrete KBOs are lumped in at the ring radius since LuPNT has no
  // ephemerides for them; several orbit beyond 44 au, so this slightly
  // overestimates their potential.
  static constexpr double GM_KUIPER_BELT = 1.1232e13;  // [m^3/s^2]
  static constexpr double A_KUIPER_BELT = 44.0 * AU;   // [m] DE440 ring radius
  /// Default mean total solar irradiance at 1 AU [W/m^2].
  ///
  /// 1360.8 +/- 0.5 is the measured TSI (Kopp & Lean 2011, SORCE/TIM); the value varies by
  /// ~0.1% over the solar cycle. The widely quoted 1367 predates TIM-era radiometry and is
  /// ~0.5% high. Missions usually specify their own figure -- override it per force model
  /// with `NBodyDynamics::SetSolarFlux()` rather than relying on this default.
  static constexpr double SOLAR_FLUX_AU = 1361;
  static constexpr double C = 299792458;  // Light speed [m/s]
  /// Solar radiation pressure at 1 AU [N/m^2] for `SOLAR_FLUX_AU`.
  static constexpr double P_SUN = SOLAR_FLUX_AU / C;

  /**
   * @brief Common physical constants scaled into one coherent unit system.
   *
   * The legacy global constants, such as `GM_EARTH` and `R_MOON`, remain in SI
   * units. Use `GetPhysicalConstants(KM_S_KG_UNITS)` when propagating states in
   * kilometers, seconds, and kilograms.
   */
  struct PhysicalConstants {
    double GM_SUN;
    double GM_MERCURY;
    double GM_VENUS;
    double GM_EARTH;
    double GM_MOON;
    double GM_MARS_SYSTEM;
    double GM_JUPITER_SYSTEM;
    double GM_SATURN_SYSTEM;
    double GM_URANUS_SYSTEM;
    double GM_NEPTUNE_SYSTEM;
    double GM_PLUTO_SYSTEM;
    double GM_MARS;
    double GM_JUPITER;
    double GM_SATURN;
    double GM_URANUS;
    double GM_NEPTUNE;
    double GM_CERES;
    double GM_VESTA;

    double D_EARTH_MOON;
    double D_EARTH_EMB;
    double D_MOON_EMB;
    double R_EARTH;
    double R_MOON;
    double R_SUN;
    double R_MERCURY;
    double R_VENUS;
    double R_MARS;
    double R_JUPITER;
    double R_SATURN;
    double R_URANUS;
    double R_NEPTUNE;
    double R_PLUTO;
    double WGS84_A;

    double OMEGA_EARTH_MOON;
    double OMEGA_SUN;
    double OMEGA_MERCURY;
    double OMEGA_VENUS;
    double OMEGA_EARTH;
    double OMEGA_MOON;
    double OMEGA_MARS;
    double OMEGA_JUPITER;
    double OMEGA_SATURN;
    double OMEGA_URANUS;
    double OMEGA_NEPTUNE;

    double AU;
    double C;
    double P_SUN;

    CoordinateScale coordinate_scale = CoordinateScale::TDB;
  };

  /// @brief Build a `PhysicalConstants` bundle (GM/radii/rotation rates/AU/c/...) with
  /// every dimensional quantity rescaled from SI into `units`, at the default
  /// `CoordinateScale::TDB`.
  ///
  /// Called by orbit/clock dynamics models (e.g.
  /// `NumericalOrbitDynamics::CalcContrib`, `JointOrbitClockDynamics`,
  /// `Clock`) at setup time to get gravitational parameters and body radii
  /// consistent with the dynamics' working unit system (e.g. `KM_S_KG_UNITS`), rather
  /// than repeatedly converting the SI-valued `GM_EARTH`/`R_MOON`/etc. constants.
  ///
  /// @param units Target coherent unit system (default: `SI_UNITS`)
  /// @return      Physical constants with each dimensional quantity scaled from SI to
  ///              `units`, `coordinate_scale == CoordinateScale::TDB`
  constexpr PhysicalConstants GetPhysicalConstants(const UnitSystem& units = SI_UNITS) {
    return {
        .GM_SUN = units.GravitationalParameter(GM_SUN),
        .GM_MERCURY = units.GravitationalParameter(GM_MERCURY),
        .GM_VENUS = units.GravitationalParameter(GM_VENUS),
        .GM_EARTH = units.GravitationalParameter(GM_EARTH),
        .GM_MOON = units.GravitationalParameter(GM_MOON),
        .GM_MARS_SYSTEM = units.GravitationalParameter(GM_MARS_SYSTEM),
        .GM_JUPITER_SYSTEM = units.GravitationalParameter(GM_JUPITER_SYSTEM),
        .GM_SATURN_SYSTEM = units.GravitationalParameter(GM_SATURN_SYSTEM),
        .GM_URANUS_SYSTEM = units.GravitationalParameter(GM_URANUS_SYSTEM),
        .GM_NEPTUNE_SYSTEM = units.GravitationalParameter(GM_NEPTUNE_SYSTEM),
        .GM_PLUTO_SYSTEM = units.GravitationalParameter(GM_PLUTO_SYSTEM),
        .GM_MARS = units.GravitationalParameter(GM_MARS),
        .GM_JUPITER = units.GravitationalParameter(GM_JUPITER),
        .GM_SATURN = units.GravitationalParameter(GM_SATURN),
        .GM_URANUS = units.GravitationalParameter(GM_URANUS),
        .GM_NEPTUNE = units.GravitationalParameter(GM_NEPTUNE),
        .GM_CERES = units.GravitationalParameter(GM_CERES),
        .GM_VESTA = units.GravitationalParameter(GM_VESTA),
        .D_EARTH_MOON = units.Length(D_EARTH_MOON),
        .D_EARTH_EMB = units.Length(D_EARTH_EMB),
        .D_MOON_EMB = units.Length(D_MOON_EMB),
        .R_EARTH = units.Length(R_EARTH),
        .R_MOON = units.Length(R_MOON),
        .R_SUN = units.Length(R_SUN),
        .R_MERCURY = units.Length(R_MERCURY),
        .R_VENUS = units.Length(R_VENUS),
        .R_MARS = units.Length(R_MARS),
        .R_JUPITER = units.Length(R_JUPITER),
        .R_SATURN = units.Length(R_SATURN),
        .R_URANUS = units.Length(R_URANUS),
        .R_NEPTUNE = units.Length(R_NEPTUNE),
        .R_PLUTO = units.Length(R_PLUTO),
        .WGS84_A = units.Length(WGS84_A),
        .OMEGA_EARTH_MOON = units.AngularVelocity(OMEGA_EARTH_MOON),
        .OMEGA_SUN = units.AngularVelocity(OMEGA_SUN),
        .OMEGA_MERCURY = units.AngularVelocity(OMEGA_MERCURY),
        .OMEGA_VENUS = units.AngularVelocity(OMEGA_VENUS),
        .OMEGA_EARTH = units.AngularVelocity(OMEGA_EARTH),
        .OMEGA_MOON = units.AngularVelocity(OMEGA_MOON),
        .OMEGA_MARS = units.AngularVelocity(OMEGA_MARS),
        .OMEGA_JUPITER = units.AngularVelocity(OMEGA_JUPITER),
        .OMEGA_SATURN = units.AngularVelocity(OMEGA_SATURN),
        .OMEGA_URANUS = units.AngularVelocity(OMEGA_URANUS),
        .OMEGA_NEPTUNE = units.AngularVelocity(OMEGA_NEPTUNE),
        .AU = units.Length(AU),
        .C = units.Velocity(C),
        .P_SUN = units.Pressure(P_SUN),
        .coordinate_scale = CoordinateScale::TDB,
    };
  }

  /// @brief Rescale a `PhysicalConstants` bundle from one relativistic
  /// `CoordinateScale` to another.
  ///
  /// The stored DE-derived constants are treated as TDB-compatible. Constant
  /// coordinate scaling is therefore supported for the barycentric TDB/TCB pair.
  /// Requests for geocentric or lunicentric scales throw because those require a
  /// reference-system transformation, not a single scalar applied to this bundle.
  ///
  /// Length-like quantities and gravitational parameters scale by the IAU scale ratio
  /// (`CheckedCoordinateScaleRatio`). Angular rates scale by the inverse ratio because
  /// coordinate time and coordinate length share the same scale factor. Light speed is
  /// invariant. Called by `GetPhysicalConstants(units, scale)` to produce constants in
  /// a non-default coordinate scale.
  ///
  /// @param constants Physical constants currently expressed in `from`'s scale (with
  ///                   `constants.coordinate_scale == from`)
  /// @param from      Source coordinate scale
  /// @param to        Target coordinate scale
  /// @return          Rescaled constants with `coordinate_scale == to`
  /// @throws std::invalid_argument if `from` and `to` are not related by a constant scale factor
  inline PhysicalConstants ScalePhysicalConstantsForCoordinateScale(PhysicalConstants constants,
                                                                    CoordinateScale from,
                                                                    CoordinateScale to) {
    double ratio = CheckedCoordinateScaleRatio(from, to);
    double inverse_ratio = 1.0 / ratio;

    constants.GM_SUN *= ratio;
    constants.GM_MERCURY *= ratio;
    constants.GM_VENUS *= ratio;
    constants.GM_EARTH *= ratio;
    constants.GM_MOON *= ratio;
    constants.GM_MARS_SYSTEM *= ratio;
    constants.GM_JUPITER_SYSTEM *= ratio;
    constants.GM_SATURN_SYSTEM *= ratio;
    constants.GM_URANUS_SYSTEM *= ratio;
    constants.GM_NEPTUNE_SYSTEM *= ratio;
    constants.GM_PLUTO_SYSTEM *= ratio;
    constants.GM_MARS *= ratio;
    constants.GM_JUPITER *= ratio;
    constants.GM_SATURN *= ratio;
    constants.GM_URANUS *= ratio;
    constants.GM_NEPTUNE *= ratio;
    constants.GM_CERES *= ratio;
    constants.GM_VESTA *= ratio;

    constants.D_EARTH_MOON *= ratio;
    constants.D_EARTH_EMB *= ratio;
    constants.D_MOON_EMB *= ratio;
    constants.R_EARTH *= ratio;
    constants.R_MOON *= ratio;
    constants.R_SUN *= ratio;
    constants.R_MERCURY *= ratio;
    constants.R_VENUS *= ratio;
    constants.R_MARS *= ratio;
    constants.R_JUPITER *= ratio;
    constants.R_SATURN *= ratio;
    constants.R_URANUS *= ratio;
    constants.R_NEPTUNE *= ratio;
    constants.R_PLUTO *= ratio;
    constants.WGS84_A *= ratio;
    constants.AU *= ratio;

    constants.OMEGA_EARTH_MOON *= inverse_ratio;
    constants.OMEGA_SUN *= inverse_ratio;
    constants.OMEGA_MERCURY *= inverse_ratio;
    constants.OMEGA_VENUS *= inverse_ratio;
    constants.OMEGA_EARTH *= inverse_ratio;
    constants.OMEGA_MOON *= inverse_ratio;
    constants.OMEGA_MARS *= inverse_ratio;
    constants.OMEGA_JUPITER *= inverse_ratio;
    constants.OMEGA_SATURN *= inverse_ratio;
    constants.OMEGA_URANUS *= inverse_ratio;
    constants.OMEGA_NEPTUNE *= inverse_ratio;

    constants.coordinate_scale = to;
    return constants;
  }

  /// @brief Overload of `GetPhysicalConstants` that additionally rescales the result
  /// from the default `CoordinateScale::TDB` to `scale` via
  /// `ScalePhysicalConstantsForCoordinateScale`.
  /// @throws std::invalid_argument if `scale` is not related to TDB by a constant scale factor
  inline PhysicalConstants GetPhysicalConstants(const UnitSystem& units, CoordinateScale scale) {
    return ScalePhysicalConstantsForCoordinateScale(GetPhysicalConstants(units),
                                                    CoordinateScale::TDB, scale);
  }

  /// @brief Overload of `GetPhysicalConstants` using `SI_UNITS` and rescaling from
  /// `CoordinateScale::TDB` to `scale`.
  /// @throws std::invalid_argument if `scale` is not related to TDB by a constant scale factor
  inline PhysicalConstants GetPhysicalConstants(CoordinateScale scale) {
    return GetPhysicalConstants(SI_UNITS, scale);
  }

  // Ionosphere
  static const double IONOSPHERIC_CONSTANT = 40.3e16;  // Adjusted constant for units
  static const double L1_FREQ = 1.57542e9;             // L1 frequency [Hz]

  // File Path *******************************************************************
  static constexpr std::string_view TAI_UTC_FILENAME = "tai-utc.dat";
  static constexpr std::string_view EOP_FILENAME = "eopc04_08.62-now";
  // IERS finals ("Bulletin A") cache, written by LoadLatestEopFinalsFromIers and the default
  // path for EopSource::Finals. Unlike EOP_FILENAME this is not bundled -- it has to be
  // downloaded or supplied, since predictions go stale by construction.
  static constexpr std::string_view EOP_FINALS_FILENAME = "finals.all.iau1980.txt";
  static constexpr std::string_view IAU_SOFA_FILENAME = "IAU_SOFA.DAT";

  static constexpr std::string_view UNDEFINED = "Unknown";  // Used for unknown values

  // NAIF Intefer ID codes
  // Reference:
  // https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/naif_ids.html
  enum class BodyId {
    SSB = 0,
    SOLAR_SYSTEM_BARYCENTER = SSB,
    MERCURY_BARYCENTER = 1,
    VENUS_BARYCENTER = 2,
    EMB = 3,
    EARTH_MOON_BARYCENTER = EMB,
    MARS_BARYCENTER = 4,
    JUPITER_BARYCENTER = 5,
    SATURN_BARYCENTER = 6,
    URANUS_BARYCENTER = 7,
    NEPTUNE_BARYCENTER = 8,
    PLUTO_BARYCENTER = 9,
    SUN = 10,
    MERCURY = 199,
    VENUS = 299,
    EARTH = 399,
    MOON = 301,
    MARS = 499,
    PHOBOS = 401,
    DEIMOS = 402,
    JUPITER = 599,
    SATURN = 699,
    URANUS = 799,
    NEPTUNE = 899,
  };

  /// @brief Get the mean equatorial radius of a solar-system body [m, SI units].
  ///
  /// Used by occultation/eclipse checks (e.g. `lupnt/environment/occultation.cc`) and
  /// plotting helpers (`lupnt/interfaces/matplot.h`) that need a body's size to test
  /// line-of-sight blockage or to draw the body to scale. See also the
  /// `(BodyId, UnitSystem)` overload in `lupnt/environment/body.h` for results in a
  /// non-SI unit system.
  ///
  /// @param body Solar-system body identifier
  /// @return     Mean equatorial radius [m]
  double GetBodyRadius(BodyId body);

  enum class Time {
    UT1,  // Universal Time 1
    UTC,  // Coordinated Universal Time
    TAI,  // International Atomic Time
    TDB,  // Barycentric Dynamical Time
    TT,   // Terrestrial Time
    TCG,  // Geocentric Coordinate Time
    TCB,  // Barycentric Coordinate Time
    GPS,  // Global Positioning System Time
    // A Julian Date is a *representation* of an instant, not a time scale, so it
    // does not belong here. Use JdToTime()/TimeToJd() with the relevant scale.
    TCL,  // Lunar Coordinate Time
    LT,   // Lunar Time
  };

}  // namespace lupnt

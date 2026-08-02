/**
 * @file agent.h
 * @author Stanford NAV LAB
 * @brief List of agents
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <string>

#include "lupnt/conversions/frame_converter.h"
#include "lupnt/core/constants.h"

namespace lupnt {

  /**
   * @brief Physical and frame metadata for a solar-system body.
   */
  struct BodyData {
    BodyId id;
    std::string name;
    Real GM;
    Real R;
    Frame fixed_frame;
    Frame inertial_frame;
    Real flattening;
    Real omega;  // sidereal rotation rate
    UnitSystem units = SI_UNITS;
  };

  /**
   * @brief Spherical-harmonic gravity field coefficients and metadata.
   *
   * The gravitational parameter and reference radius are stored in the requested
   * UnitSystem. Coefficients are unnormalized after loading.
   */
  template <typename T = double> struct GravityField {
    int n_max, m_max;  // Maximum degree and order
    int n, m;          // Degree and order

    T GM;                            // Gravitational constant [m^3/s^2]
    T R;                             // Reference radius [m]
    Matrix<T, Dynamic, Dynamic> CS;  // Unnormalized coefficients
    UnitSystem units = SI_UNITS;
  };

  /**
   * @brief A celestial body used by environment and dynamics models.
   *
   * A Body stores the body's gravity parameter, reference radius, rotation rate,
   * fixed/inertial frames, optional spherical-harmonic gravity field, and unit
   * system. Factory methods create common solar-system bodies with SI units by
   * default or with a caller-provided UnitSystem.
   */
  struct Body {
    BodyId id;
    std::string name;

    Real GM;     // Gravitational constant
    Real R;      // Reference radius
    Real omega;  // sidereal rotation rate
    UnitSystem units = SI_UNITS;

    Frame fixed_frame;
    Frame inertial_frame;

    bool use_gravity_field;
    GravityField<Real> gravity_field;

    /// @brief Create the Sun as a point-mass `Body` (no spherical-harmonic gravity field), in SI
    /// units.
    static Body Sun();
    /// @brief Vec3-units overload of `Sun()`: create the Sun in the requested unit system.
    static Body Sun(const UnitSystem& units);

    /// @brief Create Earth as a `Body`, in SI units, optionally loading an `n`x`m`
    /// spherical-harmonic gravity field from `gravity_file`.
    ///
    /// Called by `Body::CreateBody`/`NBodyDynamics::AddBody` to assemble the
    /// central or perturbing-body list for `NumericalOrbitDynamics`; `n`/`m` set
    /// `GravityField::n_max`/`m_max` truncation used by `AccelarationGravityField`.
    ///
    /// @param n            Max degree of the gravity field to load (0 disables the field,
    /// point-mass only)
    /// @param m            Max order of the gravity field to load (0 disables the field, point-mass
    /// only)
    /// @param gravity_file Spherical-harmonic coefficient file name (default "EGM96.cof")
    /// @return Earth `Body` with GM, R, omega, frames, and (if requested) gravity field populated
    static Body Earth(int n = 0, int m = 0, std::string gravity_file = "EGM96.cof");
    /// @brief Unit-system overload of `Earth()`: create Earth in the requested unit system.
    static Body Earth(const UnitSystem& units, int n = 0, int m = 0,
                      std::string gravity_file = "EGM96.cof");

    /// @brief Create the Moon as a `Body`, in SI units, optionally loading an
    /// `n`x`m` spherical-harmonic gravity field (default GRGM1200B) from `gravity_file`.
    ///
    /// See `Earth()` for the role of the gravity-field parameters; used the same
    /// way by `Body::CreateBody`/`NBodyDynamics::AddBody` for Moon-centered or
    /// Moon-perturbation dynamics.
    static Body Moon(int n = 0, int m = 0, std::string gravity_file = "grgm1200b.cof");
    /// @brief Unit-system overload of `Moon()`: create the Moon in the requested unit system.
    static Body Moon(const UnitSystem& units, int n = 0, int m = 0,
                     std::string gravity_file = "grgm1200b.cof");

    /// @brief Create Venus as a `Body`, in SI units, optionally loading an `n`x`m`
    /// spherical-harmonic gravity field (default MGN75HSAAP) from `gravity_file`.
    static Body Venus(int n = 0, int m = 0, std::string gravity_file = "MGN75HSAAP.cof");
    /// @brief Unit-system overload of `Venus()`: create Venus in the requested unit system.
    static Body Venus(const UnitSystem& units, int n = 0, int m = 0,
                      std::string gravity_file = "MGN75HSAAP.cof");

    /// @brief Create Mars as a `Body`, in SI units, optionally loading an `n`x`m`
    /// spherical-harmonic gravity field (default GMM1) from `gravity_file`.
    static Body Mars(int n = 0, int m = 0, std::string gravity_file = "GMM1.cof");
    /// @brief Unit-system overload of `Mars()`: create Mars in the requested unit system.
    static Body Mars(const UnitSystem& units, int n = 0, int m = 0,
                     std::string gravity_file = "GMM1.cof");

    /// @brief Create Jupiter as a point-mass `Body` (no gravity field), in SI units.
    static Body Jupiter();
    /// @brief Unit-system overload of `Jupiter()`: create Jupiter in the requested unit system.
    static Body Jupiter(const UnitSystem& units);

    /// @brief Create Saturn as a point-mass `Body` (no gravity field), in SI units.
    static Body Saturn();
    /// @brief Unit-system overload of `Saturn()`: create Saturn in the requested unit system.
    static Body Saturn(const UnitSystem& units);

    /// @brief Create Uranus as a point-mass `Body` (no gravity field), in SI units.
    static Body Uranus();
    /// @brief Unit-system overload of `Uranus()`: create Uranus in the requested unit system.
    static Body Uranus(const UnitSystem& units);

    /// @brief Create Neptune as a point-mass `Body` (no gravity field), in SI units.
    static Body Neptune();
    /// @brief Unit-system overload of `Neptune()`: create Neptune in the requested unit system.
    static Body Neptune(const UnitSystem& units);

    /// @brief Dispatch to the matching per-planet factory (`Sun`/`Earth`/`Moon`/...)
    /// based on `body`'s `BodyId`, in SI units.
    ///
    /// Called by `NBodyDynamics::AddBody` (via `Body::CreateBody(body_id, units_, n,
    /// m, gravity_file)`) to build each perturbing/central body listed in a
    /// dynamics-model config without per-body switch statements at the call site.
    ///
    /// @param body         SPICE/LuPNT body identifier selecting which planet factory to call
    /// @param n            Max degree of the gravity field to load (0 disables the field)
    /// @param m            Max order of the gravity field to load (0 disables the field)
    /// @param gravity_file Spherical-harmonic coefficient file name (empty = use the per-body
    /// default)
    /// @return `Body` populated by the corresponding factory method
    static Body CreateBody(BodyId body, int n = 0, int m = 0, std::string gravity_file = "");
    /// @brief Unit-system overload of `CreateBody()`: create a body from its
    /// identifier in the requested unit system.
    static Body CreateBody(BodyId body, const UnitSystem& units, int n = 0, int m = 0,
                           std::string gravity_file = "");
  };

  /// @brief Load and unnormalize a spherical-harmonic gravity field of degree/order
  /// `n`x`m` from `filename`, in SI units.
  ///
  /// Used by `Body::Earth`/`Moon`/`Mars`/`Venus` (and directly in examples/tests) to
  /// populate `GravityField::CS`, `GM`, and `R` for high-fidelity gravity
  /// acceleration via `AccelarationGravityField`. Parses a GFC-style coefficient
  /// file (`POTFIELD`/`RECOEF` records); if `normalized` is true, coefficients are
  /// converted to unnormalized form on load.
  ///
  /// @param filename   Coefficient file name, resolved via the LuPNT data path (e.g. "EGM96.cof",
  /// "grgm1200b.cof")
  /// @param n          Max degree to retain (must be <= the field's stored n_max)
  /// @param m          Max order to retain (must be <= the field's stored m_max)
  /// @param normalized True if the file's stored coefficients are normalized and should be
  /// converted to unnormalized form
  /// @return           `GravityField<T>` with GM [m^3/s^2], R [m], and unnormalized CS coefficients
  /// up to degree/order `n`/`m`
  template <typename T> GravityField<T> ReadHarmonicGravityField(const std::string& filename, int n,
                                                                 int m, bool normalized);
  /// @brief Unit-system overload of `ReadHarmonicGravityField()`: load a
  /// spherical-harmonic gravity field and scale GM/R into the requested unit system.
  template <typename T> GravityField<T> ReadHarmonicGravityField(const std::string& filename, int n,
                                                                 int m, bool normalized,
                                                                 const UnitSystem& units);

  /// @brief Look up physical and frame metadata (GM, radius, flattening, sidereal
  /// rate, fixed/inertial frames) for a solar-system body, in SI units.
  ///
  /// Used throughout the simulator -- e.g. by `Occultation::ComputeOccultation` to
  /// get an occulting body's radius/name, by `GroundStation` and
  /// `JointOrbitClockDynamics::RelativisticRate` to get the central body's GM and
  /// frames, and by frame-conversion/plotting code -- as the single source of truth
  /// for body constants instead of hardcoded globals.
  ///
  /// @param id    Body identifier (e.g. `BodyId::EARTH`, `BodyId::MOON`)
  /// @return      `BodyData` with GM [m^3/s^2], R [m], flattening, omega [rad/s], and
  /// fixed/inertial `Frame`s
  BodyData GetBodyData(BodyId id);
  /// @brief Unit-system overload of `GetBodyData()`: return physical and frame
  /// metadata for a body, with GM/R/omega scaled into the requested unit system.
  BodyData GetBodyData(BodyId id, const UnitSystem& units);

  /// @brief Return the body's mean equatorial reference radius [m] (SI units).
  ///
  /// Convenience accessor over `GetBodyData` used e.g. by `Occultation` for
  /// eclipse/occultation radius checks and by plotting code (`maplot.cc`) to size
  /// drawn bodies.
  ///
  /// @param body Body identifier
  /// @return     Reference radius [m]
  double GetBodyRadius(BodyId body);
  /// @brief Unit-system overload of `GetBodyRadius()`: return the body's reference
  /// radius in the requested distance unit.
  double GetBodyRadius(BodyId body, const UnitSystem& units);

  /// @brief Return the body's gravitational parameter GM [m^3/s^2] (SI units).
  ///
  /// Convenience accessor over `GetBodyData` used e.g. by
  /// `time_conversions.cc` for relativistic time-scale corrections and by
  /// `StateConverter` to convert orbital-element states using the correct
  /// central-body GM.
  ///
  /// @param body Body identifier
  /// @return     Gravitational parameter [m^3/s^2]
  double GetBodyGM(BodyId body);
  /// @brief Unit-system overload of `GetBodyGM()`: return the body's gravitational
  /// parameter in the requested unit system.
  double GetBodyGM(BodyId body, const UnitSystem& units);

  /// @brief Return the body's sidereal rotation rate [rad/s] (SI units).
  /// @param body Body identifier
  /// @return     Sidereal rotation rate [rad/s]
  double GetBodyOmega(BodyId body);
  /// @brief Unit-system overload of `GetBodyOmega()`: return the body's sidereal
  /// rotation rate in the requested time unit.
  double GetBodyOmega(BodyId body, const UnitSystem& units);

  /// @brief Return the body's oblateness (flattening) coefficient [-] (SI units, dimensionless).
  double GetBodyFlattening(BodyId body);

  /// @brief Return the body's display name (e.g. "EARTH", "MOON").
  std::string GetBodyName(BodyId body);

  /// @brief Return the body's inertial reference `Frame` (e.g. GCRF for Earth, MOON_CI for the
  /// Moon).
  Frame GetInertialFrameName(BodyId body);

  /// @brief Return the body's body-fixed `Frame` (e.g. ITRF for Earth, MOON_PA for the Moon).
  Frame GetBodyFixedFrameName(BodyId body);

  /// @brief Free-function wrapper for `Body::CreateBody(body, n, m)`: create a body
  /// from its identifier in SI units.
  Body CreateBody(BodyId body, int n = 0, int m = 0);
  /// @brief Unit-system overload of `CreateBody()`: create a body from its
  /// identifier in the requested unit system.
  Body CreateBody(BodyId body, const UnitSystem& units, int n = 0, int m = 0);

}  // namespace lupnt

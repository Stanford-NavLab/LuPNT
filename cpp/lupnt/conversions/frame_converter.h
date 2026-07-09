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
#include <ostream>

#include "lupnt/core/constants.h"

namespace lupnt {

  class State;
  enum class Axis { X = 0, Y = 1, Z = 2 };

  enum class Frame {
    // Earth
    UNDEFINED,    // No frame
    ICRF,         // International Celestial Reference System
    ITRF,         // International Terrestrial Reference Frame
    ECEF = ITRF,  // Earth-Centered Earth-Fixed
    GCRF,         // Geocentric Reference System
    EME,          // Earth-Centered mean equator and equinox at J2000
    ECI = EME,    // Earth-Centered Inertial
    SER,          // Sun-Earth Rotating Frame
    GSE,          // Geocentric Solar Ecliptic
    MOD,          // Mean of date equatorial system
    TOD,          // True of date equatorial system
    EMR,          // Earth-Moon Rotating Frame
    // Moon
    MOON_CI,  // Moon-centered Inertial Frame (Axis aligened with ICRF)
    MOON_PA,  // Moon-Fixed with principal axes
    MOON_ME,  // Moon-Fixed with mean-Earth / polar axes
    MOON_OP,  // Earth Orbit Frame
    // Solar System
    MERCURY_FIXED,  // Mercury fixed frame
    VENUS_FIXED,    // Venus fixed frame
    MARS_FIXED,     // Mars fixed frame
    JUPITER_FIXED,  // Jupiter fixed frame
    SATURN_FIXED,   // Saturn fixed frame
    URANUS_FIXED,   // Uranus fixed frame
    NEPTUNE_FIXED,  // Neptune fixed frame
    // Inertial
    MERCURY_CI,  // Mercury-centered Inertial Frame
    VENUS_CI,    // Venus-centered Inertial Frame
    MARS_CI,     // Mars-centered Inertial Frame
    JUPITER_CI,  // Jupiter-centered Inertial Frame
    SATURN_CI,   // Saturn-centered Inertial Frame
    URANUS_CI,   // Uranus-centered Inertial Frame
    NEPTUNE_CI,  // Neptune-centered Inertial Frame
  };

  /// @brief Writes the short name of `frame` (e.g. "GCRF", "MOON_PA") to an output
  /// stream, for logging/debugging reference-frame values.
  std::ostream &operator<<(std::ostream &os, Frame frame);

  extern const std::map<Frame, BodyId> frame_centers;

  /// @brief Returns the natural central body (NAIF body ID) of a reference frame,
  /// e.g. `Frame::GCRF`/`Frame::ITRF` -> `BodyId::EARTH`, `Frame::MOON_CI` ->
  /// `BodyId::MOON`, `Frame::ICRF` -> `BodyId::SOLAR_SYSTEM_BARYCENTER`.
  ///
  /// Used wherever code needs to know which body a frame is centered on without
  /// hard-coding the mapping -- e.g. `NumericalOrbitDynamics` checks that a state's
  /// frame and its body-fixed frame share the same center, `JointOrbitClockDynamics`
  /// uses it to determine the gravitational parameter for the propagated body, and
  /// `Clock`/`StateConverter` use it to pick the relativity-correction or
  /// gravitational-parameter reference body for a state's frame.
  ///
  /// @param frame  Reference frame
  /// @return       NAIF body identifier for the frame's origin
  BodyId GetFrameCenter(Frame frame);

  /// @brief True if `frame` is a solar-system planet body-fixed frame (e.g.
  /// `MARS_FIXED`).
  bool IsPlanetFixedFrame(Frame frame);
  /// @brief True if `frame` is a solar-system planet-centered inertial frame
  /// (e.g. `MARS_CI`), i.e. ICRF-aligned axes centered on the planet.
  bool IsPlanetCiFrame(Frame frame);
  /// @brief True if `frame` is either a planet `*_CI` or `*_FIXED` frame.
  bool IsPlanetFrame(Frame frame);

  extern std::map<std::pair<Frame, Frame>, std::function<Vec6(Real, const Vec6 &rv)>>
      frame_conversions;

  /// @brief Converts a position vector at epoch `t_tdb` from `frame_in` to
  /// `frame_out`.
  ///
  /// This is the position-only (Vec3) overload of the central frame-conversion
  /// dispatcher used throughout the simulator (dynamics, measurements, agents,
  /// environment) to express a state in whatever frame the caller needs; it
  /// dispatches through `ConvertFrameBase`, which walks the frame graph (see the
  /// adjacency diagram on the Vec6 `ConvertFrame` overload below) applying the
  /// individual `*To*` conversions from frame_conversions.h until `frame_out` is
  /// reached.
  ///
  /// @param t_tdb        Epoch [s, TDB since J2000]
  /// @param r_in         Position vector in `frame_in` [m]
  /// @param frame_in     Input reference frame
  /// @param frame_out    Output reference frame
  /// @param rotate_only  If true, apply only the rotation (and skip the
  ///                      frame-origin translation) -- e.g. to transform a
  ///                      direction/offset vector rather than an absolute position
  /// @return              Position vector in `frame_out` [m]
  Vec3 ConvertFrame(Real t_tdb, const Vec3 &r_in, Frame frame_in, Frame frame_out,
                    bool rotate_only = false);

  /// @brief Converts a Cartesian position+velocity state at epoch `t_tdb` from
  /// `frame_in` to `frame_out`.
  ///
  /// This is the primary entry point used throughout the simulator (dynamics
  /// propagators, measurement models, agent state access, environment models)
  /// whenever a Cartesian state must be re-expressed in a different reference
  /// frame. It dispatches through `ConvertFrameBase`, which recursively walks the
  /// frame graph (ITRF/TOD/MOD/CIRS/EME/SER/GSE/EMR <-> GCRF <-> ICRF, and
  /// GCRF <-> MoonCi <-> MoonPa/MoonMe/MoonOp for the Moon) applying the
  /// individual `*To*` conversion functions declared in frame_conversions.h until
  /// `frame_out` is reached.
  ///
  /// @param t_tdb        Epoch [s, TDB since J2000]
  /// @param rv_in        Position+velocity state in `frame_in` [m, m/s]
  /// @param frame_in     Input reference frame
  /// @param frame_out    Output reference frame
  /// @param rotate_only  If true, apply only the rotation (and skip the
  ///                      frame-origin translation/velocity offset) -- e.g. to
  ///                      transform a relative state rather than an absolute one
  /// @return              Position+velocity state in `frame_out` [m, m/s]
  Vec6 ConvertFrame(Real t_tdb, const Vec6 &rv_in, Frame frame_in, Frame frame_out,
                    bool rotate_only = false);

  /// @brief Row-wise Cartesian-state overload of ConvertFrame(): converts each row
  /// (a 6-element position+velocity state) of `rv_in` at the single epoch `t_tdb`
  /// from `frame_in` to `frame_out`.
  MatX6 ConvertFrame(Real t_tdb, const MatX6 &rv_in, Frame frame_in, Frame frame_out,
                     bool rotate_only = false);

  /// @brief Row-wise position overload of ConvertFrame(): converts each row (a
  /// 3-element position vector) of `r_in` at the single epoch `t_tdb` from
  /// `frame_in` to `frame_out`.
  MatX3 ConvertFrame(Real t_tdb, const MatX3 &r_in, Frame frame_in, Frame frame_out,
                     bool rotate_only = false);

  /// @brief Multi-epoch overload of ConvertFrame(): converts the single
  /// Cartesian state `rv_in` from `frame_in` to `frame_out` independently at each
  /// epoch in `t_tdb`, returning one output row per epoch.
  MatX6 ConvertFrame(const VecX &t_tdb, const Vec6 &rv_in, Frame frame_in, Frame frame_out,
                     bool rotate_only = false);

  /// @brief Multi-epoch overload of ConvertFrame(): converts the single position
  /// vector `r_in` from `frame_in` to `frame_out` independently at each epoch in
  /// `t_tdb`, returning one output row per epoch.
  MatX3 ConvertFrame(const VecX &t_tdb, const Vec3 &r_in, Frame frame_in, Frame frame_out,
                     bool rotate_only = false);

  /// @brief Time-tagged-rows overload of ConvertFrame(): converts row `i` of
  /// `rv_in` from `frame_in` to `frame_out` at epoch `t_tdb(i)`, for each `i`
  /// (`t_tdb.size()` must equal `rv_in.rows()`).
  MatX6 ConvertFrame(const VecX &t_tdb, const MatX6 &rv_in, Frame frame_in, Frame frame_out,
                     bool rotate_only = false);

  /// @brief Time-tagged-rows overload of ConvertFrame(): converts row `i` of
  /// `r_in` from `frame_in` to `frame_out` at epoch `t_tdb(i)`, for each `i`
  /// (`t_tdb.size()` must equal `r_in.rows()`).
  MatX3 ConvertFrame(const VecX &t_tdb, const MatX3 &r_in, Frame frame_in, Frame frame_out,
                     bool rotate_only = false);

  /// @brief Converts a typed LuPNT `State` (e.g. `Cart6`, `ClassicalOE`, ...) to the
  /// frame `frame_out` at epoch `t_tdb`, preserving its element-type
  /// representation.
  ///
  /// Internally converts `x_in` to Cartesian (`Cart6`) via `ConvertState`, applies
  /// the Vec6 `ConvertFrame` overload to change frames, then converts back to
  /// `x_in`'s original `StateType`. Used by agent/dynamics code that stores a state
  /// in a non-Cartesian representation (e.g. classical orbital elements) but needs
  /// it in a different reference frame.
  ///
  /// @param t_tdb        Epoch [s, TDB since J2000]
  /// @param x_in         Input state (any `StateType`), tagged with its current `Frame`
  /// @param frame_out    Output reference frame
  /// @param rotate_only  If true, apply only the rotation and skip frame-origin translation
  /// @return              `x_in` re-expressed in `frame_out`, in its original `StateType`
  State ConvertFrame(Real t_tdb, const State &x_in, Frame frame_out, bool rotate_only = false);

  /// @brief Returns the affine position transform `(R, t)` such that
  /// `r_to = R * r_from + t` converts a position from `from_frame` to `to_frame` at
  /// epoch `t_tdb`.
  ///
  /// Used by the Vec3/Vec6 `ConvertFrame` overloads' `rotate_only` path, and by
  /// `NumericalOrbitDynamics` to obtain the rotation/translation between a
  /// propagated body's inertial frame and its body-fixed frame (e.g. for
  /// evaluating a body-fixed gravity field at an inertial-frame position) without
  /// repeatedly re-deriving it from the full frame-conversion chain.
  ///
  /// @param t_tdb        Epoch [s, TDB since J2000]
  /// @param from_frame   Input reference frame
  /// @param to_frame     Output reference frame
  /// @return              Pair `(R, t)`: 3x3 rotation matrix R (dimensionless) and
  ///                      3-vector translation t [m]
  std::pair<Mat3, Vec3> GetFrameRotationTranslation(Real t_tdb, Frame from_frame, Frame to_frame);

  /// @brief Returns the affine Cartesian-state transform `(R6, t6)` such that
  /// `x_to = R6 * x_from + t6` converts a 6-element position+velocity state from
  /// `from_frame` to `to_frame` at epoch `t_tdb`.
  ///
  /// `t6` is obtained by transforming the zero state, and `R6` (6x6) by
  /// transforming each of the 6 basis states and subtracting `t6`; this lets
  /// callers apply a single linear(-affine) map to many states/derivatives without
  /// repeatedly calling `ConvertFrame`. Used wherever a Jacobian or batch of
  /// relative states needs to be mapped between frames at a fixed epoch.
  ///
  /// @param t_tdb        Epoch [s, TDB since J2000]
  /// @param from_frame   Input reference frame
  /// @param to_frame     Output reference frame
  /// @return              Pair `(R6, t6)`: 6x6 rotation/state matrix R6
  ///                      (dimensionless for the rotation block) and 6-vector
  ///                      translation t6 [m, m/s]
  std::pair<Mat6, Vec6> GetFrameRotationTranslationRv(Real t_tdb, Frame from_frame, Frame to_frame);

}  // namespace lupnt

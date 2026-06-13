/**
 * @file spice_interface.h
 * @author Stanford NAV LAB
 * @brief  SPICE Interface functions
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <string.h>

#include <filesystem>
#include <map>
#include <string>
#include <tuple>

#include "lupnt/conversions/frame_converter.h"
#include "lupnt/core/constants.h"
#include "lupnt/interfaces/spice_cheby.h"

namespace lupnt {

  namespace spice {
    struct GroundStationSpiceData {
      int naif_id;
      std::string name;
      BodyId body_id = BodyId::EARTH;
      std::string frame;
      Vec3d position_m = Vec3d::Zero();
      double latitude_deg = 0.0;
      double longitude_deg = 0.0;
      double altitude_m = 0.0;
    };

    /// @brief Load LuPNT's default SPICE kernel set (leap seconds, planetary
    /// ephemeris, planetary constants, and high-accuracy Earth/lunar
    /// orientation kernels), and extract the Chebyshev position-only SPK
    /// segments for `GetBodyPosVel`.
    ///
    /// Called lazily (idempotent -- returns immediately if kernels are
    /// already loaded) by every other function in this namespace before it
    /// touches the SPICE kernel pool, e.g. `GetBodyPosVel`,
    /// `GetFrameConversionMat`, `GetGroundStationDataSpice`,
    /// `StringToTdb`/`StringToTai`. The high-accuracy Earth (`*.bpc`) and
    /// lunar (`moon_pa_de*.bpc` + `moon_de*.tf`) orientation kernels are the
    /// data source consumed by `ComputeEopFromSpice` and
    /// `InitFrameConversionFromSpice` (see frame_conversions.h) to calibrate
    /// LuPNT's analytic Earth/lunar orientation models. New high-accuracy
    /// kernels are downloaded from the NAIF generic-kernels archive (falling
    /// back to cached/bundled copies if offline, or if
    /// `LUPNT_SKIP_SPICE_KERNEL_DOWNLOAD` is set).
    void LoadSpiceKernel(void);

    /// @brief Load the default kernel set (if not already loaded), then load
    /// an additional SPICE kernel file from `filepath`.
    ///
    /// Used by `GroundStation::ConfigureFromData` (and similar agent-setup
    /// code) to load mission-specific SPICE kernels (e.g. spacecraft SPK/CK
    /// files or station-specific frame kernels) given in a simulation
    /// YAML config's `spice_kernel(s)` entries.
    void LoadSpiceKernel(const std::filesystem::path& filepath);

    /// @brief Experimental/debug helper that extracts PCK orientation
    /// (Euler-angle) coefficients for the Moon directly via low-level CSPICE
    /// PCK-reader calls (`pckeul_`, `pcklof_c`, `pcksfs_`, `pckr02_`).
    ///
    /// Not called anywhere in the simulation pipeline; retained as a
    /// reference/scratch routine for inspecting the structure of a binary PCK
    /// (hardcoded to `../data/ephemeris/moon_pa_de440_200625.bpc`, body 301 =
    /// Moon) and prints diagnostics to stdout.
    void ExtractPckCoeffs(void);

    /// @brief Resolve a body/spacecraft/station name to its NAIF integer ID
    /// code (e.g. `"MOON"` -> 301, `"EARTH"` -> 399).
    ///
    /// Used wherever a SPICE call needs an integer body code, e.g. by
    /// `GetGroundStationDataSpice(t_tdb, station_name, ...)` to resolve a
    /// named ground station to its NAIF ID before the position lookup.
    /// Throws if `body_name` is not found in the loaded kernel pool.
    int GetNaifId(const std::string& body_name);

    /// @brief Resolve a NAIF integer ID code to its body/spacecraft/station
    /// name (e.g. 301 -> `"MOON"`).
    ///
    /// Used by `GetGroundStationDataSpice(t_tdb, station_id, ...)` to fill in
    /// `GroundStationSpiceData::name` from a numeric station ID. Returns the
    /// numeric ID as a string if no name is found (does not throw).
    std::string GetNaifName(int naif_id);

    /// @brief True if `body_name` resolves to a NAIF body/spacecraft/station
    /// ID in the currently-loaded SPICE kernel pool.
    bool HasNaifBody(const std::string& body_name);

    /// @brief True if `frame_name` (e.g. `"ITRF93"`, `"IAU_MOON"`) is a
    /// recognized SPICE reference frame in the currently-loaded kernel pool.
    ///
    /// Used to check, before calling `GetFrameConversionMat`, whether a
    /// high-accuracy frame (e.g. `"ITRF93"`, `"IAU_MOON"`) is actually
    /// available from the loaded kernels.
    bool HasFrame(const std::string& frame_name);

    /// @brief Get the 6x6 state (position+velocity) rotation matrix from
    /// `from_frame` to `to_frame` at epoch `t_tdb`, via SPICE's `sxform_c`.
    ///
    /// This is the core building block of `spice::ConvertFrameSpice`
    /// (frame_converter_spice.cc), used e.g. to rotate ECEF <-> GCRF
    /// (`"ITRF93"` <-> `"J2000"`) and Moon-body-fixed <-> Moon-inertial
    /// (`"IAU_MOON"` <-> `"J2000"`) state vectors using SPICE's
    /// high-accuracy orientation kernels (see `LoadSpiceKernel`). Applying
    /// the returned matrix to a `Vec6` position+velocity state transforms it
    /// between frames, including the frame's angular-velocity (Coriolis)
    /// contribution to velocity.
    ///
    /// @param t_tdb       Epoch [s, TDB seconds since J2000]
    /// @param from_frame  SPICE name of the source frame (e.g. `"ITRF93"`)
    /// @param to_frame    SPICE name of the target frame (e.g. `"J2000"`)
    /// @return            6x6 state-transformation matrix `M` such that
    ///                    `rv_to = M * rv_from`
    Mat6d GetFrameConversionMat(Real t_tdb, const std::string& from_frame,
                                const std::string& to_frame);

    /// @brief Get a planet's orientation as the (alpha0, delta0, W) IAU
    /// rotation-pole right ascension/declination and prime-meridian angle
    /// [rad], derived from SPICE's `"IAU_<body>"` -> `"J2000"` rotation at
    /// `t_tdb`.
    ///
    /// Provides a SPICE-based cross-check / alternative source for the
    /// (alpha0, delta0, W) IAU body-orientation angles used elsewhere to
    /// build body-fixed <-> inertial rotation matrices for non-Earth/Moon
    /// bodies.
    ///
    /// @param id    Body whose orientation is requested (resolved to
    ///              `"IAU_<bodyname>"` via `bodc2n_c`)
    /// @param t_tdb Epoch [s, TDB seconds since J2000]
    /// @return      `(alpha0, delta0, W)` [rad]: right ascension and
    ///              declination of the body's rotation pole in J2000, and the
    ///              prime-meridian rotation angle
    Vec3d GetPlanetOrientation(BodyId id, Real t_tdb);

    /// @brief Parse a calendar/Julian-date string into ephemeris time (TDB
    /// seconds past J2000), via SPICE's `str2et_c`.
    ///
    /// Accepts the wide range of UTC/calendar/Julian-date string formats
    /// supported by SPICE (see the format table in spice.cc); used to convert
    /// human-readable epoch strings (e.g. from YAML scenario configs) into
    /// the `Real t_tdb` representation used throughout LuPNT.
    ///
    /// @param str Time string (UTC calendar, day-of-year, or Julian date)
    /// @return     Epoch [s, TDB seconds since J2000]
    Real StringToTdb(const std::string& str);

    /// @brief Parse a calendar/Julian-date string directly into TAI seconds
    /// (via `StringToTdb` + `ConvertTime(..., Time::TDB, Time::TAI)`).
    Real StringToTai(const std::string& str);

    /// @brief Format a TAI epoch as a UTC calendar string (`et2utc_c`,
    /// `"C"` format), via `TDBtoStringUTC`.
    ///
    /// @param t_tai TAI epoch [s, seconds since J2000]
    /// @param prec  Number of digits of fractional seconds in the output
    /// @return      UTC calendar string, e.g. `"2025 JUN 07 12:00:00.000"`
    std::string TAItoStringUTC(Real t_tai, int prec);

    /// @brief Format a TDB epoch as a UTC calendar string, via SPICE's
    /// `et2utc_c` with format `"C"`.
    ///
    /// Used by `TAItoStringUTC` and elsewhere for human-readable logging /
    /// debugging of simulation epochs.
    ///
    /// @param t_tdb TDB epoch [s, seconds since J2000]
    /// @param prec  Number of digits of fractional seconds in the output
    /// @return      UTC calendar string, e.g. `"2025 JUN 07 12:00:00.000"`
    std::string TDBtoStringUTC(Real t_tdb, int prec);

    /// @brief Convert a time value between time systems (TAI, TDB, TT, UTC,
    /// GPS, ...) using SPICE's `unitim_c`.
    ///
    /// A SPICE-based alternative to `lupnt::ConvertTime`
    /// (time_conversions.h); computes the offset between `from_time` and
    /// `to_time` at `t` via SPICE and applies it as a (autodiff-preserving)
    /// shift to `t`, so the returned `Real` retains `t`'s derivative
    /// information.
    ///
    /// @param t         Input time value [s]
    /// @param from_time Time system of `t`
    /// @param to_time   Desired output time system
    /// @return          `t` converted to `to_time` [s]
    Real ConvertTime(Real t, Time from_time, Time to_time);

    /// @brief Get the inertial (J2000) position and velocity of `target`
    /// relative to `center` at `t_tdb`, via the cached Chebyshev SPK segments
    /// extracted by `LoadSpiceKernel`.
    ///
    /// This is `lupnt`'s primary high-accuracy ephemeris source: it is the
    /// function backing `spice::ConvertFrameSpice`'s body-relative frame
    /// shifts (GCRF<->ICRF, GCRF<->MOON_CI, GCRF<->EMR, ...) and is the
    /// reference sampled by `InitFrameConversionFromSpice` /
    /// `ComputeEopFromSpice`. Handles barycenter-vs-body indirection
    /// internally (e.g. Earth/Moon/Mercury/Venus positions are obtained via
    /// their system barycenters plus a body-relative-to-barycenter Chebyshev
    /// segment) by delegating to the internal `GetBodyPosVelBase` helper.
    /// Note: returned in the SPK's native km / km/s units (unlike most other
    /// `lupnt` state vectors, which are in meters); `ConvertFrameSpice`
    /// operates consistently in these units internally.
    ///
    /// @param t_tdb   Epoch [s, TDB seconds since J2000]
    /// @param center  Origin body of the returned state
    /// @param target  Body whose state is returned
    /// @return        `[r; v]` position [km] and velocity [km/s] of `target`
    ///                relative to `center`, in the J2000 inertial frame
    Vec6 GetBodyPosVel(const Real t_tdb, BodyId center, BodyId target);

    /// @brief Vectorized (per-epoch) overload of `GetBodyPosVel`.
    ///
    /// @param t_tdb  Epochs [s, TDB seconds since J2000]
    /// @param center Origin body of the returned states
    /// @param target Body whose state is returned
    /// @return       `[N x 6]` position [km] / velocity [km/s] of `target`
    ///               relative to `center` at each epoch, in J2000
    MatX6 GetBodyPosVel(const VecX& t_tdb, BodyId center, BodyId target);

    /// @brief Get the position of `target` relative to `obs` at `t_tdb`
    /// directly from SPICE (`spkpos_c`), without LuPNT's Chebyshev cache.
    ///
    /// A lower-level, double-precision alternative to `GetBodyPosVel`,
    /// allowing arbitrary SPICE reference frames and aberration corrections
    /// (e.g. light-time/stellar-aberration corrections for measurement
    /// modeling) not expressible via the J2000-only `GetBodyPosVel`.
    ///
    /// @param t_tdb        Epoch [s, TDB seconds since J2000]
    /// @param obs          Observer body (SPICE `obs`)
    /// @param target       Target body (SPICE `targ`)
    /// @param refFrame     SPICE reference frame for the returned vector
    ///                     (default `"J2000"`)
    /// @param abCorrection SPICE aberration correction (e.g. `"NONE"`, `"LT"`)
    /// @return             Position of `target` relative to `obs` [km], in
    ///                     `refFrame`
    Vec3d GetBodyPosSpice(Real t_tdb, BodyId obs, BodyId target,
                          const std::string& refFrame = "J2000",
                          const std::string& abCorrection = "NONE");

    /// @brief Get the position and velocity of `target` relative to `obs` at
    /// `t_tdb` directly from SPICE (`spkezr_c`), without LuPNT's Chebyshev
    /// cache.
    ///
    /// Vec6 (position+velocity) counterpart of `GetBodyPosSpice`; see that
    /// function for the role of `refFrame`/`abCorrection`.
    ///
    /// @return `[r; v]` position [km] / velocity [km/s] of `target` relative
    ///         to `obs`, in `refFrame`
    Vec6d GetBodyPosVelSpice(Real t_tdb, BodyId obs, BodyId target,
                             const std::string& refFrame = "J2000",
                             const std::string& abCorrection = "NONE");

    /// @brief Look up a ground station's position and geodetic coordinates
    /// from SPICE by station name.
    ///
    /// Resolves `station_name` to its NAIF ID via `GetNaifId`, then delegates
    /// to the `station_id` overload. Used by `GroundStation`
    /// (agents/ground_station.h/.cc) to initialize a ground station's ECEF
    /// position and latitude/longitude/altitude from a SPICE station-frame
    /// kernel (e.g. DSN complexes) given a station name in a YAML config.
    ///
    /// @param t_tdb     Epoch [s, TDB seconds since J2000]
    /// @param station_name SPICE body/station name (e.g. a DSN station name)
    /// @param center    Center body for the returned position (default Earth)
    /// @param refFrame  SPICE reference frame for the position (default
    ///                  `"ITRF93"`, i.e. ECEF)
    /// @param abCorrection SPICE aberration correction (default `"NONE"`)
    /// @return          Station NAIF ID/name, body-fixed position [m], and
    ///                  geodetic latitude/longitude [deg] / altitude [m]
    GroundStationSpiceData GetGroundStationDataSpice(Real t_tdb, const std::string& station_name,
                                                     BodyId center = BodyId::EARTH,
                                                     const std::string& refFrame = "ITRF93",
                                                     const std::string& abCorrection = "NONE");

    /// @brief Look up a ground station's position and geodetic coordinates
    /// from SPICE by NAIF station ID.
    ///
    /// Computes the station's position via `spkpos_c` in `refFrame` relative
    /// to `center`, then converts it to geodetic latitude/longitude/altitude
    /// using `center`'s body shape (`GetBodyData`/`CartToLatLonAlt`). This is
    /// the underlying implementation used by the `station_name` overload, and
    /// by `GroundStation::ConfigureFromData` / config-based ground-station
    /// setup.
    ///
    /// @param t_tdb     Epoch [s, TDB seconds since J2000]
    /// @param station_id NAIF integer ID of the station
    /// @param center    Center body for the returned position (default Earth)
    /// @param refFrame  SPICE reference frame for the position (default
    ///                  `"ITRF93"`, i.e. ECEF)
    /// @param abCorrection SPICE aberration correction (default `"NONE"`)
    /// @return          Station NAIF ID/name, body-fixed position [m], and
    ///                  geodetic latitude/longitude [deg] / altitude [m]
    GroundStationSpiceData GetGroundStationDataSpice(Real t_tdb, int station_id,
                                                     BodyId center = BodyId::EARTH,
                                                     const std::string& refFrame = "ITRF93",
                                                     const std::string& abCorrection = "NONE");

  }  // namespace spice
}  // namespace lupnt

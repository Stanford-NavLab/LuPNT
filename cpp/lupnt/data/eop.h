#pragma once

#include <filesystem>

#include "lupnt/core/definitions.h"

namespace lupnt {

  // Struct to hold EOP data
  struct EopFileData {
    VecXi years;
    VecXi months;
    VecXi days;
    VecXd mjds_utc;
    VecXd x;
    VecXd y;
    VecXd ut1_utc;
    VecXd lod;
    VecXd dpsi;
    VecXd deps;
    VecXd xErr;
    VecXd yErr;
    VecXd ut1_utc_err;
    VecXd lod_err;
    VecXd dpsi_err;
    VecXd deps_err;
  };

  struct EopData {
    Real x_pole;
    Real y_pole;
    Real ut1_utc;
    Real lod;
    Real dpsi;
    Real deps;
    Real dx_pole;
    Real dy_pole;
    Real tai_utc;
  };

  /// @brief Parse an IERS EOP 14 C04 (IAU1980) time-series file into the global EOP table.
  ///
  /// Called lazily (via GetEopFileData/GetEopData) the first time Earth orientation data is
  /// needed, using the bundled file at GetFilePath(EOP_FILENAME); also called by
  /// LoadLatestEopFromIers after downloading a fresh file. Subsequent calls are no-ops unless
  /// `force` is set, since all later frame conversions (GcrfToItrf/ItrfToGcrf via
  /// RotPolarMotion/RotSideralMotion) and time conversions (UT1-UTC) read from this table.
  ///
  /// @param filepath Path to an IERS EOP 14 C04 (IAU1980) text file (14 header lines followed by
  ///                  whitespace-separated columns: year, month, day, MJD (UTC), x_pole, y_pole,
  ///                  UT1-UTC, LOD, dpsi, deps, and their formal errors).
  /// @param force    If true, reload and replace any previously-loaded EOP table even if data is
  ///                  already present (default: false).
  void LoadEopFileData(const std::filesystem::path& filepath, bool force = false);

  /// @brief Download the latest IERS EOP 14 C04 (IAU1980) time series and load it, replacing any
  /// previously-loaded EOP data.
  ///
  /// Used at simulation startup (or on demand) to refresh the polar-motion/UT1-UTC table used by
  /// GetEopData/GetUt1UtcDifference with up-to-date IERS measurements instead of the bundled
  /// snapshot. Caches the downloaded file under
  /// GetDataPath()/"planetary_coeff"/"eopc04_iers_latest.txt". On download failure, ensures the
  /// bundled EOP file is loaded (if not already) and returns false.
  ///
  /// @param force If true, reload even if EOP data has already been loaded (default: true).
  /// @return      True if the latest IERS data was downloaded and loaded successfully, false if
  ///              the download failed (in which case the bundled data is used as a fallback).
  bool LoadLatestEopFromIers(bool force = true);

  /// @brief Interpolate Earth orientation parameters (polar motion, UT1-UTC, LOD, nutation
  /// corrections) to a requested epoch from the loaded IERS EOP table.
  ///
  /// Called by frame-conversion routines such as RotPolarMotion, RotSideralMotion, and
  /// RotSideralMotionDot (see frame_conversions.cc) to obtain the EOP values needed for the
  /// GCRF<->ITRF rotation at a given epoch. Lazily loads the bundled EOP file via
  /// LoadEopFileData/GetFilePath(EOP_FILENAME) on first use. For epochs outside the table's
  /// range, the nearest endpoint's values are held constant (no extrapolation); otherwise a
  /// 3rd-order Lagrange interpolation is used.
  ///
  /// @param mjd_utc Epoch at which to evaluate the EOP [MJD, UTC]
  /// @return        Interpolated EOP values (x_pole/y_pole [rad], ut1_utc [s], lod [s],
  ///                dpsi/deps [rad], dx_pole/dy_pole [rad], tai_utc [s])
  EopData GetEopData(Real mjd_utc);

  /// @brief Return a pointer to the raw loaded IERS EOP table, loading the bundled file first if
  /// no EOP data has been loaded yet.
  ///
  /// Used where direct access to the full raw time series (e.g. for plotting or diagnostics) is
  /// needed, rather than the interpolated single-epoch values returned by GetEopData.
  ///
  /// @return Pointer to the global EopFileData table (owned internally; non-null after the call).
  EopFileData* GetEopFileData();

  /// @brief Get the UT1-UTC time difference at a given epoch.
  ///
  /// Used by time-conversion routines (see conversions/time_conversions.cc) to convert between
  /// UTC and UT1 (e.g. for Earth rotation angle / sidereal time computations). Thin wrapper
  /// around GetEopData that extracts just the ut1_utc field.
  ///
  /// @param mjd_utc Epoch at which to evaluate UT1-UTC [MJD, UTC]
  /// @return        UT1 - UTC [s]
  Real GetUt1UtcDifference(Real mjd_utc);

}  // namespace lupnt

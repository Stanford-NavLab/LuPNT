#pragma once

#include <filesystem>

#include "lupnt/core/definitions.h"

namespace lupnt {

  // Struct to hold EOP data
  struct IauSofaFileData {
    VecXd jd_tt;
    VecXd X;
    VecXd Y;
    VecXd s;
  };

  struct IauSofaData {
    Real X;
    Real Y;
    Real s;
  };

  /// @brief Parse a tabulated IAU SOFA CIP coordinates / CIO locator file (IAU_SOFA.DAT) into the
  /// global precession-nutation series table.
  ///
  /// Called lazily by GetIauSofaData (via GetFilePath(IAU_SOFA_FILENAME)) the first time the
  /// IAU 2006/2000A precession-nutation matrix is needed. Subsequent calls are no-ops once data
  /// is loaded.
  ///
  /// @param filepath Path to a whitespace-separated text file with columns: Julian Date (TT),
  ///                  CIP coordinate X [arcsec], CIP coordinate Y [arcsec], CIO locator s
  ///                  [arcsec], one row per day.
  void LoadIauSofaFileData(const std::filesystem::path& filepath);

  /// @brief Interpolate the CIP coordinates (X, Y) and CIO locator (s) to a requested epoch from
  /// the loaded IAU SOFA series table.
  ///
  /// Called by RotPrecessionNutation (frame_conversions.cc) to build the IAU 2006/2000A
  /// precession-nutation rotation matrix used throughout GCRF<->ITRF frame conversions
  /// (GcrfToItrf/ItrfToGcrf). Lazily loads the bundled table via LoadIauSofaFileData on first
  /// use. For epochs outside the table's range, the nearest endpoint's values are held constant;
  /// otherwise a 9th-order Lagrange interpolation is used.
  ///
  /// @param jd_tt Epoch at which to evaluate the series [Julian Date, TT]
  /// @return      Interpolated CIP coordinates and CIO locator, X/Y/s [arcsec] (multiply by
  ///              RAD_ARCSEC to convert to radians, as done in RotPrecessionNutation)
  IauSofaData GetIauSofaData(Real jd_tt);

}  // namespace lupnt

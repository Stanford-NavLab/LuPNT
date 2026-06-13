#include "lupnt/core/definitions.h"

namespace lupnt {
  struct TaiUtcFileData {
    VecXd jd;
    VecXd tai_utc;
    VecXd mjd0;
    VecXd scale;
  };

  /// @brief Get the TAI-UTC leap-second offset at a given epoch.
  ///
  /// Used by time-conversion routines (e.g. UtcToTai/TaiToUtc in conversions/time_conversions.cc)
  /// to convert between UTC and TAI, which underlies all higher time-scale conversions (TT, TDB,
  /// TCB, GPS time) used throughout the simulator. Lazily loads the bundled leap-second table
  /// (GetFilePath(TAI_UTC_FILENAME)) on first use. For epochs before the table's first entry,
  /// returns 0; for epochs at or after the last entry, extrapolates using that entry's linear
  /// drift rate; otherwise evaluates the piecewise-linear/step formula of the matching table
  /// entry.
  ///
  /// @param mjd Epoch at which to evaluate the offset [MJD, UTC]
  /// @return    TAI - UTC [s]
  double GetTaiUtcDifference(double mjd);
}  // namespace lupnt

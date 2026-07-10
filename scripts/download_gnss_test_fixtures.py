"""Download the SP3/BRDC GNSS fixture files the C++ test suite expects on disk.

`pixi run test-cpp` includes several tests (interfaces.sp3_loader,
interfaces.rinex_nav_loader, agents.gnss_constellation.setup_from_files,
measurements.gnss_measurements.visibility) that parse pre-downloaded CDDIS
products rather than downloading them itself -- see the `GnssFilesDir()`
helper in cpp/test/interfaces/test_sp3_loader.cc, which resolves to
`<repo>/output/gnss_files`. This script performs that one-time download using
LuPNT's C++ loaders (`pnt.Sp3Loader.download_file_for_epoch` /
`pnt.RinexNavLoader.download_file_for_epoch`), whose default cache directory
(`GetOutputDir("gnss_files")`) is exactly that location, for the epochs the C++
tests hardcode (2026-01-14 and the following day).

Requires a NASA Earthdata Login configured in ~/.netrc (see README.md,
"Prerequisites") or the EARTHDATA_USERNAME / EARTHDATA_PASSWORD environment
variables.
"""

from __future__ import annotations

import pylupnt as pnt


def _epoch_tai(year: int, month: int, day: int) -> float:
    """Mid-day TAI epoch for a calendar day (avoids midnight day-rollover ambiguity)."""
    t_tdb = pnt.gregorian_to_time(year, month, day, 12, 0, 0)
    return pnt.convert_time(t_tdb, pnt.Time.TDB, pnt.Time.TAI)


def main() -> None:
    # The C++ tests parse the 2026-01-14 product plus the following day (the 1-day
    # SP3 span in test_sp3_loader.cc spills into 2026-01-15).
    print("SP3 files:")
    for day in (14, 15):
        path = pnt.Sp3Loader.download_file_for_epoch(_epoch_tai(2026, 1, day), pnt.Time.TAI)
        print(f"  {path}")

    print("BRDC files:")
    path = pnt.RinexNavLoader.download_file_for_epoch(_epoch_tai(2026, 1, 14), pnt.Time.TAI)
    print(f"  {path}")


if __name__ == "__main__":
    main()

"""Download the SP3/BRDC GNSS fixture files the C++ test suite expects on disk.

`pixi run test-cpp` includes several tests (interfaces.sp3_loader,
interfaces.rinex_nav_loader, agents.gnss_constellation.setup_from_files,
measurements.gnss_measurements.visibility) that parse pre-downloaded CDDIS
products rather than downloading them itself -- see the `GnssFilesDir()`
helper in cpp/test/interfaces/test_sp3_loader.cc. This script performs that
one-time download via the same Python SP3Loader/BRDCLoader used by
python/examples/ex_gnss/download_sp3_epochs.py, for the exact epoch the C++
tests hardcode (2026-01-14/15).

Requires a NASA Earthdata Login configured in ~/.netrc (see README.md,
"Prerequisites").
"""

from __future__ import annotations

from datetime import datetime, timezone

import pylupnt as pnt
from pylupnt.interfaces.gnss_file_loader import BRDCLoader, SP3Loader


def main() -> None:
    epoch_utc = datetime(2026, 1, 14, 0, 0, 0, tzinfo=timezone.utc)

    sp3 = SP3Loader(target_dt=epoch_utc, sim_t=86400, dt_timesys=pnt.UTC)
    print("SP3 files:")
    for filename in sp3.filenames:
        print(f"  {filename}")

    brdc = BRDCLoader(target_dt=epoch_utc, sim_t=0, dt_timesys=pnt.UTC)
    print("BRDC files:")
    for filename in brdc.filenames:
        print(f"  {filename}")


if __name__ == "__main__":
    main()

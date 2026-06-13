"""Download and verify SP3 precise-orbit files for several epochs.

This example uses the Python SP3Loader, which downloads missing files from the
NASA CDDIS archive and stores them under ``output/gnss_files/sp3``. Configure
Earthdata Login credentials before running if CDDIS asks for authentication.
"""

from __future__ import annotations

from datetime import datetime, timezone

import pylupnt as pnt
from pylupnt.interfaces.gnss_file_loader import SP3Loader


def main() -> None:
    epochs_utc = [
        datetime(2025, 1, 1, 0, 0, 0, tzinfo=timezone.utc),
        datetime(2025, 1, 2, 0, 0, 0, tzinfo=timezone.utc),
        datetime(2025, 1, 3, 12, 0, 0, tzinfo=timezone.utc),
    ]

    for epoch_utc in epochs_utc:
        print(f"\n=== {epoch_utc.isoformat()} ===")

        # The Python interface loader handles CDDIS downloads.  The C++ binding
        # `pnt.Sp3Loader` is the parser/interpolator for files that already exist.
        sp3 = SP3Loader(target_dt=epoch_utc, sim_t=0, dt_timesys=pnt.UTC)
        print("Downloaded/loaded files:")
        for filename in sp3.filenames:
            print(f"  {filename}")

        # Parsed metadata confirms the downloaded file can be used by LuPNT.
        print(f"Number of satellites: {len(sp3.sats)}")
        print(f"First satellites: {', '.join(sp3.sats[:8])}")
        if len(sp3.epochs) > 0:
            print(f"Epoch coverage: {sp3.epochs[0]} to {sp3.epochs[-1]} TAI seconds")


if __name__ == "__main__":
    main()

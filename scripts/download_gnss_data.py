"""Download the SP3/BRDC GNSS products the examples and scenarios expect on disk.

`configs/lunar_gnss_odts.yaml` and `python/examples/ex6_gnss_odts*` read precise
orbits (SP3) and broadcast ephemerides (BRDC) from

    data/LuPNT_data/gnss/sp3
    data/LuPNT_data/gnss/brdc

alongside the other GNSS inputs already there (`igs20.atx`, `gps_table.csv`).
Neither directory ships with the repo -- `data/` is gitignored -- so this script
performs the one-time fetch using LuPNT's own loaders, which pull the COD MGEX
final products from CDDIS.

This is separate from `download_gnss_test_fixtures.py`, which populates
`output/gnss_files` for the C++ test suite at its own hardcoded epochs.

Requires a NASA Earthdata Login in ~/.netrc (see README.md, "Prerequisites") or
the EARTHDATA_USERNAME / EARTHDATA_PASSWORD environment variables.

Usage
-----
    pixi run python scripts/download_gnss_data.py
    pixi run python scripts/download_gnss_data.py --epoch 2026-01-01T02:00:00 --days 2
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pylupnt as pnt
import yaml

REPO = Path(__file__).resolve().parents[1]
GNSS_DIR = REPO / "data" / "LuPNT_data" / "gnss"
SP3_DIR = GNSS_DIR / "sp3"
BRDC_DIR = GNSS_DIR / "brdc"
CONFIG_YAML = REPO / "configs" / "lunar_gnss_odts.yaml"


def _epoch_from_config() -> str:
    """The scenario epoch, so the default download always matches what runs."""
    scen = yaml.safe_load(open(CONFIG_YAML))
    return scen["epoch"]


def _tai_seconds(utc_iso: str, day_offset: int = 0) -> float:
    """Mid-day TAI epoch for the calendar day `day_offset` days after `utc_iso`.

    Products are daily; anchoring at 12:00 avoids the midnight rollover picking
    the neighbouring file.
    """
    e = pnt.Epoch.from_gregorian(utc_iso, pnt.Time.UTC)
    day_start = pnt.Epoch.from_seconds(float(e.to_seconds()) + day_offset * 86400.0, pnt.Time.UTC)
    y, mo, d = (int(x) for x in day_start.to_gregorian_string(0).split("T")[0].split("-"))
    return float(pnt.convert_time(pnt.gregorian_to_time(y, mo, d, 12, 0, 0), pnt.UTC, pnt.TAI))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--epoch", default=None, help="UTC ISO epoch (default: configs YAML)")
    ap.add_argument(
        "--days",
        type=int,
        default=2,
        help="Number of consecutive daily products to fetch (default 2: a run "
        "starting mid-day can span into the following product).",
    )
    args = ap.parse_args()

    epoch = args.epoch or _epoch_from_config()
    SP3_DIR.mkdir(parents=True, exist_ok=True)
    BRDC_DIR.mkdir(parents=True, exist_ok=True)
    print(f"epoch    : {epoch}  (+{args.days} daily products)")
    print(f"sp3 dir  : {SP3_DIR}")
    print(f"brdc dir : {BRDC_DIR}\n")

    failed = False
    for day in range(args.days):
        t_tai = _tai_seconds(epoch, day)
        for label, fn, dest in (
            ("SP3 ", pnt.Sp3Loader.download_file_for_epoch, SP3_DIR),
            ("BRDC", pnt.RinexNavLoader.download_file_for_epoch, BRDC_DIR),
        ):
            try:
                path = fn(t_tai, pnt.TAI, str(dest))
                print(f"  {label} day+{day}: {Path(path).name}")
            except Exception as exc:  # noqa: BLE001 - report and continue
                failed = True
                print(f"  {label} day+{day}: FAILED - {exc}", file=sys.stderr)

    if failed:
        print(
            "\nOne or more downloads failed. CDDIS requires a NASA Earthdata Login:\n"
            "  * put your credentials in ~/.netrc, or\n"
            "  * set EARTHDATA_USERNAME / EARTHDATA_PASSWORD\n"
            "See README.md, 'Prerequisites'.",
            file=sys.stderr,
        )
        return 1

    print("\nDone. ex6_gnss_odts and configs/lunar_gnss_odts.yaml read these directly.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

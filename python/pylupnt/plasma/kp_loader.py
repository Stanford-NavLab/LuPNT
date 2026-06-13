"""
Kp-index loader for the plasma module.

Downloads and converts Kp index data from GFZ Potsdam.
Data is cached under $PECSIMPY_BASE_PATH/data/kp/ (set in pixi activation env
to $LUPNT_DATA_PATH/plasma/data/kp/).
"""

import json
import os
from datetime import datetime

import numpy as np
import requests


def _get_base_dir():
    """Return the plasma data base directory.

    Resolution order:

    #. ``PECSIMPY_BASE_PATH`` environment variable (set by pixi activation).
    #. ``LUPNT_DATA_PATH`` + ``/plasma`` (fallback when the notebook kernel
       does not inherit the full pixi activation env).
    """
    base = os.environ.get("PECSIMPY_BASE_PATH", "")
    if base:
        return base
    lupnt = os.environ.get("LUPNT_DATA_PATH", "")
    if lupnt:
        candidate = os.path.join(lupnt, "plasma")
        if os.path.isdir(candidate):
            return candidate
    raise RuntimeError(
        "Could not locate the plasma data directory. "
        "Set PECSIMPY_BASE_PATH (e.g. $LUPNT_DATA_PATH/plasma) or activate "
        "the pixi environment with 'pixi shell'."
    )


def load_table(base_dir, end_yyyy, end_mm, end_dd):
    """Download Kp index JSON for *end_yyyy* if not already up to date."""
    kp_json_dir = os.path.join(base_dir, "data", "kp", "json")
    os.makedirs(kp_json_dir, exist_ok=True)

    first_date_file = os.path.join(kp_json_dir, "first_date.json")
    last_date_file = os.path.join(kp_json_dir, "last_date.json")
    _sentinel = {"yyyy": 0, "mm": 0, "dd": 0}

    if not os.path.exists(first_date_file):
        with open(first_date_file, "w") as f:
            json.dump(_sentinel, f)
    if not os.path.exists(last_date_file):
        with open(last_date_file, "w") as f:
            json.dump(_sentinel, f)

    with open(first_date_file) as f:
        first_date = json.load(f)
    with open(last_date_file) as f:
        last_date = json.load(f)

    update_first_date = (end_yyyy, 1, 1) < (
        first_date["yyyy"],
        first_date["mm"],
        first_date["dd"],
    )
    update_last_date = (last_date["yyyy"], last_date["mm"], last_date["dd"]) < (
        end_yyyy,
        end_mm,
        end_dd,
    )

    if not (update_first_date or update_last_date):
        print(f"Kp index already up to date (through {end_yyyy}-{end_mm:02d}-{end_dd:02d}).")
        return

    url = (
        f"https://kp.gfz-potsdam.de/app/json/"
        f"?start={end_yyyy}-01-01T00%3A00%3A00Z"
        f"&end={end_yyyy}-{end_mm:02d}-{end_dd:02d}T23%3A59%3A59Z"
        f"&index=Kp#kpdatadownload-143"
    )

    response = requests.get(url)
    if response.status_code == 200:
        target = os.path.join(kp_json_dir, f"kp_{end_yyyy}.json")
        with open(target, "wb") as fh:
            fh.write(response.content)
        print(f"Downloaded Kp index for {end_yyyy}.")
    else:
        print(f"Download failed (HTTP {response.status_code}).")

    if update_last_date:
        last_date.update({"yyyy": end_yyyy, "mm": end_mm, "dd": end_dd})
        with open(last_date_file, "w") as f:
            json.dump(last_date, f)

    if update_first_date:
        first_date.update({"yyyy": end_yyyy, "mm": 1, "dd": 1})
        with open(first_date_file, "w") as f:
            json.dump(first_date, f)


def update_kp_table(base_dir=None, start_year: int = 1995):
    """Download Kp JSON files from *start_year* to today."""
    if not base_dir:
        base_dir = _get_base_dir()
    now = datetime.now()
    for year in range(start_year, now.year):
        load_table(base_dir, year, 12, 31)
    load_table(base_dir, now.year, now.month, now.day)


def convert_to_csv(base_dir=None):
    """Convert all kp_YYYY.json files to kp_YYYY.csv."""
    if not base_dir:
        base_dir = _get_base_dir()
    json_dir = os.path.join(base_dir, "data", "kp", "json")
    csv_dir = os.path.join(base_dir, "data", "kp", "csv")
    os.makedirs(csv_dir, exist_ok=True)

    for fname in os.listdir(json_dir):
        if not (fname.endswith(".json") and fname.startswith("kp_")):
            continue
        with open(os.path.join(json_dir, fname)) as f:
            data = json.load(f)

        datetimes = data["datetime"]
        kp_vals = data["Kp"]
        rows = []
        for dt_str, kp in zip(datetimes, kp_vals):
            date_part, time_part = dt_str.split("T")
            y, m, d = date_part.split("-")
            h = time_part.split(":")[0]
            rows.append(f"{y},{m},{d},{h},{kp}")

        csv_path = os.path.join(csv_dir, fname.replace(".json", ".csv"))
        with open(csv_path, "w") as f:
            f.write("yyyy,mm,dd,hh,kp\n")
            f.write("\n".join(rows) + "\n")

    print("Converted Kp JSON files to CSV.")


def update_kp(base_dir=None, start_year: int = 1995):
    """Download Kp tables and convert to CSV in one step.

    Parameters
    ----------
    base_dir:
        Plasma data base directory.  Defaults to ``$PECSIMPY_BASE_PATH``.
        An empty string is treated the same as ``None``.
    start_year:
        Earliest year to fetch (default 1995).
    """
    if not base_dir:
        base_dir = _get_base_dir()
    print(f"Plasma data directory: {base_dir}")
    update_kp_table(base_dir, start_year)
    convert_to_csv(base_dir)
    print("Kp index updated and converted to CSV successfully.")

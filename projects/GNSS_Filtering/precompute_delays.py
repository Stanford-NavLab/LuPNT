#!/usr/bin/env python3
from __future__ import annotations

import argparse
import multiprocessing as mp
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


FREQ_HZ = {
    "L1": 1575.42e6,
    "L2": 1227.60e6,
    "L5": 1176.45e6,
    "E1": 1575.42e6,
    "E6": 1278.75e6,
    "E5": 1191.795e6,
    "E5a": 1176.45e6,
    "E5b": 1207.14e6,
}

EARTH_RADIUS_KM = 6371.2

IDENTITY_INT_COLUMNS = ["link_id", "epoch_index", "gnss_const_id", "prn", "frequency_id"]
IDENTITY_FLOAT_COLUMNS = [
    "t_tdb",
    "rx_x_ecef_m",
    "rx_y_ecef_m",
    "rx_z_ecef_m",
    "tx_x_ecef_m",
    "tx_y_ecef_m",
    "tx_z_ecef_m",
]


def tangent_radius_altitude_m(tx_m: np.ndarray, rx_m: np.ndarray) -> tuple[float, float]:
    line = tx_m - rx_m
    line_norm2 = float(np.dot(line, line))
    u = 0.0 if line_norm2 <= 0.0 else float(np.clip(-np.dot(rx_m, line) / line_norm2, 0.0, 1.0))
    radius_m = float(np.linalg.norm(rx_m + u * line))
    return radius_m, radius_m - EARTH_RADIUS_KM * 1000.0


def _resolve(path: str | None, fallback: Path) -> Path:
    if not path:
        return fallback.resolve()
    p = Path(path)
    return p if p.is_absolute() else p.resolve()


def load_pipeline_config(config_path: Path) -> tuple[Path, Path, dict]:
    with config_path.open() as f:
        cfg = yaml.safe_load(f)
    output_dir = _resolve(
        cfg.get("simulation", {}).get("output_dir"), Path("output/gnss_filtering")
    )
    pipeline = cfg.get("pipeline", {})
    links = _resolve(pipeline.get("links_file"), output_dir / "precomputed_links.csv")
    delays = _resolve(pipeline.get("delays_file"), output_dir / "precomputed_delays.csv")
    return links, delays, cfg.get("plasma", {})


def append_checkpoint(path: Path, row: dict) -> None:
    pd.DataFrame([row]).to_csv(path, mode="a", header=not path.exists(), index=False)


def read_existing(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path)


def link_identity(row: dict | pd.Series) -> tuple:
    """Return a stable identity for a link row.

    Reusing by link_id alone is unsafe because a regenerated link table can keep
    the same numeric IDs while changing PRN, epoch, or geometry.  Millimeter
    rounding is tight enough for cache validation while avoiding CSV formatting
    noise from pandas/C++ streams.
    """

    ints = tuple(int(row[col]) for col in IDENTITY_INT_COLUMNS)
    floats = tuple(round(float(row[col]), 3) for col in IDENTITY_FLOAT_COLUMNS)
    return ints + floats


def current_link_identities(links: pd.DataFrame) -> dict[int, tuple]:
    return {int(row["link_id"]): link_identity(row) for _, row in links.iterrows()}


def matching_existing_link_ids(existing: pd.DataFrame, current: dict[int, tuple]) -> set[int]:
    if existing.empty or "link_id" not in existing:
        return set()
    missing = [col for col in IDENTITY_INT_COLUMNS + IDENTITY_FLOAT_COLUMNS if col not in existing]
    if missing:
        return set()
    done: set[int] = set()
    for _, row in existing.iterrows():
        link_id = int(row["link_id"])
        if link_id in current and link_identity(row) == current[link_id]:
            done.add(link_id)
    return done


def filter_current_links(result: pd.DataFrame, current: dict[int, tuple]) -> pd.DataFrame:
    if result.empty or "link_id" not in result:
        return result
    missing = [col for col in IDENTITY_INT_COLUMNS + IDENTITY_FLOAT_COLUMNS if col not in result]
    if missing:
        return pd.DataFrame()
    keep = []
    for _, row in result.iterrows():
        link_id = int(row["link_id"])
        keep.append(link_id in current and link_identity(row) == current[link_id])
    return result.loc[keep].copy()


def raytrace_one(args):
    _, row, plasma = args
    import pylupnt as pnt

    pnt.set_iri_model(str(plasma.get("raytrace", {}).get("iri_model", "IRI2007")))

    config = pnt.RayTraceConfig()
    config.freq_Hz = FREQ_HZ.get(str(row["frequency"]), FREQ_HZ["L1"])
    config.step_size = float(plasma.get("raytrace", {}).get("step_size_km", 100.0))
    config.correction = bool(plasma.get("raytrace", {}).get("correction", False))
    config.fine_correction = bool(plasma.get("raytrace", {}).get("fine_correction", False))
    config.cutoff_r = (
        float(plasma.get("raytrace", {}).get("cutoff_radius_re", 4.0)) * EARTH_RADIUS_KM
    )
    config.gradn_dx = float(plasma.get("raytrace", {}).get("gradient_step_km", 1.0))
    config.integ_method = str(plasma.get("raytrace", {}).get("integrator", "RK4"))
    config.correction_method = str(
        plasma.get("raytrace", {}).get("correction_method", "neldermead")
    )
    config.kp = float(plasma.get("raytrace", {}).get("kp", 3.0))
    config.use_fortran_gcpm = bool(plasma.get("raytrace", {}).get("use_fortran_gcpm", True))
    config.corr_tol = float(plasma.get("raytrace", {}).get("correction_tolerance_m", 100.0))
    config.compute_higher_order = bool(plasma.get("raytrace", {}).get("compute_higher_order", True))
    config.use_adaptive_step = bool(plasma.get("raytrace", {}).get("use_adaptive_step", True))
    config.straight_ray = bool(plasma.get("raytrace", {}).get("straight_ray", True))

    tx_m = np.array([row["tx_x_ecef_m"], row["tx_y_ecef_m"], row["tx_z_ecef_m"]], dtype=float)
    rx_m = np.array([row["rx_x_ecef_m"], row["rx_y_ecef_m"], row["rx_z_ecef_m"]], dtype=float)
    tangent_radius_m, tangent_altitude_m = tangent_radius_altitude_m(tx_m, rx_m)
    tx = tx_m / 1000.0
    rx = rx_m / 1000.0
    profile = pnt.trace_ray(float(row["t_utc"]), tx, rx, config, False, False)

    iono_delay = float(profile.tec_delay_m)
    if config.compute_higher_order:
        iono_delay += float(profile.second_delay_m) + float(profile.third_delay_m)
    final_pos_err_m = float(np.linalg.norm(np.asarray(profile.corr_final_pos_err, dtype=float)))

    out = dict(row)
    out.update(
        {
            "ionosphere_plasma_delay_m": iono_delay,
            "tangent_radius_earth_center_m": float(
                row.get("tangent_radius_earth_center_m", tangent_radius_m)
            ),
            "tangent_altitude_m": float(row.get("tangent_altitude_m", tangent_altitude_m)),
            "tecu": float(profile.tecu),
            "tec_delay_m": float(profile.tec_delay_m),
            "second_delay_m": float(profile.second_delay_m),
            "third_delay_m": float(profile.third_delay_m),
            "dist_bend_m": float(profile.dist_bend_m),
            "tec_delay_bend_m": float(profile.tec_delay_bend_m),
            "max_sep_line_m": float(profile.max_sep_line_m),
            "final_pos_err_m": final_pos_err_m,
            "total_delay_m": float(profile.total_delay_m),
        }
    )
    return out


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Precompute GNSS ionosphere/plasma delays with GCPM."
    )
    parser.add_argument(
        "config", nargs="?", default="projects/GNSS_Filtering/gnss_filtering_config.yaml"
    )
    parser.add_argument("--links", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--workers", type=int, default=max(1, mp.cpu_count() - 1))
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--serial", action="store_true")
    args = parser.parse_args()

    config_path = Path(args.config)
    links_path, output_path, plasma = load_pipeline_config(config_path)
    if args.links:
        links_path = args.links
    if args.output:
        output_path = args.output

    links = pd.read_csv(links_path, comment="#")
    if "tangent_altitude_m" in links:
        occulted = links["tangent_altitude_m"].astype(float) < 0.0
        if occulted.any():
            print(
                f"Skipping {int(occulted.sum())} Earth-occulted links with negative tangent altitude"
            )
            links = links.loc[~occulted].copy()
    current_identities = current_link_identities(links)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = output_path.with_suffix(output_path.suffix + ".tmp")

    if args.overwrite and temp_path.exists():
        temp_path.unlink()

    jobs = [(idx, row.to_dict(), plasma) for idx, row in links.iterrows()]
    if not args.overwrite:
        done = set()
        for existing_path in (output_path, temp_path):
            existing = read_existing(existing_path)
            done.update(matching_existing_link_ids(existing, current_identities))
        jobs = [job for job in jobs if int(job[1]["link_id"]) not in done]

    print(f"Loaded {len(links)} links from {links_path}")
    print(
        f"Computing {len(jobs)} missing delays with {1 if args.serial else args.workers} worker(s)"
    )
    print(f"Checkpointing completed rays to {temp_path}")

    rows = []
    if args.serial or args.workers == 1:
        for job in jobs:
            row = raytrace_one(job)
            rows.append(row)
            append_checkpoint(temp_path, row)
    else:
        try:
            mp.set_start_method("spawn")
        except RuntimeError:
            pass
        with mp.Pool(processes=args.workers) as pool:
            for i, row in enumerate(pool.imap_unordered(raytrace_one, jobs), start=1):
                rows.append(row)
                append_checkpoint(temp_path, row)
                if i % 25 == 0 or i == len(jobs):
                    print(f"  completed {i}/{len(jobs)}", flush=True)

    frames = []
    if output_path.exists() and not args.overwrite:
        frames.append(read_existing(output_path))
    if temp_path.exists():
        frames.append(read_existing(temp_path))
    result = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(rows)
    result = filter_current_links(result, current_identities)
    if result.empty:
        result = links.copy()
    if "link_id" in result:
        result = result.drop_duplicates("link_id", keep="last").sort_values("link_id")
    result.to_csv(output_path, index=False)
    print(f"Wrote delays to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

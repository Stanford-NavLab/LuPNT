"""Stage 1 + Stage 2 precomputation for the ex6 lunar GNSS ODTS example.

These two stages are the slow part of the example and are painful to run inside a notebook.
Run this once from a terminal instead::

    pixi run python python/examples/ex6_precompute.py            # generate the cache
    pixi run python python/examples/ex6_precompute.py --force    # wipe + regenerate
    pixi run python python/examples/ex6_precompute.py --workers 8

It writes into ``output/python_examples/ex6_gnss_odts_data/`` (under the repo root):

  * Stage 1 -- ``precomputed_links.csv``: GNSS link geometry / C-N0 (``sim.precompute()``).
  * Stage 2 -- ``precomputed_delays.csv``: GCPM/IRI2007 plasmaspheric delay, ray-traced on a
    coarse cadence and interpolated to 1 Hz.

The notebook builds the same baseline config (``ex6_gnss_odts_config.build_config``), so after
this script runs its Stage 1 cell is an instant cache hit and its Stage 2 cell just loads the CSV.
The C++ link cache is intentionally keyed only on measurement-generation inputs; filter tuning
and pseudorange/Doppler/TDCP/ionosphere-free selection can change without regenerating the
precomputed measurements. ``compute_delays`` is also importable, so the notebook can fall back to
generating the delays itself if the cache is missing.
"""

import argparse
import hashlib
import json
import multiprocessing as mp
import os
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd

# Importing the config module runs the pylupnt import guard (correct sys.path), so importing
# pnt afterwards resolves to this repo's build.
from ex6_gnss_odts_config import (
    DATA_DIR,
    OUTPUT_DIR,
    PLASMA_DELAY_DT_S,
    PLASMA_DELAY_PROGRESS_RAYS,
    build_config,
    describe,
)
import pylupnt as pnt  # noqa: E402

FREQ_HZ = {
    "L1": 1575.42e6,
    "L2": 1227.60e6,
    "L5": 1176.45e6,
    "E1": 1575.42e6,
    "E5a": 1176.45e6,
}
PRIMARY_FREQUENCIES = {"L1", "E1"}
EARTH_RADIUS_KM = 6371.2
DELAY_COLS = ["tecu", "tec_delay_m", "second_delay_m", "third_delay_m", "ionosphere_plasma_delay_m"]
KEY_COLS = ["gnss_const_id", "prn", "frequency_id"]  # identifies one satellite-frequency link

_RT_CONFIG = None  # per-worker RayTraceConfig, built in _init_worker


# ---------------------------------------------------------------------------------------------
# Cache tagging: the delays are reused only when the simulation settings they depend on match.
# ---------------------------------------------------------------------------------------------
def _links_fingerprint(cfg):
    """The C++ Stage-1 link-cache fingerprint, written by precompute() next to the links CSV."""
    meta = Path(str(cfg.links_file) + ".meta")
    return meta.read_text().strip() if meta.exists() else ""


def delays_fingerprint(cfg, plasma_dt_s=PLASMA_DELAY_DT_S):
    """A tag capturing everything the Stage-2 delays depend on: the Stage-1 link geometry (via
    its fingerprint) plus the plasma ray-trace settings and interpolation cadence. The cached
    delays are reused only when this fingerprint matches the current config."""
    p = cfg.plasma
    tag = {
        "links": _links_fingerprint(cfg),
        "plasma_dt_s": plasma_dt_s,
        "iri_model": "IRI2007",
        "simulate_truth": p.simulate_truth,
        "kp": p.raytrace_kp,
        "rz12": p.raytrace_rz12,
        "step_size_km": p.raytrace_step_size_km,
        "correction": p.raytrace_correction,
        "fine_correction": p.raytrace_fine_correction,
        "cutoff_radius_re": p.raytrace_cutoff_radius_re,
        "gradient_step_km": p.raytrace_gradient_step_km,
        "integrator": p.raytrace_integrator,
        "correction_method": p.raytrace_correction_method,
        "correction_tolerance_m": p.raytrace_correction_tolerance_m,
        "compute_higher_order": p.raytrace_compute_higher_order,
        "use_adaptive_step": p.raytrace_use_adaptive_step,
        "straight_ray": p.raytrace_straight_ray,
        "use_fortran_gcpm": p.raytrace_use_fortran_gcpm,
    }
    return hashlib.md5(json.dumps(tag, sort_keys=True).encode()).hexdigest()


def _delays_meta_path(cfg):
    return Path(str(cfg.delays_file) + ".meta")


def delays_cache_valid(cfg, plasma_dt_s=PLASMA_DELAY_DT_S):
    """True iff a cached delays CSV exists AND was generated with the current settings."""
    if not (Path(cfg.delays_file).exists() and _delays_meta_path(cfg).exists()):
        return False
    return _delays_meta_path(cfg).read_text().strip() == delays_fingerprint(cfg, plasma_dt_s)


def _rt_params(cfg):
    """Picklable dict of RayTraceConfig attribute values, sent to each worker."""
    p = cfg.plasma
    return dict(
        step_size=p.raytrace_step_size_km,
        correction=p.raytrace_correction,
        fine_correction=p.raytrace_fine_correction,
        cutoff_r=p.raytrace_cutoff_radius_re * EARTH_RADIUS_KM,
        gradn_dx=p.raytrace_gradient_step_km,
        integ_method=p.raytrace_integrator,
        correction_method=p.raytrace_correction_method,
        kp=p.raytrace_kp,
        rz12=p.raytrace_rz12,
        use_fortran_gcpm=p.raytrace_use_fortran_gcpm,
        corr_tol=p.raytrace_correction_tolerance_m,
        compute_higher_order=p.raytrace_compute_higher_order,
        use_adaptive_step=p.raytrace_use_adaptive_step,
        straight_ray=p.raytrace_straight_ray,
    )


def _build_rt_config(params):
    rt = pnt.RayTraceConfig()
    for key, value in params.items():
        setattr(rt, key, value)
    return rt


def _init_worker(params):
    """Runs in each worker (fork or spawn): select the IRI model and build its RayTraceConfig."""
    global _RT_CONFIG
    pnt.set_iri_model("IRI2007")
    _RT_CONFIG = _build_rt_config(params)


def _make_task(row):
    """Serializable ray-trace inputs for one link row."""
    return (
        float(row["t_utc"]),
        np.array([row["tx_x_ecef_m"], row["tx_y_ecef_m"], row["tx_z_ecef_m"]]) / 1e3,
        np.array([row["rx_x_ecef_m"], row["rx_y_ecef_m"], row["rx_z_ecef_m"]]) / 1e3,
        FREQ_HZ.get(str(row["frequency"]), FREQ_HZ["L1"]),
    )


def _trace_worker(task):
    """Ray-trace one link through GCPM/IRI2007 -> dict of delay columns (uses _RT_CONFIG)."""
    t_utc, tx_km, rx_km, freq_hz = task
    _RT_CONFIG.freq_Hz = freq_hz
    profile = pnt.trace_ray(t_utc, tx_km, rx_km, _RT_CONFIG, False, False)
    tec_m = float(profile.tec_delay_m)
    higher_m = (
        float(profile.second_delay_m) + float(profile.third_delay_m)
        if _RT_CONFIG.compute_higher_order
        else 0.0
    )
    return {
        "tecu": float(profile.tecu),
        "tec_delay_m": tec_m,
        "second_delay_m": float(profile.second_delay_m),
        "third_delay_m": float(profile.third_delay_m),
        "ionosphere_plasma_delay_m": tec_m + higher_m,
    }


def _trace_links_parallel(pool, rows, label, progress_every=10):
    """Trace every row of `rows` on `pool`, returning {row_index: delay_dict}."""
    progress_every = max(1, int(progress_every))
    idxs, tasks = [], []
    for idx, row in rows.iterrows():
        idxs.append(idx)
        tasks.append(_make_task(row))
    records, t_start = {}, time.time()
    for k, (idx, rec) in enumerate(
        zip(idxs, pool.imap(_trace_worker, tasks, chunksize=1)), start=1
    ):
        records[idx] = rec
        if k % progress_every == 0 or k == len(tasks):
            print(f"  {k}/{len(tasks)} {label} rays traced ({time.time() - t_start:.0f}s elapsed)")
    return records


def _chunk_ranges(n_items, n_chunks):
    """Even half-open ranges covering [0, n_items)."""
    n_chunks = max(1, min(n_chunks, n_items))
    base, extra = divmod(n_items, n_chunks)
    ranges, start = [], 0
    for i in range(n_chunks):
        stop = start + base + (1 if i < extra else 0)
        if start < stop:
            ranges.append((start, stop))
        start = stop
    return ranges


def _link_range_worker(args):
    """Worker process entry point for Stage 1 link chunks."""
    epoch_begin, epoch_end, links_file, threads_per_worker = args
    worker_cfg = build_config()
    worker_cfg.links_file = links_file
    worker_cfg.precompute_num_threads = max(1, int(threads_per_worker))
    pnt.precompute_lunar_gnss_odts_links_range(worker_cfg, epoch_begin, epoch_end)
    return links_file


def _remove_file(path):
    path = Path(path)
    if path.exists():
        path.unlink()


def _invalidate_link_dependents(cfg):
    _remove_file(cfg.links_file)
    _remove_file(str(cfg.links_file) + ".meta")
    _remove_file(cfg.delays_file)
    _remove_file(str(cfg.delays_file) + ".meta")


def _merge_link_chunks(chunk_files, output_file):
    frames = [pd.read_csv(path) for path in chunk_files]
    links = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if not links.empty:
        links["link_id"] = np.arange(len(links), dtype=int)
    links.to_csv(output_file, index=False)
    return links


def precompute_links(cfg, n_workers=None, threads_per_worker=1):
    """Stage 1: build reusable link geometry/measurement metadata. Returns (sim, links_df)."""
    if pnt.lunar_gnss_odts_link_cache_valid(cfg):
        print(
            f"Stage 1: cached links match the current measurement-generation settings "
            f"({cfg.links_file}); skipping."
        )
        return None, pd.read_csv(cfg.links_file)

    n_epochs = pnt.lunar_gnss_odts_precompute_epoch_count(cfg)
    n_workers = n_workers or (os.cpu_count() or 1)
    n_workers = max(1, min(int(n_workers), n_epochs))

    t0 = time.time()
    if n_workers == 1:
        # Full-range link precompute + cache-metadata write (the former
        # LunarGnssODTSSimulation.precompute()), via the retained free functions.
        _invalidate_link_dependents(cfg)
        pnt.precompute_lunar_gnss_odts_links_range(cfg, 0, n_epochs)
        pnt.finalize_lunar_gnss_odts_link_cache(cfg)
        print(f"Stage 1: precompute() took {time.time() - t0:.1f} s")
    else:
        ranges = _chunk_ranges(n_epochs, n_workers)
        chunk_dir = OUTPUT_DIR / "_link_chunks"
        if chunk_dir.exists():
            shutil.rmtree(chunk_dir)
        chunk_dir.mkdir(parents=True, exist_ok=True)
        _invalidate_link_dependents(cfg)

        tasks = [
            (start, stop, str(chunk_dir / f"precomputed_links_{i:04d}.csv"), threads_per_worker)
            for i, (start, stop) in enumerate(ranges)
        ]
        print(
            f"Stage 1: precomputing {n_epochs} epochs in {len(tasks)} Python worker "
            f"processes ({threads_per_worker} C++ thread(s)/worker)..."
        )
        ctx = mp.get_context("spawn")
        with ctx.Pool(len(tasks)) as pool:
            chunk_files = pool.map(_link_range_worker, tasks)
        links = _merge_link_chunks(chunk_files, cfg.links_file)
        pnt.finalize_lunar_gnss_odts_link_cache(cfg)
        shutil.rmtree(chunk_dir)
        print(f"Stage 1: parallel precompute took {time.time() - t0:.1f} s")
        print(
            f"  {len(links)} links across {links['epoch_index'].nunique()} epochs "
            f"({links.groupby('epoch_index').size().mean():.1f} links/epoch)"
        )
        return None, links

    links = pd.read_csv(cfg.links_file)
    print(
        f"  {len(links)} links across {links['epoch_index'].nunique()} epochs "
        f"({links.groupby('epoch_index').size().mean():.1f} links/epoch)"
    )
    return sim, links


def compute_delays(
    cfg,
    links,
    plasma_dt_s=PLASMA_DELAY_DT_S,
    n_workers=None,
    progress_every=PLASMA_DELAY_PROGRESS_RAYS,
):
    """Stage 2: ray-trace primary-frequency delays, interpolate to 1 Hz, write CSV.

    Only non-occulted primary-frequency links (GPS L1, Galileo E1) are traced. The C++
    ionosphere-free pseudorange derives secondary-frequency delay from the matching primary
    delay, and TDCP also uses the primary-frequency delay. Returns the full per-link delays
    DataFrame.
    """
    pnt.set_iri_model("IRI2007")
    params = _rt_params(cfg)

    occulted = links["tangent_altitude_m"] < 0.0
    traceable = (~occulted) & links["frequency"].isin(PRIMARY_FREQUENCIES)
    if occulted.any():
        print(f"Skipping {int(occulted.sum())} Earth-occulted links")

    # Anchor epochs: first epoch of each plasma_dt_s window (plus the final epoch), so every
    # 1 Hz link is bracketed or edge-clamped by traced samples of the same satellite link.
    epoch_time = links.groupby("epoch_index")["t_utc"].first()
    window = np.floor((epoch_time - epoch_time.min()) / plasma_dt_s).astype(int)
    anchor_epochs = set(epoch_time.groupby(window).idxmin().tolist())
    anchor_epochs.add(int(epoch_time.index.max()))
    anchor_links = links.loc[links["epoch_index"].isin(anchor_epochs) & traceable]

    n_workers = n_workers or (os.cpu_count() or 1)
    n_workers = max(1, min(n_workers, len(anchor_links)))
    print(
        f"Ray-tracing {len(anchor_links)} anchor links across {len(anchor_epochs)} epochs "
        f"({plasma_dt_s:.0f}s cadence) on {n_workers} processes; interpolating the rest to "
        f"{cfg.receiver_app.rate_hz:.0f} Hz..."
    )
    ctx = mp.get_context("fork")

    t0 = time.time()
    with ctx.Pool(n_workers, initializer=_init_worker, initargs=(params,)) as pool:
        anchor_records = _trace_links_parallel(pool, anchor_links, "anchor", progress_every)

        anchors = pd.DataFrame.from_dict(anchor_records, orient="index")
        anchors[KEY_COLS + ["t_utc"]] = anchor_links.loc[anchors.index, KEY_COLS + ["t_utc"]]

        delays = links.copy()
        for col in DELAY_COLS:
            delays[col] = 0.0

        anchors_by_key = {key: g.sort_values("t_utc") for key, g in anchors.groupby(KEY_COLS)}
        fallback_idx = []
        for key, grp in delays.loc[traceable].groupby(KEY_COLS):
            a = anchors_by_key.get(key)
            if a is None:
                fallback_idx.extend(grp.index.tolist())  # short pass never sampled at an anchor
                continue
            for col in DELAY_COLS:
                delays.loc[grp.index, col] = np.interp(
                    grp["t_utc"].to_numpy(), a["t_utc"].to_numpy(), a[col].to_numpy()
                )

        if fallback_idx:
            print(
                f"Ray-tracing {len(fallback_idx)} links for satellites never seen at an "
                f"anchor epoch (passes shorter than {plasma_dt_s:.0f}s)..."
            )
            fb = _trace_links_parallel(pool, links.loc[fallback_idx], "fallback", progress_every)
            for idx, rec in fb.items():
                for col, val in rec.items():
                    delays.loc[idx, col] = val

    n_traced = len(anchor_links) + len(fallback_idx)
    n_links = int(traceable.sum())
    print(
        f"Ray-traced {n_traced}/{n_links} non-occulted primary-frequency links "
        f"({100.0 * n_traced / max(n_links, 1):.1f}%); interpolated the remaining "
        f"{n_links - n_traced}."
    )
    delays_path = Path(cfg.delays_file)
    delays.to_csv(delays_path, index=False)
    # Tag the cache with the settings it was generated from so it is reused only on a match.
    _delays_meta_path(cfg).write_text(delays_fingerprint(cfg, plasma_dt_s) + "\n")
    print(f"Wrote {delays_path} in {time.time() - t0:.1f} s")
    return delays


def main():
    ap = argparse.ArgumentParser(description="Precompute ex6 GNSS ODTS Stage 1 + Stage 2.")
    ap.add_argument(
        "--force",
        action="store_true",
        help="delete the cached output directory and regenerate everything",
    )
    ap.add_argument(
        "--workers",
        type=int,
        default=None,
        help="default worker process count for Stage 1 links and Stage 2 delays",
    )
    ap.add_argument(
        "--link-workers",
        type=int,
        default=None,
        help="Stage 1 link worker processes (default: --workers or os.cpu_count())",
    )
    ap.add_argument(
        "--delay-workers",
        type=int,
        default=None,
        help="Stage 2 ray-trace worker processes (default: --workers or os.cpu_count())",
    )
    ap.add_argument(
        "--link-threads-per-worker",
        type=int,
        default=1,
        help="C++ OpenMP constellation threads inside each Stage 1 worker "
        "(default: 1 to avoid oversubscribing cores)",
    )
    ap.add_argument(
        "--delay-progress-rays",
        type=int,
        default=None,
        help=f"print Stage 2 ray-trace progress every N completed rays "
        f"(default from ex6_gnss_odts_config.py: "
        f"{PLASMA_DELAY_PROGRESS_RAYS})",
    )
    args = ap.parse_args()

    cfg = build_config()
    if args.force and OUTPUT_DIR.exists():
        print(f"--force: removing {OUTPUT_DIR}")
        shutil.rmtree(OUTPUT_DIR)
        OUTPUT_DIR.mkdir(exist_ok=True)
    print(describe(cfg))
    print("-" * 70)

    link_workers = args.link_workers if args.link_workers is not None else args.workers
    delay_workers = args.delay_workers if args.delay_workers is not None else args.workers
    delay_progress_rays = (
        args.delay_progress_rays
        if args.delay_progress_rays is not None
        else PLASMA_DELAY_PROGRESS_RAYS
    )
    _, links = precompute_links(cfg, link_workers, args.link_threads_per_worker)

    delays_path = Path(cfg.delays_file)
    if delays_cache_valid(cfg, PLASMA_DELAY_DT_S):
        print(
            f"Stage 2: cached delays match the current settings ({delays_path}); skipping "
            f"(use --force to redo)."
        )
    else:
        if delays_path.exists():
            print("Stage 2: cached delays are for different settings -- regenerating.")
        print("Stage 2: ray-tracing plasmaspheric delays...")
        compute_delays(cfg, links, PLASMA_DELAY_DT_S, delay_workers, delay_progress_rays)

    print("-" * 70)
    print("Done. Open ex6_gnss_odts.ipynb and run it -- Stages 1 & 2 are now cache hits.")


if __name__ == "__main__":
    main()

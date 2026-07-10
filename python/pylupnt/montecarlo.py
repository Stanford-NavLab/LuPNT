"""Run Monte-Carlo trials of a ``pnt.Simulation`` config in parallel.

A Monte-Carlo trial is one independent ``pnt.Simulation`` run, built from a deep
copy of the scenario config with a distinct top-level ``seed`` (which the engine
propagates to every agent/app) and a distinct per-run ``output_dir`` so file
outputs never collide. Monte-Carlo settings therefore live at the *simulation*
level, not on an agent or app: the config carries a base ``seed``, and the number
of trials + worker count are arguments here.

The C++ engine stays serial; parallelism is across **processes** (``multiprocessing``),
so the GIL and the thread-local C++ RNG are not shared between trials. Each worker
rebuilds the simulation from the (picklable) config dict — ``pnt.Simulation`` objects
themselves are not picklable — runs it, and returns whatever the picklable
``extract`` callable pulls out of the finished sim.

Example
-------
>>> import pylupnt as pnt
>>> def final_error(sim):                       # must be a top-level (picklable) function
...     app = sim.get_agent("gs_manager").get_application()
...     return float(app.est_central()[-1, 0])
>>> summaries = pnt.run_monte_carlo(config, n_runs=16, extract=final_error, n_workers=8)
"""

from __future__ import annotations

import copy
import multiprocessing as mp
import os
from typing import Any, Callable, List, Optional

__all__ = ["run_monte_carlo"]


def _set_output_dirs(node: Any, run_dir: str) -> None:
    """Point every nested ``output_dir`` key at this run's private directory."""
    if isinstance(node, dict):
        for k, v in node.items():
            if k == "output_dir":
                node[k] = run_dir
            else:
                _set_output_dirs(v, run_dir)
    elif isinstance(node, (list, tuple)):
        for v in node:
            _set_output_dirs(v, run_dir)


def _run_trial(payload):
    config, seed, run_index, output_root, extract = payload
    import pylupnt as pnt  # imported inside the worker process

    cfg = copy.deepcopy(config)
    cfg["seed"] = seed  # the engine propagates this to all agents/apps
    if output_root is not None:
        run_dir = os.path.join(output_root, f"run_{run_index:03d}")
        os.makedirs(run_dir, exist_ok=True)
        _set_output_dirs(cfg, run_dir)

    sim = pnt.Simulation(cfg)
    sim.run()
    return extract(sim) if extract is not None else None


def run_monte_carlo(
    config: dict,
    n_runs: int,
    extract: Optional[Callable[[Any], Any]] = None,
    n_workers: Optional[int] = None,
    seed0: int = 0,
    output_root: Optional[str] = None,
) -> List[Any]:
    """Run ``n_runs`` Monte-Carlo trials of ``config`` and return the per-run results.

    Parameters
    ----------
    config : dict
        The scenario config (as passed to ``pnt.Simulation``). Deep-copied per run.
    n_runs : int
        Number of independent trials; trial ``i`` uses ``seed = seed0 + i``.
    extract : callable, optional
        ``extract(sim) -> picklable`` pulled from each finished ``pnt.Simulation``.
        Must be a top-level function (multiprocessing pickles it). If ``None``,
        results are ``None`` (run for side effects / file outputs only).
    n_workers : int, optional
        Process-pool size (default: ``min(n_runs, cpu_count())``). ``1`` runs serially
        in-process (handy for debugging).
    seed0 : int
        Base seed for trial 0.
    output_root : str, optional
        If given, each trial's file ``output_dir`` keys are redirected to
        ``<output_root>/run_<i>`` so parallel trials don't overwrite each other.

    Returns
    -------
    list
        ``[extract(sim_0), extract(sim_1), ...]`` in trial order.
    """
    if n_runs < 1:
        raise ValueError("n_runs must be >= 1")
    payloads = [(config, seed0 + i, i, output_root, extract) for i in range(n_runs)]

    if n_workers is None:
        n_workers = min(n_runs, os.cpu_count() or 1)
    if n_workers <= 1:
        return [_run_trial(p) for p in payloads]

    # `spawn` gives each worker a clean interpreter (mirrors ex6_precompute.py) so the
    # thread-local C++ RandomEngine and any global state start fresh and independent.
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=n_workers) as pool:
        return pool.map(_run_trial, payloads)

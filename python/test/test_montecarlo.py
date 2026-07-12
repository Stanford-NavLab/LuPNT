"""Tests for pylupnt.run_monte_carlo (the simulation-level Monte-Carlo helper)."""

from pathlib import Path

import yaml
import pytest

import pylupnt as pnt


def _extract_ok(sim):
    """Top-level (picklable) extractor: returns a trivial per-run marker."""
    return sim is not None


def _ephemeris_config():
    for base in [Path.cwd(), *Path.cwd().parents]:
        p = base / "configs" / "ephemeris.yaml"
        if p.exists():
            return yaml.safe_load(open(p))
    return None


def test_run_monte_carlo_rejects_nonpositive_runs():
    with pytest.raises(ValueError):
        pnt.run_monte_carlo({}, n_runs=0)


def test_run_monte_carlo_serial_runs_each_trial(tmp_path):
    cfg = _ephemeris_config()
    if cfg is None:
        pytest.skip("configs/ephemeris.yaml not found")
    # n_workers=1 runs each trial in-process (no multiprocessing), fast for the one-shot
    # ephemeris study, and per-run output dirs keep the trials isolated.
    results = pnt.run_monte_carlo(
        cfg, n_runs=2, extract=_extract_ok, n_workers=1, seed0=7, output_root=str(tmp_path)
    )
    assert results == [True, True]
    assert (tmp_path / "run_000").exists()
    assert (tmp_path / "run_001").exists()

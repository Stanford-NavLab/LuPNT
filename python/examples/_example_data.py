"""Seed slow-to-generate example outputs from the shipped precomputed data.

Examples ex6/ex7/ex9/ex10/ex11 each need a precomputed ``output/python_examples/exN_*``
directory that's normally produced by running the matching ``exN_run_*.py``/``exN_precompute.py``
script (minutes, sometimes much more, for ex6's plasma ray-trace). A copy of each of these
directories ships in ``data/LuPNT_data/examples/`` (part of the LuPNT_data mirror — see
cmake/FetchLuPNTData.cmake) so a fresh checkout doesn't have to pay that cost just to view the
tutorial notebooks.

``seed_example_output`` copies the shipped copy into place the first time it's needed, so every
existing "does output/python_examples/exN_* already exist" check downstream (notebook cells, the
exN_run_*.py scripts, ex6's C++ link-precompute cache fingerprint) finds it exactly as if it had
been freshly computed, and skips recomputation. Delete the seeded output/python_examples/exN_*
directory to force a real recompute even when shipped data exists.
"""

from __future__ import annotations

import shutil
from pathlib import Path


def _repo_root() -> Path:
    for base in [Path.cwd(), *Path.cwd().parents]:
        if (base / "python" / "pylupnt" / "__init__.py").exists():
            return base
    return Path.cwd()


def seed_example_output(name: str, markers: tuple[str, ...] = ("meta.json",)) -> Path:
    """Return ``output/python_examples/<name>``, copying it from
    ``data/LuPNT_data/examples/<name>`` first if it's missing any of ``markers`` and the
    shipped copy has them."""
    repo = _repo_root()
    out_dir = repo / "output" / "python_examples" / name
    if all((out_dir / marker).exists() for marker in markers):
        return out_dir
    shipped_dir = repo / "data" / "LuPNT_data" / "examples" / name
    if all((shipped_dir / marker).exists() for marker in markers):
        out_dir.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(shipped_dir, out_dir, dirs_exist_ok=True)
        print(f"[example-data] seeded {out_dir} from shipped {shipped_dir}")
    return out_dir

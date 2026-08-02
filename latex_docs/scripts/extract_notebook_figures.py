#!/usr/bin/env python3
"""Extract embedded figures from the LuPNT example notebooks for the LaTeX manual.

Purpose
-------
The notebooks under ``python/examples/*.ipynb`` are committed **with their
outputs** (``nbstripout`` is configured in this repository to preserve them), so
every plot a tutorial produces is already sitting in the notebook JSON as a
base64 ``image/png`` output.  The LaTeX manual in ``latex_docs/`` is a text-only
tree -- its ``.gitignore`` rejects images -- so the tutorials chapter cannot
reference those plots directly.  This script bridges the two: it decodes the
embedded PNGs and writes them into the generated-artefact directory that
``preamble.tex`` exposes through the ``\\datafig`` macro.

Output layout
-------------
Each image is written to::

    $LUPNT_DATA_PATH/docs/tutorials/<example>_<n>.png

where ``<example>`` is the notebook stem with the trailing descriptive part
kept (e.g. ``ex1_propagate_orbit``) and ``<n>`` is a 1-based counter over the
images produced by that notebook, in execution order (cell order, then output
order within a cell).  A chapter then includes one with::

    \\datafig{tutorials/ex1_propagate_orbit_3.png}

Only ``image/png`` outputs are extracted.  Cells whose output is a Cesium/HTML
widget, a Plotly bundle, or plain text produce nothing.

Usage
-----
::

    # writes to $LUPNT_DATA_PATH/docs/tutorials/
    LUPNT_DATA_PATH=/path/to/LuPNT_data python latex_docs/scripts/extract_notebook_figures.py

    # explicit destination and a subset of notebooks
    python latex_docs/scripts/extract_notebook_figures.py \\
        --notebook-dir python/examples --out /tmp/figs --only ex1 --only ex15

    # list what would be written without touching the filesystem
    python latex_docs/scripts/extract_notebook_figures.py --dry-run

Options
-------
``--notebook-dir``  directory holding the ``ex*.ipynb`` files
                    (default: ``python/examples`` relative to the repo root).
``--out``           destination directory.  Defaults to
                    ``$LUPNT_DATA_PATH/docs/tutorials``.
``--only``          repeatable notebook-stem prefix filter (``--only ex2``).
``--min-bytes``     skip decoded images smaller than this (default 2000), which
                    drops colour-bar strips and 1-pixel spacers.
``--dry-run``       report what would be written and exit.

The script is idempotent: re-running it overwrites the same file names, so the
figure referenced by the manual stays stable as long as the notebook's plot
order does not change.  It has no third-party dependencies.
"""

from __future__ import annotations

import argparse
import base64
import binascii
import json
import os
import struct
import sys
from pathlib import Path


def png_size(data: bytes) -> tuple[int, int] | None:
    """Return (width, height) of a PNG buffer, or None if it is not a PNG."""
    if len(data) < 24 or data[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    width, height = struct.unpack(">II", data[16:24])
    return int(width), int(height)


def iter_png_outputs(nb: dict):
    """Yield decoded PNG payloads from a notebook, in execution order."""
    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        for out in cell.get("outputs", []):
            bundle = out.get("data") or {}
            payload = bundle.get("image/png")
            if payload is None:
                continue
            if isinstance(payload, list):
                payload = "".join(payload)
            try:
                yield base64.b64decode(payload)
            except (binascii.Error, ValueError) as exc:  # pragma: no cover
                print(f"    ! undecodable image/png output: {exc}", file=sys.stderr)


def default_out_dir() -> Path:
    data_path = os.environ.get("LUPNT_DATA_PATH")
    if not data_path:
        sys.exit(
            "LUPNT_DATA_PATH is not set and --out was not given.\n"
            "Set it to the LuPNT_data directory, e.g.\n"
            "  LUPNT_DATA_PATH=/path/to/LuPNT_data python "
            "latex_docs/scripts/extract_notebook_figures.py"
        )
    return Path(data_path) / "docs" / "tutorials"


def main(argv: list[str] | None = None) -> int:
    repo_root = Path(__file__).resolve().parents[2]

    parser = argparse.ArgumentParser(
        description="Extract embedded notebook figures for the LaTeX manual.",
    )
    parser.add_argument(
        "--notebook-dir",
        type=Path,
        default=repo_root / "python" / "examples",
        help="directory holding the ex*.ipynb notebooks",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="destination directory (default: $LUPNT_DATA_PATH/docs/tutorials)",
    )
    parser.add_argument(
        "--only",
        action="append",
        default=[],
        metavar="PREFIX",
        help="only process notebooks whose stem starts with PREFIX (repeatable)",
    )
    parser.add_argument(
        "--min-bytes",
        type=int,
        default=2000,
        help="skip decoded images smaller than this many bytes",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="report what would be written, write nothing",
    )
    args = parser.parse_args(argv)

    out_dir = args.out if args.out is not None else default_out_dir()

    notebooks = sorted(args.notebook_dir.glob("*.ipynb"))
    if args.only:
        notebooks = [nb for nb in notebooks if any(nb.stem.startswith(p) for p in args.only)]
    if not notebooks:
        sys.exit(f"no notebooks matched under {args.notebook_dir}")

    if not args.dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)

    total = 0
    for nb_path in notebooks:
        try:
            nb = json.loads(nb_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            print(f"{nb_path.name}: unreadable ({exc})", file=sys.stderr)
            continue

        index = 0
        for data in iter_png_outputs(nb):
            if len(data) < args.min_bytes:
                continue
            index += 1
            name = f"{nb_path.stem}_{index}.png"
            dims = png_size(data)
            shape = f"{dims[0]}x{dims[1]}" if dims else "unknown"
            if args.dry_run:
                print(f"  would write {name}  ({shape}, {len(data)} B)")
            else:
                (out_dir / name).write_bytes(data)
                print(f"  {name}  ({shape}, {len(data)} B)")
            total += 1
        print(f"{nb_path.name}: {index} figure(s)")

    verb = "would extract" if args.dry_run else "extracted"
    print(f"\n{verb} {total} figure(s) into {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

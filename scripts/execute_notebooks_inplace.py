#!/usr/bin/env python
"""Execute every example notebook in ``python/examples`` and write outputs back in place.

Unlike ``run_example_notebooks.py`` (which runs on an in-memory copy for CI validation),
this populates each ``ex*.ipynb`` with fresh cell outputs so the rendered docs
(``nbsphinx_execute = "never"``) show real results. Run under ``pixi run`` so the
``lupnt`` kernel env (PYTHONPATH / LUPNT_DATA_PATH / ...) is present.

    pixi run python scripts/execute_notebooks_inplace.py                 # all
    pixi run python scripts/execute_notebooks_inplace.py --only ex1 ex4  # subset
    pixi run python scripts/execute_notebooks_inplace.py --skip ex13 ex14
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES_DIR = ROOT / "python" / "examples"


def natural_key(path: Path):
    m = re.search(r"ex(\d+)", path.stem)
    return (int(m.group(1)) if m else 1 << 30, path.stem)


def discover(only, skip):
    notebooks = sorted(EXAMPLES_DIR.glob("ex*.ipynb"), key=natural_key)
    if only:
        notebooks = [n for n in notebooks if any(pat in n.stem for pat in only)]
    if skip:
        notebooks = [n for n in notebooks if not any(pat in n.stem for pat in skip)]
    return notebooks


def run_one(path: Path, kernel: str, timeout: int, allow_errors: bool):
    nb = nbformat.read(path, as_version=4)
    client = NotebookClient(
        nb,
        timeout=timeout,
        kernel_name=kernel,
        resources={"metadata": {"path": str(path.parent)}},
        allow_errors=allow_errors,
    )
    start = time.perf_counter()
    err = None
    try:
        client.execute()
    except CellExecutionError as e:
        err = str(e).strip().splitlines()[-1]
    except Exception as e:
        err = f"{type(e).__name__}: {e}"
    # Write back whatever outputs we produced (partial on failure when allow_errors).
    nbformat.write(nb, path)
    return err is None, time.perf_counter() - start, err


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--only", nargs="+", metavar="SUBSTR")
    p.add_argument("--skip", nargs="+", metavar="SUBSTR")
    p.add_argument("--timeout", type=int, default=2400)
    p.add_argument("--kernel", default="lupnt")
    p.add_argument("--allow-errors", action="store_true", help="keep going within a notebook past a failing cell")
    args = p.parse_args()

    notebooks = discover(args.only, args.skip)
    if not notebooks:
        print(f"No notebooks matched under {EXAMPLES_DIR}", file=sys.stderr)
        return 1

    print(f"Executing {len(notebooks)} notebook(s) in place "
          f"[kernel={args.kernel}, timeout={args.timeout}s, allow_errors={args.allow_errors}]\n", flush=True)

    results = []
    for i, nb in enumerate(notebooks, 1):
        print(f"[{i}/{len(notebooks)}] {nb.name} ... ", end="", flush=True)
        ok, secs, err = run_one(nb, args.kernel, args.timeout, args.allow_errors)
        results.append((nb.name, ok, secs, err))
        print(f"{'PASS' if ok else 'FAIL'} ({secs:.1f}s)", flush=True)
        if err:
            print(f"        -> {err}", flush=True)

    passed = [r for r in results if r[1]]
    failed = [r for r in results if not r[1]]
    print("\n" + "=" * 60)
    print(f"Summary: {len(passed)} passed, {len(failed)} failed")
    for name, ok, secs, err in results:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}  ({secs:.1f}s)")
    if failed:
        print("\nFailures:")
        for name, ok, secs, err in failed:
            print(f"  - {name}: {err}")
    print("=" * 60, flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

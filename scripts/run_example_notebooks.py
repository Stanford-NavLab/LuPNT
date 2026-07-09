#!/usr/bin/env python
"""Execute every example notebook in ``python/examples`` and report a pass/fail summary.

This is the notebook analogue of ``pixi run test-py``: it runs each ``ex*.ipynb`` end to
end (in a fresh kernel) so we catch examples that have silently rotted after an API change.
Outputs are *not* written back to the source notebooks -- execution happens on an in-memory
copy -- so running this never produces a giant diff of re-rendered cell outputs.

All notebooks are forced onto the ``lupnt`` kernelspec (``pixi run install-kernel``), which
bakes in ``PYTHONPATH`` / ``LUPNT_DATA_PATH`` / ... regardless of the kernel recorded in each
notebook's metadata (some say ``lupnt``, some say ``python3``). Under ``pixi run`` the
activation env is already present too, so this works from either entry point.

Usage (typically ``pixi run run-notebooks``, which supplies the env + kernel):

    python scripts/run_example_notebooks.py                 # run all
    python scripts/run_example_notebooks.py --only ex1 ex4  # substring filter
    python scripts/run_example_notebooks.py --skip ex13     # exclude by substring
    python scripts/run_example_notebooks.py --timeout 600   # per-notebook seconds
    python scripts/run_example_notebooks.py --fail-fast     # stop at first failure

Exits non-zero if any executed notebook raised.
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

# scripts/ is one level below the repo root.
ROOT = Path(__file__).resolve().parents[1]
EXAMPLES_DIR = ROOT / "python" / "examples"


def natural_key(path: Path):
    """Sort ex1, ex2, ..., ex10 numerically rather than lexicographically."""
    m = re.search(r"ex(\d+)", path.stem)
    return (int(m.group(1)) if m else 1 << 30, path.stem)


def discover(only, skip):
    notebooks = sorted(EXAMPLES_DIR.glob("ex*.ipynb"), key=natural_key)
    if only:
        notebooks = [n for n in notebooks if any(pat in n.stem for pat in only)]
    if skip:
        notebooks = [n for n in notebooks if not any(pat in n.stem for pat in skip)]
    return notebooks


def run_one(path: Path, kernel: str, timeout: int):
    nb = nbformat.read(path, as_version=4)
    client = NotebookClient(
        nb,
        timeout=timeout,
        kernel_name=kernel,
        # Run from the notebook's own directory so relative paths resolve as they do
        # interactively; don't abort the run on the first bad cell -- we surface it below.
        resources={"metadata": {"path": str(path.parent)}},
        allow_errors=False,
    )
    start = time.perf_counter()
    try:
        client.execute()
        return True, time.perf_counter() - start, None
    except CellExecutionError as e:
        return False, time.perf_counter() - start, str(e).strip().splitlines()[-1]
    except Exception as e:  # kernel death, timeout, etc.
        return False, time.perf_counter() - start, f"{type(e).__name__}: {e}"


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--only",
        nargs="+",
        metavar="SUBSTR",
        help="run only notebooks whose name contains any of these",
    )
    p.add_argument(
        "--skip",
        nargs="+",
        metavar="SUBSTR",
        help="skip notebooks whose name contains any of these",
    )
    p.add_argument(
        "--timeout",
        type=int,
        default=1800,
        help="per-notebook cell timeout in seconds (default 1800)",
    )
    p.add_argument(
        "--kernel", default="lupnt", help="kernelspec name to execute with (default 'lupnt')"
    )
    p.add_argument("--fail-fast", action="store_true", help="stop at the first failing notebook")
    args = p.parse_args()

    notebooks = discover(args.only, args.skip)
    if not notebooks:
        print(f"No notebooks matched under {EXAMPLES_DIR}", file=sys.stderr)
        return 1

    print(
        f"Running {len(notebooks)} notebook(s) from {EXAMPLES_DIR.relative_to(ROOT)} "
        f"[kernel={args.kernel}, timeout={args.timeout}s]\n"
    )

    results = []
    for i, nb in enumerate(notebooks, 1):
        print(f"[{i}/{len(notebooks)}] {nb.name} ... ", end="", flush=True)
        ok, secs, err = run_one(nb, args.kernel, args.timeout)
        results.append((nb.name, ok, secs, err))
        print(f"{'PASS' if ok else 'FAIL'} ({secs:.1f}s)")
        if err:
            print(f"        -> {err}")
        if not ok and args.fail_fast:
            break

    passed = [r for r in results if r[1]]
    failed = [r for r in results if not r[1]]

    print("\n" + "=" * 60)
    print(
        f"Summary: {len(passed)} passed, {len(failed)} failed, "
        f"{len(notebooks) - len(results)} not run"
    )
    for name, ok, secs, err in results:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}  ({secs:.1f}s)")
    if failed:
        print("\nFailures:")
        for name, ok, secs, err in failed:
            print(f"  - {name}: {err}")
    print("=" * 60)

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

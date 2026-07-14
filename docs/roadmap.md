# LuPNT Roadmap & Release Notes

This document records the release/packaging state of LuPNT and the open threads a
future maintainer should know about. It is intentionally short and factual.

## Current status (v0.1.0)

- **Version:** `0.1.0` (declared in `cpp/lupnt/version.txt`, `pixi.toml`, and
  `python/pyproject.toml`). Tagged as a source release; **not** published to PyPI.
- **Install:** from source only, via `pixi` (see the README). There is no
  `pip install pylupnt`.
- **CI:** build (Ubuntu/macOS), Python (Ubuntu/macOS), docs, examples, install, style,
  and coverage all run on push and are green. Windows is `workflow_dispatch`-only
  (not continuously verified).
- **Scope of `0.1.0`:** a large minor bump over the last PyPI release (`0.0.4`,
  ~89 commits) — the World/agent simulation architecture, the ex1–ex17 tutorials, the
  relativistic n-body dynamics (Moyer), the IAU 2006/2000A GCRF↔ITRF chain, and the
  ground-station media / solid-tide corrections all landed in this window.

## Why there is no PyPI distribution right now

`pylupnt` **was** on PyPI (`0.0.3`, `0.0.4`) with real prebuilt wheels
(manylinux2014 + macOS arm64, CPython 3.10–3.12 + PyPy) built by a `cibuildwheel`
workflow. Two things happened since:

1. Commit `fca4dca5` ("Rewrite CI around Pixi") replaced that `cibuildwheel`
   workflow with an **sdist-only** publish job.
2. The library got substantially heavier: **OpenCV**, **Matplot++**, and the
   **Fortran** plasma models (IRI/GCPM) were all added *after* `0.0.4`, on top of
   HDF5 and (for the Cesium viz) boost/crow.

An sdist forces every user to compile that entire native stack at `pip install`
time, which is not a usable distribution — so PyPI publishing has been removed
(the `.github/workflows/pypi.yml` job is deleted; the working `cibuildwheel`
configuration is preserved in git history at
`git show v0.0.4:.github/workflows/pypi.yml`).

## Restoring a Python distribution (future work)

The foundation is already correct: `python/pyproject.toml` uses
`scikit-build-core` + `pybind11`, which is exactly what `cibuildwheel` expects.
Two viable paths:

1. **Minimal wheels via `cibuildwheel`.** Add a CMake option (e.g.
   `LUPNT_PYTHON_MINIMAL`) that builds the bindings **without** the visualization-only
   deps — OpenCV, Matplot++, and crow/cesium (plotting is done Python-side with
   matplotlib/plotly anyway). Keep core numerics + SPICE + Fortran + HDF5. Then
   `CIBW_BEFORE_ALL` installs `gfortran` (+ any remaining system libs) and
   `auditwheel`/`delocate` vendor the shared libraries. This shrinks the wheel back
   toward the tractable `0.0.4` situation. Note: **Matplot++ needs `gnuplot` at
   runtime** (a system executable that cannot live in a wheel), so it must be excluded
   from any wheel build regardless.
2. **conda-forge.** For this dependency stack (OpenCV, HDF5, boost, gfortran are all
   trivial conda packages) conda-forge is arguably the better-fit primary channel and
   sidesteps the manylinux vendoring entirely.

Large data (SPICE kernels, LOLA DEM, plasma model tables) must **not** ship in the
package — they are downloaded on demand (some require a NASA Earthdata login). This
is already the pattern.

## Toward v1.0.0

`1.0.0` is a stable-public-API (SemVer) commitment and is **not** recommended yet:

- The API churned heavily in the `0.0.4 → 0.1.0` window (agent/app architecture,
  class renames, config-schema and unit-convention changes). Let it settle across a
  `0.x` cycle or two under real usage before locking it.
- Several underlying methods cite publications still under review.

Before cutting `1.0.0`:

- Restore a real Python distribution (wheels and/or conda-forge) so `pip`/`conda`
  install "just works".
- Decide Windows: put it in continuous CI or explicitly mark it unsupported.
- Provide a zero-credential quickstart (the "hello world" path should run without a
  NASA Earthdata login).
- Add a curated `CHANGELOG` and a statement of what the *public* (supported) API is
  vs. internal, plus migration notes for the recent renames.

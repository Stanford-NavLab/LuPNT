# Running Tests

Run tests from the repository root with pixi so the compiler, library paths,
data paths, and Python path all come from the same environment.

A NASA Earthdata Login (`~/.netrc`, see the main [README](../../README.md#prerequisites)) is
required for the full suite, since several tests fetch or parse live CDDIS GNSS/EOP products. Run
`pixi run download-gnss-test-data` once beforehand to pre-fetch the SP3/BRDC fixture files that
`interfaces.sp3_loader`, `interfaces.rinex_nav_loader`, and
`agents.gnss_constellation.setup_from_files` expect on disk. Without Earthdata Login, use
`pixi run test-cpp-ci` instead (same command CI runs; excludes the CDDIS-dependent tests).

## C++ Tests

```bash
pixi run test-cpp
```

This configures `cpp/test`, builds the `lupnt_tests` executable in
`build-cpp-test-pixi/`, and runs the discovered Catch2 tests with CTest.

To build without running:

```bash
pixi run build-test-cpp
```

To rerun an already-built test binary directly:

```bash
pixi run ctest --test-dir build-cpp-test-pixi --output-on-failure
```

## Python Tests

```bash
pixi run test-py
```

This first builds and copies the `_pylupnt` extension into `python/pylupnt/`,
then runs `pytest python/test`.

## All Tests

```bash
pixi run test
```

If another conda or Anaconda environment is active and CMake finds headers from
that environment, deactivate it before running the pixi task.

## Cross-validation against Orekit and GMAT

The tests under `cpp/test/orekit/` and `cpp/test/gmat/` compare LuPNT's
time-scale conversions, frame conversions, ephemerides, and orbit propagation
against reference values pre-computed with
[Orekit](https://www.orekit.org/) and
[GMAT](https://sourceforge.net/projects/gmat/) (R2022a), checked into
`cpp/test/orekit/data/orekit_reference.json` and
`cpp/test/gmat/data/gmat_reference.json` respectively.

These tests run as part of `pixi run test-cpp` like any other test and **do
not require Orekit/Java or GMAT**. Those tools are only needed by developers
who want to *regenerate* the reference fixtures (e.g. to add new
cross-validation cases) -- see [`cpp/test/orekit/README.md`](orekit/README.md)
and [`cpp/test/gmat/README.md`](gmat/README.md) for those workflows, and the
documentation page `docs/pages/cross_validation.rst` for the detailed
comparison results.

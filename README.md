# LuPNT

[![MacOS](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/macos.yml/badge.svg?branch=development)](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/macos.yml)
[![Windows](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/windows.yml/badge.svg?branch=development)](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/windows.yml)
[![Ubuntu](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/ubuntu.yml/badge.svg?branch=development)](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/ubuntu.yml)
[![Style](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/style.yml/badge.svg?branch=development)](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/style.yml)
[![Install](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/install.yml/badge.svg?branch=development)](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/install.yml)
[![Python](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/python.yml/badge.svg?branch=development)](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/python.yml)
[![Examples](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/examples.yml/badge.svg?branch=development)](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/examples.yml)
[![Documentation Status](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/docs.yaml/badge.svg?branch=development)](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/docs.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI Release](https://img.shields.io/pypi/v/pylupnt.svg)](https://pypi.org/project/pylupnt)
[![Python Versions](https://img.shields.io/pypi/pyversions/pylupnt)](https://pypi.org/project/pylupnt)
[![codecov](https://codecov.io/gh/Stanford-NavLab/LuPNT/branch/development/graph/badge.svg)](https://codecov.io/gh/Stanford-NavLab/LuPNT)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Stanford-NavLab/LuPNT/development?labpath=python%2Fexamples%2Fex_frozen_orbits%2Fex_frozen_orbits.ipynb)
[![Open In Collab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1yhHImp3hB5P8dadLlv0CcYWQWLabNfuh)


`LuPNT` is an open-source C++/Python library for Lunar Positioning, Navigation, and Timing (PNT) research. It provides high-fidelity astrodynamics, signal propagation models, and navigation algorithms tailored for cislunar missions. This project is a product of the [Stanford NAV Lab](https://navlab.stanford.edu/).

If using this project in your own work please cite the following:

```bibtex
@inproceedings{IiyamaCasadesus2023,
  title = {LuPNT: Open-Source Simulator for Lunar Positioning, Navigation, and Timing},
  author={Iiyama, Keidai and Casadesus Vila, Guillem and Gao, Grace},
  booktitle={Proceedings of the Institute of Navigation Gnss+ conference (ION Gnss+ 2023)},
  institution = {Stanford University},
  year = {2023},
  url = {https://github.com/Stanford-NavLab/LuPNT},
}
```

---

## Features

| Module | Description |
|--------|-------------|
| **Dynamics** | Two-body, N-body, and high-fidelity force models (gravity, drag, SRP) for Earth and lunar orbits |
| **Environment** | Gravity fields, atmosphere, solar system bodies, occultation, and plasma/ionosphere models |
| **Plasma** | GCPM v2.4 plasmasphere + IRI ionosphere electron density; GNSS ray-tracing with TEC and signal delay for cislunar links (via integrated pecsim) |
| **GNSS** | GPS/GNSS signal generation, SP3/ANTEX-backed constellations, and space-user ranging |
| **Measurements** | Pseudorange, Doppler, and carrier-phase measurement models with light-time, Shapiro, and optional plasma corrections |
| **Filters** | Extended Kalman Filter (EKF), Unscented Kalman Filter (UKF), and batch least-squares estimators |
| **Conversions** | Reference frame transformations (ECI, ECEF, LVLH, Moon-centered) and time system utilities |
| **Numerics** | Numerical integration (RK4, Euler), Nelder-Mead optimizer, and matrix utilities |
| **Python bindings** | Full `pylupnt` Python package exposing the C++ library via pybind11 |

---

## Repository Structure

```
LuPNT/
├── cpp/
│   ├── lupnt/                  # C++ library source
│   │   ├── agents/             # Agent/satellite abstractions
│   │   ├── conversions/        # Frame and unit conversions
│   │   ├── core/               # Logging, constants, math utils
│   │   ├── dynamics/           # Orbital dynamics and propagators
│   │   ├── environment/        # Gravity, atmosphere, solar system, plasma
│   │   │   └── plasma/         # GCPM v2.4 + IRI + ray-tracer (pecsim)
│   │   ├── filters/            # Navigation filters (EKF, UKF, …)
│   │   ├── measurements/       # Measurement models
│   │   └── numerics/           # Numerical methods
│   └── examples/               # C++ example programs
│       ├── environment/        # Plasma, solar system examples
│       ├── dynamics/           # Orbit propagation examples
│       └── …
├── python/
│   ├── pylupnt/                # Python package (installed in-place)
│   │   ├── plasma/             # pylupnt.plasma sub-module
│   │   ├── plot/               # Plotting utilities
│   │   └── _pylupnt.so         # Compiled pybind11 extension
│   └── bindings/               # pybind11 binding source files
├── projects/                   # Research project notebooks and scripts
│   ├── GNSS_Filtering/         # Staged lunar GNSS filtering simulation
│   ├── Plasma_Examples/        # GCPM, IRI, and ray-tracing Jupyter notebooks
│   └── …
├── data/
│   └── LuPNT_data/             # Runtime data (ephemeris, TLE, plasma coefficients)
│       ├── ephemeris/
│       ├── plasma/             # IRI coefficient files, Kp index data
│       └── tle/
└── pixi.toml                   # Reproducible environment and task definitions
```

---

## Development

### Prerequisites

Install [pixi](https://pixi.sh) (manages the compiler toolchain, Python, and all C++ dependencies via conda-forge — no manual dependency installation needed):
```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

### Setup

1. **Install the environment** (downloads all dependencies into an isolated conda environment)
```bash
pixi install
```

2. **Build the C++ library**
```bash
pixi run build
```

3. **Build and deploy the Python bindings**
```bash
pixi run build-py
```

### Build modes

Three CMake build types are available, each using its own build directory so you can switch between them without reconfiguring:

| Pixi task | CMake type | Build directory | Use case |
|-----------|-----------|-----------------|----------|
| `pixi run build` | `Release` | `build/` | Production — full optimisation, no debug symbols |
| `pixi run build-reldbg` | `RelWithDebInfo` | `build-reldbg/` | Profiling — optimised with debug symbols |
| `pixi run build-debug` | `Debug` | `build-debug/` | Development — no optimisation, full debug symbols |

Corresponding Python-binding tasks copy the compiled `.so` into `python/pylupnt/`:

```bash
pixi run build-py          # Release bindings
pixi run build-py-reldbg   # RelWithDebInfo bindings
pixi run build-py-debug    # Debug bindings
```

### Daily workflow

Activate the pixi environment in a shell (sets `LUPNT_DATA_PATH`, `PECSIMPY_BASE_PATH`, and `PYTHONPATH` automatically):
```bash
pixi shell
```

Or run a single command without entering a shell:
```bash
pixi run <command>
```

### Running C++ examples

After building, examples are in `build/examples/`. Run any example with:
```bash
# From inside pixi shell
./build/examples/ex_plasma_gcpm
./build/examples/ex_plasma_raytrace
```

Environment variables are set automatically inside `pixi shell`. Outside of it, set them explicitly:
```bash
LUPNT_DATA_PATH=$PWD/data/LuPNT_data \
PECSIMPY_BASE_PATH=$PWD/data/LuPNT_data/plasma \
./build/examples/ex_plasma_gcpm
```

### Jupyter notebooks

Select the **`lupnt (pixi)`** kernel in Jupyter. The kernel has `LUPNT_DATA_PATH`, `PECSIMPY_BASE_PATH`, and `PYTHONPATH` pre-configured so all notebooks in `projects/` work out of the box.

### Running tests

```bash
pixi run test-cpp        # Build and run C++ Catch2/CTest tests
pixi run test-py         # Build Python bindings, then run pytest
pixi run test            # Run both suites
```

More details are in [cpp/test/README.md](cpp/test/README.md).

### GNSS filtering example

The lunar GNSS filtering project under
[`projects/GNSS_Filtering/`](projects/GNSS_Filtering/README.md) demonstrates
the current high-fidelity workflow:

1. Precompute receiver-GNSS links in C++ from SP3/ANTEX data.
2. Precompute ionospheric/plasmaspheric delays in Python with process-parallel
   GCPM ray tracing.
3. Run the C++ Monte Carlo filter with light-time, Shapiro delay, relativistic
   clock dynamics, OCXO clock noise, and optional SRP coefficient estimation.
4. Post-process plasma delays, tracked satellites, and RTN state errors with
   covariance bounds.

Run the full staged pipeline with:

```bash
pixi run run-gnss-pipeline
```

To skip the expensive GCPM batch when the delay table already exists, or when
using a no-plasma config:

```bash
pixi run run-gnss-pipeline --skip-delays
```

### Optional: Orekit / GMAT cross-validation (developers)

The C++ tests under [`cpp/test/orekit/`](cpp/test/orekit/README.md) and
[`cpp/test/gmat/`](cpp/test/gmat/README.md) cross-check LuPNT's time-scale
conversions, sidereal angles, frame conversions (Earth and Moon), Sun/Moon
ephemerides, and orbit propagation against reference values from
[Orekit](https://www.orekit.org/) (a JVM-based astrodynamics library) and
[GMAT](https://sourceforge.net/projects/gmat/) (NASA's General Mission
Analysis Tool). To keep LuPNT lightweight for normal users, these tests do
**not** depend on Orekit/Java or GMAT at all:

- The reference values were computed once with each tool and checked into
  `cpp/test/orekit/data/orekit_reference.json` and
  `cpp/test/gmat/data/gmat_reference.json` by the developer-only generators
  [`gen_orekit_reference.py`](cpp/test/orekit/gen_orekit_reference.py) and
  [`gen_gmat_reference.py`](cpp/test/gmat/gen_gmat_reference.py).
- The `*_orekit_reference`/`*_gmat_reference` tests load those JSON fixtures
  and compare them against LuPNT's own computations, with no runtime
  Orekit/GMAT dependency. They run automatically as part of
  `pixi run test-cpp`.
- Orekit is only declared as a dependency of the optional `dev` pixi
  environment, so `pixi install` (the default environment) never downloads
  it. GMAT is a separately installed desktop application and is not, and
  will not be, a pixi/conda dependency.

The detailed comparison results -- what is compared, the measured
differences per category, and the known modeling differences the comparison
uncovered -- are documented in
[`docs/pages/cross_validation.rst`](docs/pages/cross_validation.rst).

If you want to **regenerate** the reference vectors (e.g. to add new
cross-validation cases or update them after an intentional algorithm
change):

```bash
# Orekit (installs the optional dev environment with the orekit package)
pixi install -e dev
pixi run -e dev python cpp/test/orekit/gen_orekit_reference.py

# GMAT (requires a local GMAT install; GMAT_HOME also works, and common
# default install locations are searched automatically)
GMAT_CONSOLE=/Applications/GMAT_R2022a/bin/GmatConsole-R2022a \
    python3 cpp/test/gmat/gen_gmat_reference.py
```

Each generator overwrites its fixture; review the diff and re-run
`pixi run test-cpp` to confirm the comparison tests still pass. See
[`cpp/test/orekit/README.md`](cpp/test/orekit/README.md) and
[`cpp/test/gmat/README.md`](cpp/test/gmat/README.md) for the full data
dictionaries, tool-specific caveats, and tolerance rationale.

### Building documentation

```bash
pixi run build-docs
```

This generates the Sphinx/Doxygen documentation site from the existing headers
and Python binding, then copies the GitHub Pages-ready HTML into `build/docs/`.
If `python/pylupnt/_pylupnt*.so` is missing, run `pixi run build-py` first.
The generated folder includes `.nojekyll` so Pages serves Sphinx assets as-is.
More details are in [docs/builddocs.rst](docs/builddocs.rst).

To preview the docs page locally, run the following command

```bash
python -m http.server 8000 --directory build/docs
```
and open
```
http://localhost:8000
```

### Pre-commit

Pre-commit hooks run automatically on `git commit`. To run manually:
```bash
pixi run pre-commit run --all-files
```

---

## Third-party components

LuPNT's plasma module includes the following scientific models, which are freely
available for scientific use with attribution. The MIT license for LuPNT applies
to LuPNT's own source code only; these components retain their original terms.

| Component | Description | Source |
|-----------|-------------|--------|
| **IRI2007 / IRI2020** | International Reference Ionosphere electron density model | [irimodel.org](https://irimodel.org) |
| **GCPM v2.4** | Global Core Plasma Model (D.L. Gallagher, NASA MSFC) | [plasmasphere.nasa.gov](https://plasmasphere.nasa.gov) |
| **xform** | Coordinate transformation utilities | [plasmasphere.nasa.gov](https://plasmasphere.nasa.gov) |
| **IGRF14** | International Geomagnetic Reference Field (IAGA Working Group V-MOD / BGS) | [ngdc.noaa.gov](https://www.ngdc.noaa.gov/IAGA/vmod/igrf14.html) |

Publications using these components should cite the relevant papers:

- **IRI2020**: Bilitza, D., et al. (2022). *The International Reference Ionosphere model: A review and description of an ionospheric benchmark*. Reviews of Geophysics, 60, e2022RG000792. https://doi.org/10.1029/2022RG000792
- **GCPM**: Gallagher, D. L., et al. (2000). *An empirical model of the Earth's plasmasphere*. Advances in Space Research, 25(12), 2421–2430.
- **IGRF14**: Alken, P., et al. (2021). *International Geomagnetic Reference Field: the thirteenth generation*. Earth, Planets and Space, 73, 49. https://doi.org/10.1186/s40623-020-01288-x

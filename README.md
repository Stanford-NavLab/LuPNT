<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/lupnt_logo_horizontal_dark.svg">
    <img alt="LuPNT — Lunar Positioning, Navigation and Timing" src="assets/lupnt_logo_horizontal.svg" width="520">
  </picture>
</p>

[![MacOS](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/macos.yml/badge.svg?branch=development)](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/macos.yml)
[![Ubuntu](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/ubuntu.yml/badge.svg?branch=development)](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/ubuntu.yml)
[![Style](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/style.yml/badge.svg?branch=development)](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/style.yml)
[![Install](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/install.yml/badge.svg?branch=development)](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/install.yml)
[![Python](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/python.yml/badge.svg?branch=development)](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/python.yml)
[![Examples](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/examples.yml/badge.svg?branch=development)](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/examples.yml)
[![Documentation Status](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/docs.yaml/badge.svg?branch=development)](https://github.com/Stanford-NavLab/LuPNT/actions/workflows/docs.yaml)
[![codecov](https://codecov.io/gh/Stanford-NavLab/LuPNT/branch/development/graph/badge.svg)](https://codecov.io/gh/Stanford-NavLab/LuPNT)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


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

LuPNT is organized as a **config-driven, agent-based simulation framework**. *Agents*
(satellites, ground stations, rovers, landers, surface stations, constellations) host
*Applications* — the mission and navigation-filter logic — and run together on a shared,
event-scheduled *Simulation*. Agents, their devices, dynamics, and applications are all
assembled from YAML configuration, so new scenarios are described in config rather than code.

| Module | Description |
|--------|-------------|
| **Agents & Applications** | Config-driven agents (satellite, ground station, rover, lander, surface station, constellation) that host navigation/mission **Applications** — ODTS, ISL, surface-rover and lander navigation, ephemeris/almanac generation — built from YAML via an asset factory |
| **Simulations** | Event-scheduled `Simulation` engine plus end-to-end scenario drivers (ground-station & inter-satellite ODTS, lunar GNSS ODTS, lander/surface navigation, ephemeris fitting) |
| **Devices & sensors** | Clocks, IMUs, cameras, GNSS receivers, and communication devices attached to agents |
| **Dynamics** | Two-body, N-body, and high-fidelity force models (gravity, drag, SRP) for Earth, lunar, and arbitrary central-body orbits |
| **Environment** | Gravity fields, atmosphere, solar system bodies, occultation, and plasma/ionosphere models |
| **Plasma** | GCPM v2.4 plasmasphere + IRI ionosphere electron density; GNSS ray-tracing with TEC and signal delay for cislunar links (via integrated pecsim) |
| **GNSS** | GPS/GNSS signal generation, SP3/ANTEX-backed constellations, yaw-steering/attitude models, and space-user ranging |
| **Measurements** | Pseudorange, Doppler, and carrier-phase measurement models with light-time, Shapiro, and optional plasma corrections |
| **Filters & States** | Extended Kalman Filter (EKF), Unscented Kalman Filter (UKF), square-root information filter (SRIF)/smoother, and batch least-squares estimators over composable state/joint-state abstractions |
| **Interfaces** | Data loaders and I/O: SP3, ANTEX, RINEX nav, TLE, SPICE kernels, EOP/TAI-UTC, LOLA DEM and crater data, plus Cesium and Matplotlib export |
| **Conversions** | Reference frame transformations (ECI, ECEF, LVLH, Moon-centered, generic body-fixed/inertial) and time system utilities |
| **Numerics** | Numerical integration (RK4, Euler), Nelder-Mead optimizer, and matrix utilities |
| **Visualization** | Matplotlib/Plotly plotting plus interactive 3-D [CesiumJS](https://cesium.com/platform/cesiumjs/) scenes of constellations and surface assets (`pnt.plot.CesiumScene`) |
| **Python bindings & authoring** | Full `pylupnt` Python package exposing the C++ library via pybind11 — including subclassing `Application`, `Agent`, and `Measurement` in pure Python and registering them into a YAML-driven `Simulation` (see Example 17) |

---

## Repository Structure

```
LuPNT/
├── cpp/
│   ├── lupnt/                  # C++ library source
│   │   ├── agents/             # Agents: satellite, ground station, rover, lander, surface station, constellations
│   │   ├── applications/       # Per-agent mission/nav logic (ODTS, ISL, lander/rover nav, ephemeris/almanac gen)
│   │   ├── simulations/        # Event-scheduled Simulation engine + scenario drivers
│   │   ├── devices/            # Agent devices: clock, IMU, camera, GNSS receiver, comms
│   │   ├── dynamics/           # Orbital dynamics and propagators
│   │   ├── environment/        # Gravity, atmosphere, solar system, plasma
│   │   │   └── plasma/         # GCPM v2.4 + IRI + ray-tracer (pecsim)
│   │   ├── measurements/       # Measurement models
│   │   ├── states/             # State / joint-state / parameter abstractions
│   │   ├── transmission/       # Signal transmission modeling
│   │   ├── interfaces/         # Data loaders & I/O (SP3, ANTEX, RINEX, TLE, SPICE, EOP, DEM, Cesium, matplotlib)
│   │   ├── conversions/        # Frame and unit conversions
│   │   ├── core/               # Logging, constants, config, math utils
│   │   └── numerics/           # Numerical methods
│   └── examples/               # C++ example programs
│       ├── tutorials/          # Tutorial counterparts to the Python notebooks
│       ├── simulations/        # Full scenario applications
│       ├── environment/        # Plasma, solar system examples
│       ├── dynamics/           # Orbit propagation examples
│       └── …
├── configs/                    # YAML scenario configuration (agents, applications, dynamics, environments, datasets)
├── python/
│   ├── pylupnt/                # Python package (installed in-place)
│   │   ├── plasma/             # pylupnt.plasma sub-module
│   │   ├── plot/               # Plotting utilities
│   │   ├── interfaces/         # Python-side data interfaces
│   │   └── _pylupnt.so         # Compiled pybind11 extension
│   ├── examples/               # Tutorial notebooks (ex1–ex17) — see python/examples/README.md
│   └── bindings/               # pybind11 binding source files
├── projects/                   # Research project notebooks and scripts
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

## Tutorials & Examples

A progressive set of Jupyter notebooks in [`python/examples/`](python/examples/README.md) teaches
the `pylupnt` API and reproduces LuPNT's core navigation workflows — from single-orbit propagation
up to full orbit-determination filters, constellation design, optical navigation, and non-lunar
(Mars / LEO) PNT. See the [examples README](python/examples/README.md) for the full index and
per-notebook details.

| # | Example | Theme |
|---|---------|-------|
| 1 | [Propagating an ELFO orbit](python/examples/ex1_propagate_orbit.ipynb) | Fundamentals — elements, frames, force models |
| 2 | [Relativistic time conversions](python/examples/ex2_time_conversions.ipynb) | Fundamentals — TT/TCG/TCB/TDB/TCL/LT |
| 3 | [GNSS interface tutorial](python/examples/ex3_gnss_interface.ipynb) | Earth-GNSS data — TLEs, antenna gain, SP3 vs. broadcast |
| 4 | [Plasmasphere & ionosphere](python/examples/ex4_plasmasphere.ipynb) | GCPM v2.4 electron density, TEC, ray-traced signal delay |
| 5 | [GNSS measurement simulation](python/examples/ex5_gnss_measurement_sim.ipynb) | Lunar-receiver visibility and C/N₀ from precise ephemerides |
| 6 | [Sidelobe pseudorange/Doppler/TDCP ODTS](python/examples/ex6_gnss_odts.ipynb) | Weak-signal EKF for position, velocity, clock, and SRP |
| 7 | [Ground-station orbit determination](python/examples/ex7_groundstation_odts.ipynb) | DSN batch least squares + SRIF/smoother |
| 8 | [Distributed ISL ODTS](python/examples/ex8_isl_odts.ipynb) | 5 parallel consider-state EKFs + surface-station timing |
| 9 | [Ephemeris & almanac fitting](python/examples/ex9_ephemeris.ipynb) | Compressing trajectories into broadcast navigation models |
| 10 | [Surface rover navigation](python/examples/ex10_surface_rover.ipynb) | IMU + LCRNS pseudoranges + DEM error-state EKF |
| 11 | [Lunar lander navigation](python/examples/ex11_lander_navigation.ipynb) | Powered-descent MEKF: IMU, altimeter, craters, LunaNet |
| 12 | [Constellation design](python/examples/ex12_constellation_design.ipynb) | ELFO Walker layout, PDOP/coverage/EIRP, phasing optimization |
| 13 | [Cesium visualization](python/examples/ex13_cesium.ipynb) | Interactive 3-D CesiumJS scenes of constellations & stations |
| 14 | [Optical navigation](python/examples/ex14_opnav.ipynb) | Lunar-horizon image processing → EKF position fixes |
| 15 | [Mars PNT](python/examples/ex15_marspnt.ipynb) | Beyond the Moon — Mars gravity, frames, Walker constellation |
| 16 | [LEO PNT](python/examples/ex16_leopnt.ipynb) | Earth LEO constellation with Harris-Priester atmospheric drag |
| 17 | [New simulation in Python](python/examples/ex17_python_new_sim_example.ipynb) | Authoring a new Application + Measurement in pure Python (angles-only OD) |

C++ tutorial counterparts live in [`cpp/examples/tutorials/`](cpp/examples/tutorials/).

---

## Development

### Prerequisites

1. Install [pixi](https://pixi.sh) (manages the compiler toolchain, Python, and all C++ dependencies via conda-forge — no manual dependency installation needed):
```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

2. **NASA Earthdata Login is required.** Several C++ tests (`pixi run test-cpp`) and GNSS
   examples depend on live GNSS/EOP products (SP3 precise orbits, RINEX nav, ANTEX) fetched from
   [CDDIS](https://cddis.nasa.gov/), which requires authentication. Set this up before running
   `pixi run test-cpp` or any of the GNSS download tasks:

   1. **Create an Earthdata Login account** at
      [urs.earthdata.nasa.gov](https://urs.earthdata.nasa.gov/): click "Register for a profile",
      fill in the required fields, and verify your email address.

   2. **Create a `~/.netrc` file** with your credentials (run inside WSL on Windows, not
      PowerShell/cmd — see [Setting up on Windows (via WSL)](#setting-up-on-windows-via-wsl)
      below):
   ```bash
   touch ~/.netrc
   chmod 0600 ~/.netrc
   echo "machine urs.earthdata.nasa.gov login <your_username> password <your_password>" >> ~/.netrc
   ```
      Replace `<your_username>` and `<your_password>` with your actual credentials. The
      `chmod 0600` step is required — curl/wget refuse to use a `.netrc` with broader permissions.
      Verify with `cat ~/.netrc`. Never share or commit this file (it stays outside the repo, at
      `~/.netrc`, so there's no risk of it being checked into git).

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

> **Note:** `data/LuPNT_data` (ephemeris, GNSS antenna/clock products, plasma coefficients, TLEs)
> is downloaded and extracted automatically the first time you run step 2 or 3 above — no manual
> setup needed. It's a ~650MB one-time download (see
> [cmake/FetchLuPNTData.cmake](cmake/FetchLuPNTData.cmake)); subsequent builds skip it since the
> data is already present. To disable this (e.g. air-gapped builds), pass
> `-DLUPNT_FETCH_DATA=OFF` to the underlying `cmake` configure call, or place the data manually
> at `data/LuPNT_data` beforehand.

### Setting up on Windows (via WSL)

LuPNT does not build natively on Windows; Windows users should develop inside
**WSL (Windows Subsystem for Linux)**, which runs a real Linux environment and
follows the same setup as native Linux/macOS above.

1. **Install WSL** (PowerShell, run as Administrator). This installs the
   default Ubuntu distribution:
```powershell
wsl --install
```
   Reboot if prompted, then launch "Ubuntu" from the Start menu and create a
   Unix username/password when asked.

2. **Update packages inside the WSL Ubuntu shell**
```bash
sudo apt update && sudo apt upgrade -y
```

3. **Clone the repository inside the WSL filesystem** (not under `/mnt/c/...`,
   which is much slower and can break symlinks/permissions):
```bash
cd ~
git clone https://github.com/Stanford-NavLab/LuPNT.git
cd LuPNT
```

4. **Continue with the regular [Prerequisites](#prerequisites) and
   [Setup](#setup) steps above** (`curl -fsSL https://pixi.sh/install.sh | bash`,
   then `pixi install`, `pixi run build`, `pixi run build-py`) from inside the
   WSL shell.

5. **Editing files**: you can edit the cloned repo from Windows tools (e.g.
   VS Code with the "WSL" extension, or `code .` from inside the WSL shell)
   while all builds/tests run inside WSL.

Once set up, every command in the rest of this README (`pixi shell`,
`pixi run test`, etc.) works identically inside the WSL terminal.

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

Activate the pixi environment in a shell (sets `LUPNT_DATA_PATH`, `LUPNT_OUTPUT_PATH`, `PECSIMPY_BASE_PATH`, and `PYTHONPATH` automatically):
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

First register the Jupyter kernel (once per machine, and after moving the repo):

```bash
pixi run install-kernel
```

Then select the **`LuPNT (pixi)`** kernel in Jupyter / VS Code. `install-kernel` bakes
`PYTHONPATH`, `LUPNT_DATA_PATH`, `LUPNT_OUTPUT_PATH`, and `PECSIMPY_BASE_PATH` into the
kernelspec (derived from the repo layout, so it is correct on both macOS and Linux/WSL) —
VS Code and Jupyter launch the kernel's interpreter *without* pixi activation, so these must
live in the kernelspec itself for `import pylupnt` and the `projects/` notebooks to work out
of the box. If the picker only shows the bare pixi interpreter (no env vars), re-run the task.

### Interactive Cesium visualization

`pnt.plot.CesiumScene` renders constellations, relay satellites, and surface stations as an
interactive 3-D [CesiumJS](https://cesium.com/platform/cesiumjs/) scene — inline in Jupyter
and as a standalone `.html` you can open in any browser or host to share. See
[`python/examples/ex13_cesium.ipynb`](python/examples/ex13_cesium.ipynb) for a full tutorial
(GNSS + lunar relay constellations + surface stations).

```python
import pylupnt as pnt
from datetime import datetime

scene = pnt.plot.CesiumScene(body="MOON", name="Lunar relays", epoch=datetime(2027, 3, 1))
scene.add_trajectory("SV-1", t_tdb, rv_mci, frame_in=pnt.MOON_CI)  # inertial -> Moon-fixed
scene.add_station("Shackleton", lat=-89.9, lon=0.0)                # degrees
scene.show()                                                       # writes ./cesium_scenes/*.html
```

**Do I need a Cesium ion account?** **No.** Scenes are fully self-contained: the central body
(Earth or Moon) is drawn as an ellipsoid textured with LuPNT's own surface imagery, embedded
directly in the generated HTML — no ion access token, login, or billing. The only thing
fetched at view time is the CesiumJS library from a public CDN (`cdn.jsdelivr.net`), so an
internet connection is needed to view a scene. (If you already have an ion token and want
streamed high-resolution terrain/imagery, you can edit the generated `.html`, but nothing
here requires it.)

### Running tests

`pixi run test-cpp` includes tests (`data.eop_latest_iers`, `devices.space_comms`,
`agents.gnss_constellation.setup_from_files`, `interfaces.antex_loader`,
`interfaces.rinex_nav_loader`, `interfaces.sp3_loader`, `measurements.gnss_measurements.visibility`,
`states.tle`) that need live GNSS/EOP products (SP3, RINEX nav, ANTEX) from
[CDDIS](https://cddis.nasa.gov/) — this is why NASA Earthdata Login (see
[Prerequisites](#prerequisites)) is required. Some of these tests (`interfaces.sp3_loader`,
`interfaces.rinex_nav_loader`, `agents.gnss_constellation.setup_from_files`) additionally expect
specific SP3/BRDC files to already exist on disk rather than downloading them on the fly (see
`GnssFilesDir()` in [cpp/test/interfaces/test_sp3_loader.cc](cpp/test/interfaces/test_sp3_loader.cc)),
so fetch those once before running the suite:

```bash
pixi run download-gnss-test-data   # one-time; requires ~/.netrc (see Prerequisites)
pixi run test-cpp                  # Build and run C++ Catch2/CTest tests
pixi run test-py                   # Build Python bindings, then run pytest
pixi run test                      # Run both suites
```

If you don't have Earthdata Login set up, use `pixi run test-cpp-ci` instead — it's the same
command CI runs and excludes the tests above.

More details are in [cpp/test/README.md](cpp/test/README.md).

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

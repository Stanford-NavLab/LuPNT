# LuPNT Python Examples

A progressive set of Jupyter notebooks that teach the `pylupnt` API and reproduce the
navigation workflows LuPNT is built for. They start from single-orbit propagation and build up
to full orbit-determination filters, constellation design, optical navigation, non-lunar
(Mars / LEO) PNT systems, and authoring a whole new scenario in pure Python.

## Running the notebooks

1. Build the Python bindings once: `pixi run build-py` (from the repo root).
2. Register and select the Jupyter kernel:
   ```bash
   pixi run install-kernel   # once per machine
   ```
   Then pick the **`LuPNT (pixi)`** kernel in Jupyter or VS Code. The kernelspec bakes in
   `PYTHONPATH`, `LUPNT_DATA_PATH`, `LUPNT_OUTPUT_PATH`, and `PECSIMPY_BASE_PATH`, so
   `import pylupnt as pnt` and the data-backed examples work out of the box.

> **Data & credentials.** Examples 3–6 use live GNSS/EOP products (SP3, RINEX nav, ANTEX) from
> [CDDIS](https://cddis.nasa.gov/), which require a free
> [NASA Earthdata Login](https://urs.earthdata.nasa.gov/) and a `~/.netrc` file — see the
> [Prerequisites](../../README.md#prerequisites) in the main README. The remaining examples run
> from the bundled `data/LuPNT_data` (downloaded automatically on first build).

## The examples

### Fundamentals

| # | Notebook | What it covers |
|---|----------|----------------|
| 1 | [ex1_propagate_orbit.ipynb](ex1_propagate_orbit.ipynb) | The core LuPNT workflow: build an Elliptical Lunar Frozen Orbit (ELFO) from classical elements, convert to Cartesian, rotate between lunar frames, propagate with a configurable force model, and convert back to elements. |
| 2 | [ex2_time_conversions.ipynb](ex2_time_conversions.ipynb) | Relativistic time scales for cislunar navigation (TT, TCG, TCB, TDB, TCL, LT). Uses LuPNT's conversion routines and reproduces the secular/periodic terms with a vectorized Python model. |

### Earth-GNSS interfaces and the cislunar link

| # | Notebook | What it covers |
|---|----------|----------------|
| 3 | [ex3_gnss_interface.ipynb](ex3_gnss_interface.ipynb) | LuPNT's Earth-GNSS data interfaces: load GPS/Galileo/QZSS from TLEs, inspect transmit antenna gain patterns (main-lobe vs. sidelobe), and compare RINEX broadcast vs. IGS precise products. |
| 4 | [ex4_plasmasphere.ipynb](ex4_plasmasphere.ipynb) | The plasma environment: sample the GCPM v2.4 electron-density model on a meridional grid, then ray-trace a GPS-to-lunar link to compute total electron content (TEC) and dispersive signal delay. |
| 5 | [ex5_gnss_measurement_sim.ipynb](ex5_gnss_measurement_sim.ipynb) | Simulate which Earth GNSS signals a lunar receiver can track: propagate an ELFO receiver, load precise ephemerides with `SP3Loader`, build a `GnssConstellation`, and run `GNSSMeasurements.precompute()` for visibility and C/N₀ histories. |
| 6 | [ex6_gnss_odts.ipynb](ex6_gnss_odts.ipynb) | Sidelobe pseudorange + Doppler + TDCP orbit determination and time sync (ODTS) with a stochastic-cloning UDU EKF. Estimates position, velocity, clock bias/drift, and an SRP coefficient from weak Earth-GNSS sidelobe signals. Helper scripts: [`ex6_gnss_odts_config.py`](ex6_gnss_odts_config.py), [`ex6_precompute.py`](ex6_precompute.py), [`ex6_run_gnss_odts.py`](ex6_run_gnss_odts.py). |

### Orbit determination and time synchronization

| # | Notebook | What it covers |
|---|----------|----------------|
| 7 | [ex7_groundstation_odts.ipynb](ex7_groundstation_odts.ipynb) | Earth-based tracking of a lunar satellite from three Deep Space Network complexes (Goldstone, Canberra, Madrid). Batch least squares with an analytic (STM-chained) design matrix, refined by a square-root information filter (SRIF) and smoother. |
| 8 | [ex8_isl_odts.ipynb](ex8_isl_odts.ipynb) | **Distributed** onboard ODTS for the 5-satellite LCRNS Reference Constellation 3.1. Every satellite runs its own Schmidt (consider-state) EKF in parallel, fusing two-way crosslink range/range-rate and a one-way pseudorange from a rotating lunar surface station. |
| 9 | [ex9_ephemeris.ipynb](ex9_ephemeris.ipynb) | Compress a numerically propagated trajectory into broadcast navigation models: `pnt.CartesianEphemeris` (Chebyshev residuals) vs. `pnt.Almanac` (element polynomials + Fourier terms), trading broadcast bits against fit error. |

### Surface and terminal navigation

| # | Notebook | What it covers |
|---|----------|----------------|
| 10 | [ex10_surface_rover.ipynb](ex10_surface_rover.ipynb) | A south-pole rover using a strapdown IMU aided by LunaNet/LCRNS pseudoranges and a digital elevation model (DEM). Error-state EKF over position, velocity, attitude error, IMU biases, and clock terms; the DEM constraint sharply cuts vertical drift. |
| 11 | [ex11_lander_navigation.ipynb](ex11_lander_navigation.ipynb) | Powered-descent lander navigation with a multiplicative EKF (MEKF): quaternion attitude plus IMU, nadir altimeter, crater-bearing landmarks, and LunaNet pseudoranges, hosted as a `LanderNavApp` on a `Lander` agent. |

### Constellation design and visualization

| # | Notebook | What it covers |
|---|----------|----------------|
| 12 | [ex12_constellation_design.ipynb](ex12_constellation_design.ipynb) | Design a lunar navigation constellation for south-pole service with a reusable `LunarNavConstellation` class: symmetric ELFO Walker layout, visibility/PDOP, required EIRP, and satellite-phasing optimization. |
| 13 | [ex13_cesium.ipynb](ex13_cesium.ipynb) | Turn trajectory samples into interactive 3-D [CesiumJS](https://cesium.com/platform/cesiumjs/) scenes with `pnt.plot.CesiumScene` — Earth GNSS from TLEs, lunar LCRNS relays in `MOON_PA`, and surface stations — as a debugging and presentation tool. No Cesium ion token required. |
| 14 | [ex14_opnav.ipynb](ex14_opnav.ipynb) | A compact optical-navigation pipeline for a lunar orbiter: synthetic horizon-image generation, disk fitting to angular radius and bearing, and an EKF that turns each image into a Moon-centered position measurement. |

### Generalizing beyond the Moon

| # | Notebook | What it covers |
|---|----------|----------------|
| 15 | [ex15_marspnt.ipynb](ex15_marspnt.ipynb) | LuPNT applied to **Mars**: an 8×8 Mars gravity field (`Mars50c.cof`), native `MARS_CI`/`MARS_FIXED` frames, a 9-satellite 3-plane Walker constellation propagated with `NBodyDynamics`, and a surface-user DOP/positioning map. |
| 16 | [ex16_leopnt.ipynb](ex16_leopnt.ipynb) | The Earth companion to Example 15: a 110-satellite **LEO PNT** Walker constellation at 600 km with Harris-Priester atmospheric drag, showing why LEO providers must model drag, and mapping coverage/PDOP for a user in San Francisco. |

### Extending LuPNT in Python

| # | Notebook | What it covers |
|---|----------|----------------|
| 17 | [ex17_python_new_sim_example.ipynb](ex17_python_new_sim_example.ipynb) | Author a new simulation in **pure Python**: subclass `pnt.Application` and `pnt.Measurement`, register them with `pnt.register_application`, and run an angles-only orbit-determination scenario driven by the C++ `pnt.Simulation`. The smallest template for building your own scenario. |

---

For the C++ tutorial counterparts, see [`cpp/examples/tutorials/`](../../cpp/examples/tutorials/).

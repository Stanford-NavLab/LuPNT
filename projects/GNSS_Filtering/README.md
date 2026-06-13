# Lunar GNSS ODTS Scenario

This scenario configures the LuPNT `LunarGnssODTSSimulation` under
`cpp/lupnt/simulations/LunarGnssODTS`. The project folder intentionally only
contains scenario assets: YAML configs, the Python GCPM delay precompute,
plotting scripts, and the notebook.

The default configuration estimates receiver position, velocity, clock bias,
and clock drift from pseudorange and Doppler measurements. TDCP can be enabled
with `measurements.use_tdcp: true`; this switches the C++ filter to UDU
stochastic cloning so each TDCP row can depend on both the current and previous
receiver orbit/clock state.
It is structured as a staged pipeline so that GNSS link geometry, expensive
GCPM/IRI plasma delays, Monte Carlo filtering, and post-processing can be run,
cached, and inspected separately.

## Four-stage pipeline

1. Load receiver hardware, tracking-loop, link-budget, and measurement-noise
   parameters from a named design in `gnss_designs.yaml`, with optional per-run
   overrides in `gnss_filtering_config.yaml`.
2. Build GNSS constellations from SP3/ANTEX files. With
   `constellation.auto_select_sp3: true`, the simulation scans
   `constellation.sp3_directory` and selects the SP3 files that cover the
   configured start epoch, duration, and ephemeris margin. ANTEX PCO entries are
   selected at the same ephemeris epochs. The default uses every GPS satellite
   found in the selected SP3 file; Galileo can be added with
   `constellation.include_galileo: true` when the SP3/ANTEX files contain those
   satellites.
3. Propagate the truth receiver state on an ELFO with `NBodyDynamics` and a
   coupled `JointOrbitClockDynamics` model. The clock process uses the OCXO
   noise model from `ClockDynamics`, and relativistic clock-rate correction is
   applied with respect to the Moon for the ELFO case. Receiver app calls are scheduled at
   `receiver_app.rate_hz` in the receiver's local clock, then mapped to TDB for
   propagation and measurement generation. SRP can be enabled in truth and/or
   filter dynamics through a configurable ballistic coefficient.
4. Precompute link geometry and truth `GnssChannel` messages with light-time
   corrected transmitter states from PCO-corrected SP3 data. This stage writes
   every receiver-transmitter-frequency link with ionospheric/plasmaspheric
   delay columns set to zero.
5. Precompute ionospheric/plasmaspheric delays from the saved links with the
   GCPM ray tracer. This is a Python multiprocessing stage because the raytrace
   implementation is process-parallel rather than OpenMP/thread-parallel.
6. Generate noisy pseudorange, Doppler, and optional TDCP observations, then
   run a UDU EKF using `NBodyDynamics` for filter propagation. When TDCP is
   enabled, the filter state is `[x_k, x_{k-1}]` and the TDCP observable is
   modeled as the carrier-range difference between consecutive app ticks. The
   estimated measurement function uses the channel ephemeris coefficients and
   recomputes the light-time equation from the current estimated receiver
   state. Plasma can be left unmodeled in the filter and absorbed through
   measurement-noise inflation.
   Set `filter.estimate_srp_coefficient: true` to augment the filter state with
   the SRP ballistic coefficient `Cr * area / mass` in `m^2/kg`; otherwise the
   filter uses the fixed coefficient in `dynamics.srp_coeff_filter_m2_kg`. The
   SRP state-transition sensitivity is computed through
   `PropagateWithParams`, not finite differencing.
7. Write per-run trajectories, covariance-derived 3-sigma bounds, and a Monte
   Carlo summary under the configured output directory. The trajectory CSV
   includes both inertial and RTN position/velocity errors, RTN covariance
   bounds, and the number of tracked satellites at each receiver app tick.
8. Post-process the outputs in the notebook or the command-line plotting script
   to inspect plasma delays, tracked satellite counts, RTN state errors, and
   covariance envelopes.

Truth measurements are generated from the precomputed, light-time-corrected
SP3/ANTEX channels plus Sun-only Shapiro delay and optional precomputed
GCPM/IRI delay. TDCP, when enabled, is formed from time-differenced carrier
range in meters for channels tracked at consecutive receiver app ticks.
Estimated measurements reuse the ephemeris metadata carried in each channel and
re-solve light time from the current filter state. This keeps the truth data
product and receiver-side estimated measurement calculation separate.

The default ELFO is copied from `projects/Ephemeris/src/orbit_manager.py`:
`a=6541.4 km`, `e=0.6`, `i=65.5 deg`, `RAAN=60 deg`, `omega=90 deg`,
`M0=0 deg`, initialized in `MOON_OP` and converted to `MOON_CI` at the
configured UTC epoch.

Receiver propagation and app scheduling use a TDB coordinate timeline. GNSS
ephemeris tables are stored and interpolated in TAI seconds. Receiver states
use SI position and velocity units. Clock bias and drift are configured in
seconds in this example, while the library also supports meter and kilometer
clock-bias units for filters that benefit from range-like scaling.

## Run the staged pipeline

The single wrapper command is:

```bash
pixi run run-gnss-pipeline
```

To skip the GCPM delay batch when a delay table already exists, or when the
config has `plasma.simulate_truth: false`:

```bash
pixi run run-gnss-pipeline --skip-delays
```

The wrapper accepts `--config`, `--workers`, `--serial-delays`,
`--overwrite-delays`, `--skip-delays`, and `--no-plot`.

The same stages can be run manually:

```bash
pixi run precompute-gnss-links
pixi run precompute-gnss-delays
pixi run run-gnss-monte-carlo
pixi run plot-gnss-filtering
```

To populate SP3 files for new epochs from NASA CDDIS, configure Earthdata
Login credentials first, then run:

```bash
pixi run download-sp3-example
pixi run download-sp3-example-cpp
```

The preferred credential setup is `~/.netrc`:

```text
machine urs.earthdata.nasa.gov
  login YOUR_EARTHDATA_USERNAME
  password YOUR_EARTHDATA_PASSWORD
```

Set `chmod 0600 ~/.netrc` and never commit this file. For short local tests,
the Python and C++ loaders also accept `EARTHDATA_USERNAME` and
`EARTHDATA_PASSWORD`.
More detail is in `docs/tutorial/Python/gnss_files.rst`.

`precompute-gnss-links` builds the constellation from SP3/ANTEX, propagates the
nominal ELFO receiver trajectory, solves the light-time geometry, and writes:

```text
output/gnss_filtering/precomputed_links.csv
```

`precompute-gnss-delays` reads that file, runs GCPM ray tracing, and writes:

```text
output/gnss_filtering/precomputed_delays.csv
```

The delay precompute can also be run directly when you want to control the
number of worker processes:

```bash
python projects/GNSS_Filtering/precompute_delays.py \
  projects/GNSS_Filtering/gnss_filtering_config.yaml \
  --workers 8 \
  --overwrite
```

Use `--serial` for debugging one ray at a time. If `plasma.simulate_truth:
true`, the C++ Monte Carlo stage requires the configured delay file and stops
with an instruction to run the Python delay stage if it is missing.

This convenience task runs the full default staged workflow:

```bash
pixi run run-gnss-filtering
```

To build without running:

```bash
pixi run build-gnss-filtering
```

The default config is:

```text
projects/GNSS_Filtering/gnss_filtering_config.yaml
```

Receiver/link designs live in:

```text
projects/GNSS_Filtering/gnss_designs.yaml
```

You can also run the C++ stages directly with another config path:

```bash
./build-gnss-filtering/examples/ex_lunar_gnss_odts --config path/to/config.yaml --precompute
./build-gnss-filtering/examples/ex_lunar_gnss_odts --config path/to/config.yaml --run
```

## Outputs

The default output directory is `output/gnss_filtering`.

Each Monte Carlo run writes:

```text
trajectory_mc<N>.csv
```

The aggregate file is:

```text
summary.csv
```

The staged link and delay products are:

```text
precomputed_links.csv
precomputed_delays.csv
```

To inspect the Monte Carlo output, open:

```text
projects/GNSS_Filtering/plot_gnss_filtering_results.ipynb
```

The notebook plots ionospheric/plasmaspheric delay, RTN position/velocity
errors, clock errors, SRP coefficient error when enabled, and tracked satellite
count. To regenerate PNGs from the command line:

```bash
pixi run plot-gnss-filtering
```

The default config is a 3-minute, 1 Hz run using all GPS satellites available
in the selected SP3 file. For a faster debug run, set
`constellation.use_all_gps: false` and list a small subset in
`constellation.gps_prns`. Set `constellation.include_galileo: true` to add
Galileo when the SP3/ANTEX inputs include it. Set `plasma.simulate_truth: true`
to make the Monte Carlo stage consume the precomputed delay file. SRP is
enabled as a fixed truth/filter force by default; set
`filter.estimate_srp_coefficient: true` to estimate the coefficient as an
additional state. Set `measurements.use_tdcp: true` to add TDCP rows; the first
epoch has only current-epoch measurements, and later epochs add one TDCP row
per continuously tracked channel.

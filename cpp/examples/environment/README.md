# Environment Examples

This folder contains C++ examples for the LuPNT environment models. Most of the
files exercise the plasma stack imported from [pecsim](https://github.com/kdricemt/pecsim), which integrates GCPM v2.4, IRI-2007/IRI-2020, IGRF14, Kp loading, and ray tracing for GNSS TEC delays.

## Setup

The easiest path is to use pixi, because it sets the data paths needed by the
plasma runtime:

```bash
pixi shell
pixi run build
```

Outside `pixi shell`, set these paths before running the examples:

```bash
export LUPNT_DATA_PATH=$PWD/data/LuPNT_data
export PECSIMPY_BASE_PATH=$PWD/data/LuPNT_data/plasma
```

The plasma examples read model data from `$PECSIMPY_BASE_PATH/data` and write
generated CSV files below `$PECSIMPY_BASE_PATH/output`.

## Build

Build all examples:

```bash
cmake --build build --target all_examples
```

Build a smaller plasma subset:

```bash
cmake --build build --target \
  ex_plasma_gcpm \
  ex_plasma_iri \
  ex_plasma_igrf \
  ex_plasma_kp \
  ex_plasma_raytrace
```

## Quick Smoke Tests

These examples run quickly and are good first checks after changing the plasma
interfaces:

```bash
./build/examples/ex_plasma_igrf
./build/examples/ex_plasma_iri
./build/examples/ex_plasma_kp
./build/examples/ex_plasma_gcpm
./build/examples/ex_plasma_neldermead
./build/examples/ex_plasma_gcpm_repeat 0
./build/examples/ex_plasma_gcpm_mid 4 0 2
./build/examples/ex_plasma_openmp 8 1
```

## Density Slice Examples

The slice examples produce CSV grids for plotting electron density. They are
heavier than the smoke tests because they evaluate GCPM over a 201 by 201 grid.

```bash
./build/examples/ex_plasma_equatorial_slice
./build/examples/ex_plasma_meridianal_slice
./build/examples/ex_plasma_field_aligned
```

The configurable variants read simple `key = value` files from
`cpp/examples/environment/config` unless a path is provided on the command line:

```bash
./build/examples/ex_plasma_equatorial_slice_custom
./build/examples/ex_plasma_meridianal_slice_custom
./build/examples/ex_plasma_meridianal_slice_custom cpp/examples/environment/config/meridianal_slice_custom.cfg
```

Set `run_fortran = true` in the config only when comparing with the original
Fortran GCPM wrapper. The C++ GCPM implementation is preferred for repeatable
simulation runs because the original wrapper can preserve internal static state
between calls.

## Ray Tracing Example

`ex_plasma_raytrace` creates candidate GNSS-to-receiver links, solves a
light-time transmit position, and runs the ray tracer to compute TEC and group
delay terms.

```bash
# Args:
#   sim_nums cutoff_RE freq_idx step_size_km kp use_moon
#   min_tangent_alt_km max_tangent_alt_km mainlobe_angle_deg
./build/examples/ex_plasma_raytrace \
  2 4.0 1 50.0 -1.0 1 0.0 20000 20.0
```

The example expects `gps_2025_01_01.txt` to be available under the plasma data
directory. If it cannot be found, check `PECSIMPY_BASE_PATH` first.

## Notes

- `ex_plasma_sim_r12` generates multiple equatorial and meridianal slices while
  sweeping Rz12/solar-activity values.
- `ex_plasma_openmp` is a progress-bar and OpenMP smoke test; pass small
  arguments such as `8 1` when running it in CI-like environments.
- `ex_solar_system` is a plotting example for SPICE body queries and lunar
  mantle data. It builds with the other examples, but it opens matplot windows,
  so it is best run interactively.
- The fixed-date slice examples intentionally compare C++ GCPM with the
  original Fortran wrapper. The custom examples default to C++ only.

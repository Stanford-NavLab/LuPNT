# Python Examples

These examples are small entry points for exercising the `pylupnt` bindings.
Run them from the repository root after building the Python extension:

```bash
pixi run build-py
pixi run python python/examples/ex_states.py
```

## Core Examples

- `ex_states.py`: converts classical orbital elements to Cartesian form using SI units.
- `ex_coord_sys.py`: converts one Cartesian state from GCRF to ITRF.
- `ex_frame_converter.py`: converts a batch of Cartesian states between frames.
- `ex_spice_interface.py`: demonstrates SPICE-backed top-level time, ephemeris, and frame calls.


## Data And Credentials

`ex_gnss/download_sp3_epochs.py` downloads CDDIS SP3 products for several UTC
epochs and verifies that the files parse. Configure Earthdata Login credentials
before running it. The recommended setup is a `~/.netrc` entry:

```text
machine urs.earthdata.nasa.gov
  login YOUR_EARTHDATA_USERNAME
  password YOUR_EARTHDATA_PASSWORD
```

Then run:

```bash
pixi run download-sp3-example
```

For a longer setup guide, see `docs/pages/sp3_download.rst`.

## Plotting Examples

`test.py` plots built-in antenna gain patterns and opens a Matplotlib window
when an interactive backend is available. To smoke-test it without a window:

```bash
MPLBACKEND=Agg pixi run python python/examples/test.py
```

The notebooks in this directory are interactive demonstrations. They were not
rewritten in this cleanup pass because several depend on generated data files,
plotting backends, or external GNSS products.

SP3 Downloads
===================================================================

LuPNT can download precise GNSS SP3 orbit products from NASA CDDIS and cache
them under ``output/gnss_files/sp3``. The Python and C++ loaders use the same
cache directory, so a file downloaded from either interface can be reused by the
other.

This page describes the setup required before downloading, then shows how to
run the included Python and C++ examples.

Before You Download
-------------------------------------------------------------------

CDDIS access may require NASA Earthdata Login authentication. Before running
the examples, do the following:

1. Create an Earthdata Login account:

   .. code-block:: text

      https://urs.earthdata.nasa.gov

2. Review NASA's Earthdata Login API and integration page:

   .. code-block:: text

      https://www.earthdata.nasa.gov/engage/open-data-services-software/earthdata-developer-portal/earthdata-login-api

3. Configure command-line credentials with ``~/.netrc``.

   The file is named ``.netrc`` and lives in your home directory. It is not
   ``./netrc`` in the repository.

   .. code-block:: bash

      touch ~/.netrc
      chmod 0600 ~/.netrc
      cat >> ~/.netrc <<'EOF'
      machine urs.earthdata.nasa.gov
        login YOUR_EARTHDATA_USERNAME
        password YOUR_EARTHDATA_PASSWORD
      EOF

   Keep this file private and never commit it. The ``0600`` permission matters:
   many HTTP tools refuse to use a world-readable ``.netrc`` file.

4. If this is your first scripted CDDIS access, sign in through a browser once
   and approve any Earthdata application-access prompt for CDDIS.

Optional Environment Variables
-------------------------------------------------------------------

For short local tests, both LuPNT loaders also accept credentials from
environment variables:

.. code-block:: bash

   export EARTHDATA_USERNAME="YOUR_EARTHDATA_USERNAME"
   export EARTHDATA_PASSWORD="YOUR_EARTHDATA_PASSWORD"

Prefer ``~/.netrc`` for normal use so credentials are not exposed in shell
history, process environments, or terminal logs.

Run The Examples
-------------------------------------------------------------------

The Python example downloads or reuses SP3 files for several UTC epochs:

.. code-block:: bash

   pixi run download-sp3-example

The C++ example exercises the same CDDIS filename/URL logic and cache:

.. code-block:: bash

   pixi run download-sp3-example-cpp

Both examples print the selected SP3 product, local cache path, satellite IDs,
and TAI epoch coverage parsed from the file.

Use From Python
-------------------------------------------------------------------

The Python API is the C++ loader exposed through pybind11 (``pnt.Sp3Loader``).
Its static ``download_file_for_epoch`` helper computes the CDDIS product name,
downloads and caches missing files (under ``output/gnss_files/sp3``), and returns
the local path; the constructor then parses one or more files.

.. code-block:: python

   import pylupnt as pnt

   # A TAI epoch whose daily product we want.
   t_tai = pnt.convert_time(
       pnt.gregorian_to_time(2025, 1, 1, 0, 0, 0), pnt.Time.TDB, pnt.Time.TAI
   )

   # Download + cache the COD MGEX SP3 (Earthdata Login required), then parse it.
   sp3_path = pnt.Sp3Loader.download_file_for_epoch(t_tai, pnt.Time.TAI)
   sp3 = pnt.Sp3Loader(sp3_path)

   print(sp3.get_satellites()[:8])
   rv_ecef, clock_bias_s = sp3.get_pos_vel_clock("G01", t_tai)  # ECEF [m, m/s], bias [s]

The broadcast-ephemeris counterpart is ``pnt.RinexNavLoader`` (with the same
``download_file_for_epoch`` helper), and antenna phase-center offsets come from
``pnt.AntexLoader``.

Use From C++
-------------------------------------------------------------------

The C++ loader provides static helpers for epoch-based download/cache and then
normal parsing:

.. code-block:: cpp

   #include <lupnt/lupnt.h>

   using namespace lupnt;

   Real epoch_utc = GregorianToTime(2025, 1, 1, 0, 0, 0.0);
   std::filesystem::path sp3_path =
       Sp3Loader::DownloadFileForEpoch(epoch_utc, Time::UTC);

   Sp3Loader loader(sp3_path);
   const auto& sats = loader.GetSatellites();

Troubleshooting
-------------------------------------------------------------------

If LuPNT reports that CDDIS returned an Earthdata Login HTML page, the request
reached the login service but did not receive the SP3 file. Check:

- ``~/.netrc`` exists in your home directory, not in the repository.
- The machine entry is exactly ``urs.earthdata.nasa.gov``.
- The file permissions are ``0600``.
- Your Earthdata username and password are correct.
- You have approved CDDIS access in a browser if prompted.

You can test the same credential path with curl:

.. code-block:: bash

   curl -L -n -O "https://cddis.nasa.gov/archive/gnss/products/<week>/<file>.SP3.gz"

The flags matter: ``-L`` follows Earthdata redirects and ``-n`` uses
``~/.netrc``.

Cache Location
-------------------------------------------------------------------

By default, SP3 files are stored in:

.. code-block:: text

   output/gnss_files/sp3

When running through pixi, LuPNT environment variables are configured from the
workspace. If you run binaries manually, make sure ``LUPNT_DATA_PATH`` and any
custom output path are set consistently with the rest of your workflow.

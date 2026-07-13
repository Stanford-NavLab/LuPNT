.. _cpp_ex5_gnss_measurement_sim:

Example 5: GNSS Measurement Simulation
======================================

Simulate which Earth-GNSS signals a lunar receiver can track: propagate an
Elliptical Lunar Frozen Orbit (ELFO) receiver, load the GPS + Galileo
constellations from SP3 precise ephemerides, and run visibility + link-budget
(C/N\ :sub:`0`) at every epoch.

Mirrors :doc:`../Python/ex5_gnss_measurement_sim`.

The C++ program propagates the lunar ELFO receiver with ``NBodyDynamics``, builds
the GPS L1 + Galileo E1 constellations with
``GnssConstellation::SetupSatelliteStatesFromFiles``, configures a
``GNSSMeasurements`` model (Earth + Moon occultation, ``moongpsr`` receiver
antenna), and ``Precompute``\ s visibility + C/N\ :sub:`0` over one orbital
period, printing the per-epoch visible-satellite counts and C/N\ :sub:`0`
statistics. The full receiver-side model (light-time iteration, relativistic and
plasma corrections, online vs. batch ``Precompute``) is documented on the
:doc:`gnss_measurements` page.

.. note::

   This example is data-driven: it reads SP3 precise-ephemeris products, which
   require a free NASA Earthdata login and a one-time download (see
   :doc:`../../pages/sp3_download`). ``pixi run build-examples`` only *compiles*
   the program; the data is loaded at run time, and the program guards missing
   inputs — printing an informative note and exiting cleanly when the SP3 files
   or Earthdata credentials are absent.

.. literalinclude:: ../../../cpp/examples/tutorials/ex5_gnss_measurement_sim.cc
   :language: cpp

.. code-block:: bash

   pixi run build-examples
   ./build-examples-pixi/ex5_gnss_measurement_sim

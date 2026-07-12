.. _cpp_gnss_tutorials:

GNSS: Interface, Measurements, and ODTS
=======================================

Three of the Python example notebooks form the GNSS chain — the constellation
interface, receiver-side measurement simulation, and lunar orbit determination.
Each has a standalone C++ counterpart under ``cpp/examples/tutorials`` (embedded
in full below), mirroring the other tutorials: the notebooks additionally plot,
while the console programs summarise the same results numerically.

These examples are driven by SP3 precise-ephemeris products, which require an
Earthdata login and a one-time download (see :doc:`../../pages/sp3_download`).
``pixi run build-examples`` only *compiles* the programs; the data is loaded at
run time. Each program is written to guard missing data gracefully — it prints
an informative note and exits cleanly when the SP3/BRDC/ANTEX inputs (or
Earthdata credentials) are absent — so it always builds and runs, and produces
its full analysis once the data is present.

Example 3 — GNSS Interface
--------------------------

Python: :doc:`../Python/ex3_gnss_interface` — GNSS constellation, satellite
antenna gain patterns (ANTEX/ACE), and a broadcast-vs-precise ephemeris
comparison (SP3 vs BRDC with phase-centre-offset corrections, RTN decomposition).

In C++ the constellation catalog is ``GnssConstellation`` (``lupnt/agents``),
antenna patterns are loaded by ``Antenna`` and phase-centre offsets by
``AntexLoader`` (``lupnt/interfaces``), and the precise / broadcast ephemerides
are parsed by ``Sp3Loader`` / ``RinexNavLoader``. The program prints the
transmit-antenna peak gains, the ANTEX phase-centre offsets, and a per-PRN
radial/along-track/cross-track + clock error summary of BRDC against
(SP3 + PCO).

.. literalinclude:: ../../../cpp/examples/tutorials/ex3_gnss_interface.cc
   :language: cpp

.. code-block:: bash

   pixi run build-examples
   ./build-examples-pixi/ex3_gnss_interface

Example 5 — GNSS Measurement Simulation
---------------------------------------

Python: :doc:`../Python/ex5_gnss_measurement_sim` — propagate an ELFO, load the
GPS + Galileo constellations from SP3, and run visibility + link-budget (C/N0)
at every epoch.

The C++ program propagates the lunar ELFO receiver with N-body dynamics, builds
the GPS L1 + Galileo E1 constellations with
``GnssConstellation::SetupSatelliteStatesFromFiles``, configures a
``GNSSMeasurements`` model (Earth + Moon occultation, moongpsr receiver antenna),
and ``Precompute``\s visibility + C/N0 over one orbital period, printing the
per-epoch visible-satellite counts and C/N0 statistics. The full receiver-side
model (light-time iteration, relativistic and plasma corrections, online vs.
batch ``Precompute``) is documented on the :doc:`gnss_measurements` page.

.. literalinclude:: ../../../cpp/examples/tutorials/ex5_gnss_measurement_sim.cc
   :language: cpp

.. code-block:: bash

   pixi run build-examples
   ./build-examples-pixi/ex5_gnss_measurement_sim

Example 6 — GNSS Sidelobe ODTS (EKF)
------------------------------------

Python: :doc:`../Python/ex6_gnss_odts` — sidelobe pseudorange + TDCP orbit
determination and timing with an EKF, including a GCPM/IRI plasmaspheric-delay
ray-trace stage.

The C++ program builds a ``Simulation`` from ``configs/lunar_gnss_odts.yaml`` and
runs it: a physical ``receiver`` Spacecraft hosts the ``LunarGnssOdtsApp`` (the
release-facing ODTS engine under ``cpp/lupnt/applications/lunar_gnss_odts``),
whose scheduled step runs the whole Monte-Carlo EKF. The program then reads the
per-seed final / RMS position, velocity, and clock errors off the app. A full run
first needs the SP3 products and the staged link + plasma-delay precompute (see
``python/examples/ex6_precompute.py``):

.. literalinclude:: ../../../cpp/examples/tutorials/ex6_gnss_odts.cc
   :language: cpp

.. code-block:: bash

   pixi run python python/examples/ex6_precompute.py   # Stage 1 links + Stage 2 plasma delays
   pixi run build-examples
   ./build-examples-pixi/ex6_gnss_odts

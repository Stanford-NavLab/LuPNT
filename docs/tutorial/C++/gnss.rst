.. _cpp_gnss_tutorials:

GNSS: Interface, Measurements, and ODTS
=======================================

Three of the Python example notebooks form the GNSS chain — the constellation
interface, receiver-side measurement simulation, and lunar orbit determination.
Because these are driven by SP3 precise-ephemeris products (which require an
Earthdata login and a one-time download, see :doc:`../../pages/sp3_download`),
their C++ counterparts are the maintained, data-driven components below rather
than standalone single-file examples.

Example 3 — GNSS Interface
--------------------------

Python: :doc:`../Python/ex3_gnss_interface` — GNSS constellation visualisation,
satellite antenna gain patterns (ANTEX/ACE), and a broadcast-vs-precise ephemeris
comparison (SP3 vs BRDC with phase-centre-offset corrections, RTN decomposition).

In C++ the constellation catalog is ``GnssConstellation`` (``lupnt/agents``),
antenna patterns are loaded by ``AntexLoader`` (``lupnt/interfaces``), and the
broadcast ephemeris/almanac models are ``CartesianEphemeris`` / ``Almanac``
(``lupnt/applications``; see also :doc:`ex9_ephemeris`).

Example 5 — GNSS Measurement Simulation
---------------------------------------

Python: :doc:`../Python/ex5_gnss_measurement_sim` — propagate an ELFO, load the
GPS constellation from SP3, and run visibility + link-budget (C/N0) at every
epoch.

The full receiver-side C++ flow — light-time iteration, visibility, relativistic
and plasma corrections, and the online vs. batch ``Precompute`` paths of
``GNSSMeasurements`` — is documented in detail on the
:doc:`gnss_measurements` page.

Example 6 — GNSS Sidelobe ODTS (EKF)
------------------------------------

Python: :doc:`../Python/ex6_gnss_odts` — sidelobe pseudorange + TDCP orbit
determination and timing with an EKF, including a GCPM/IRI plasmaspheric-delay
ray-trace stage.

The release-facing C++ implementation is the ``LunarGnssOdtsApp`` under
``cpp/lupnt/applications/lunar_gnss_odts``, driven from a YAML config by
``pnt.Simulation`` (see ``python/examples/ex6_run_gnss_odts.py``) and the staged
pipeline it runs:

.. code-block:: bash

   pixi run run-gnss-pipeline          # precompute links, GCPM delays, run the filter

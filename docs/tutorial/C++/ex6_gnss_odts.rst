.. _cpp_ex6_gnss_odts:

Example 6: GNSS Sidelobe ODTS (EKF)
===================================

Sidelobe pseudorange + Doppler + TDCP orbit determination and time
synchronization (ODTS) with a stochastic-cloning UDU EKF, including a GCPM / IRI
plasmaspheric-delay ray-trace stage. It estimates position, velocity, clock
bias / drift, and an SRP coefficient from the weak Earth-GNSS sidelobe signals a
lunar receiver can track.

Mirrors :doc:`../Python/ex6_gnss_odts`.

The C++ program builds a ``Simulation`` from ``configs/lunar_gnss_odts.yaml`` and
runs it: a physical ``receiver`` Spacecraft hosts the ``LunarGnssOdtsApp`` (the
release-facing ODTS engine under ``cpp/lupnt/applications/lunar_gnss_odts``),
whose scheduled step runs the whole Monte-Carlo EKF. The program then reads the
per-seed final / RMS position, velocity, and clock errors off the app.

.. note::

   A full run first needs the SP3 precise-ephemeris products and the staged link
   + plasma-delay precompute produced by ``python/examples/ex6_precompute.py``
   (Stage 1 link geometry, Stage 2 GCPM/IRI plasma delays). As with the other
   GNSS examples the program guards missing data and exits cleanly if it is
   absent. See :doc:`../../pages/sp3_download` for the Earthdata / SP3 setup.

.. literalinclude:: ../../../cpp/examples/tutorials/ex6_gnss_odts.cc
   :language: cpp

.. code-block:: bash

   pixi run python python/examples/ex6_precompute.py   # Stage 1 links + Stage 2 plasma delays
   pixi run build-examples
   ./build-examples-pixi/ex6_gnss_odts

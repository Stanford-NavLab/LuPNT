.. _cpp_ex9_ephemeris:

Example 9: Ephemeris and Almanac Fitting
========================================

Study the accuracy-vs-broadcast-datasize trade-off of two lunar-satellite orbit
models fit to a numerically propagated truth trajectory: a precise, short-validity
``CartesianEphemeris`` (Chebyshev) and a coarse, long-validity ``Almanac``
(Algorithm 2 of Iiyama & Gao). ``EphemerisSimulation`` fits both over a sweep of
fitting-window lengths, reporting RMS / 95th-percentile RTN position and velocity
error and the minimum broadcast bit budget for a target position accuracy.

Mirrors :doc:`../Python/ex9_ephemeris`.

.. literalinclude:: ../../../cpp/examples/tutorials/ex9_ephemeris.cc
   :language: cpp

.. code-block:: bash

   pixi run build-examples
   ./build-examples-pixi/ex9_ephemeris

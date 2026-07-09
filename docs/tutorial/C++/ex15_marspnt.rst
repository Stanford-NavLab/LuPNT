.. _cpp_ex15_marspnt:

Example 15: Mars PNT
====================

A 9/3/1 Walker-delta navigation constellation at ~15,000 km altitude around Mars,
propagated for two Martian sidereal days with LuPNT's N-body dynamics (Mars 8x8
spherical-harmonic gravity field). The example decomposes the acceleration budget
on a representative satellite, builds the constellation in a Mars-equatorial frame
mapped to ``MARS_CI``, rotates the propagated states into the Mars body-fixed
frame, evaluates PDOP and visibility over a lat/lon grid, and computes a
single-point positioning fix at Jezero Crater (the Perseverance site). Where the
Python notebook renders the 3-D constellation, ground tracks and PDOP maps, this
program prints the equivalent summary statistics.

Mirrors :doc:`../Python/ex15_marspnt`.

.. literalinclude:: ../../../cpp/examples/tutorials/ex15_marspnt.cc
   :language: cpp

.. code-block:: bash

   pixi run build-examples
   ./build-examples-pixi/ex15_marspnt

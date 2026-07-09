.. _cpp_ex16_leopnt:

Example 16: LEO PNT
===================

The Earth companion to the Mars tutorial (Example 15). A dense 110/10/1 Walker
constellation at 600 km and 55° inclination — the kind of low-Earth-orbit fleet
proposed to augment GPS with strong, close-range signals — is propagated with an
EGM96 gravity field **and** Harris-Priester atmospheric drag (via
``NBodyDynamics::SetDragCoeff``). The example shows the acceleration budget with
drag, quantifies how far an unmodeled-drag ephemeris drifts over five days, maps
coverage and PDOP over a lat/lon grid, computes a single-point positioning fix in
San Francisco, and finally runs an on-board orbit- and time-determination (ODTS)
EKF that pulls the satellite's own state in from GPS pseudoranges. The Python
notebook renders the constellation, ground tracks, PDOP maps and filter
convergence; this program prints the equivalent summary statistics.

Mirrors :doc:`../Python/ex16_leopnt`.

.. literalinclude:: ../../../cpp/examples/tutorials/ex16_leopnt.cc
   :language: cpp

.. code-block:: bash

   pixi run build-examples
   ./build-examples-pixi/ex16_leopnt

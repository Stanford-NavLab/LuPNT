.. _cpp_ex8_isl_odts:

Example 8: Inter-Satellite-Link + GNSS ODTS
===========================================

A five-satellite lunar relay/navigation constellation (NASA LCRNS Reference
Constellation 3.1) where one hub satellite maintains a two-way crosslink to each
of the other four. The ``IslOdtsCoordinatorApp`` propagates every satellite's truth
trajectory, simulates two-way range/Doppler crosslink measurements, and runs an
onboard Schmidt-EKF on the hub — estimating its own ``[r, v, clock_bias,
clock_drift]`` while carrying each neighbour as a consider state. Optional
Earth-GPS pseudorange aiding adds absolute position and clock observability.

Mirrors :doc:`../Python/ex8_isl_odts`.

.. literalinclude:: ../../../cpp/examples/tutorials/ex8_isl_odts.cc
   :language: cpp

.. code-block:: bash

   pixi run build-examples
   ./build-examples-pixi/ex8_isl_odts

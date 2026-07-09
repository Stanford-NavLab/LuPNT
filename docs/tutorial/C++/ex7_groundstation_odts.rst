.. _cpp_ex7_groundstation_odts:

Example 7: Ground-Station Orbit Determination
=============================================

A lunar satellite is tracked by the three 70 m Deep Space Network antennas
(Goldstone, Canberra, Madrid) via two-way range and range-rate.
``GroundStationOdtsSimulation`` propagates the truth trajectory, runs an
elevation-mask visibility analysis, simulates noisy measurements over the
visible passes, and recovers the orbit from a perturbed initial guess with an
iterative batch (weighted least-squares) filter whose design matrix is built
analytically from the autodiff state-transition matrix.

Mirrors :doc:`../Python/ex7_groundstation_odts`.

.. literalinclude:: ../../../cpp/examples/tutorials/ex7_groundstation_odts.cc
   :language: cpp

.. code-block:: bash

   pixi run build-examples
   ./build-examples-pixi/ex7_groundstation_odts

.. _cpp_ex10_surface_rover:

Example 10: Surface Rover Navigation
====================================

A lunar-surface rover near the south pole fuses a full strapdown IMU
(accelerometer + gyroscope, Kalibr noise model), LCRNS/LANS pseudoranges from a
five-satellite relay constellation, and a LOLA-DEM altitude constraint, with the
IMU biases estimated online. ``RunSurfaceNav`` loads the DEM for the site
(downloading from NASA PGDA if not cached), synthesizes the measurements, runs
the strapdown-INS EKF, and returns the logged truth/estimate/covariance series.
The example re-runs with the DEM constraint disabled to show how the terrain
height ties down the weakly-observable vertical channel.

Mirrors :doc:`../Python/ex10_surface_rover`.

.. literalinclude:: ../../../cpp/examples/tutorials/ex10_surface_rover.cc
   :language: cpp

.. code-block:: bash

   pixi run build-examples
   ./build-examples-pixi/ex10_surface_rover

.. _cpp_ex10_surface_rover:

Example 10: Surface Rover Navigation
====================================

A lunar-surface rover near the south pole fuses a full strapdown IMU
(accelerometer + gyroscope, Kalibr noise model), LCRNS/LANS pseudoranges from a
five-satellite relay constellation, and a LOLA-DEM altitude constraint, with the
IMU biases estimated online. The scenario is agent-based
(``configs/surface_rover_nav.yaml``): a thin ``Rover`` agent hosts a
``SurfaceRoverNavApp`` (the strapdown-INS EKF), and a ``LunarNavConstellation``
agent stands in for the five navigation satellites — one config entry expanding
into five self-propagating ``Spacecraft`` — which the rover queries for the relay
geometry instead of propagating the orbits itself. The shared ``world:`` block
provides the LOLA DEM for the site (downloaded from NASA PGDA if not cached). The
example re-runs with the DEM constraint disabled to show how the terrain height
ties down the weakly-observable vertical channel.

Mirrors :doc:`../Python/ex10_surface_rover`.

.. literalinclude:: ../../../cpp/examples/tutorials/ex10_surface_rover.cc
   :language: cpp

.. code-block:: bash

   pixi run build-examples
   ./build-examples-pixi/ex10_surface_rover

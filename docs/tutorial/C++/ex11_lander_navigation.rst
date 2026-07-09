.. _cpp_ex11_lander_navigation:

Example 11: Lunar Lander Navigation
===================================

A lunar lander on powered descent to a south-pole site navigates with an
error-state INS EKF that fuses a full IMU (Kalibr noise model), a nadir radar
altimeter (height above DEM terrain), crater-landmark bearings (terrain-relative
navigation against a synthetic crater map), and LunaNet/LANS pseudoranges, with
the IMU biases estimated online. ``RunLanderNav`` loads the DEM, builds the
crater map, generates the powered-descent truth trajectory, and drives the
descent EKF. The example ends with a sensor-ablation study — disabling craters,
altimeter, and LunaNet in turn — to show what each aiding source buys.

Mirrors :doc:`../Python/ex11_lander_navigation`.

.. literalinclude:: ../../../cpp/examples/tutorials/ex11_lander_navigation.cc
   :language: cpp

.. code-block:: bash

   pixi run build-examples
   ./build-examples-pixi/ex11_lander_navigation

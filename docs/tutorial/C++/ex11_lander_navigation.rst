.. _cpp_ex11_lander_navigation:

Example 11: Lunar Lander Navigation
===================================

A lunar lander on powered descent to a south-pole site navigates with an
error-state INS MEKF that fuses a full IMU (Kalibr noise model), a nadir radar
altimeter (height above DEM terrain), crater-landmark bearings (terrain-relative
navigation against a synthetic crater map), and LunaNet/LANS pseudoranges, with
the IMU biases estimated online. The scenario is agent-based
(``configs/lander_nav.yaml``): the one ``Lander`` agent hosts **two applications**
— a guidance ``LanderGncApp`` that owns the powered-descent truth trajectory
(built-in smoothstep or an injected reference) and writes the lander's truth state
each epoch, and a navigation ``LanderNavApp`` that reads that truth to synthesize
its measurements and drive the descent MEKF. The shared ``world:`` block provides
the DEM. The example ends with a sensor-ablation study — disabling craters,
altimeter, and LunaNet in turn — to show what each aiding source buys.

Mirrors :doc:`../Python/ex11_lander_navigation`.

.. literalinclude:: ../../../cpp/examples/tutorials/ex11_lander_navigation.cc
   :language: cpp

.. code-block:: bash

   pixi run build-examples
   ./build-examples-pixi/ex11_lander_navigation

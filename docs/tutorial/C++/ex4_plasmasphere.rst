.. _cpp_ex4_plasmasphere:

Example 4: Plasmasphere and Ionosphere
======================================

Sample the GCPM v2.4 electron-density model along the equatorial radial
direction, then ray-trace a GNSS signal through the ionosphere/plasmasphere and
report the total electron content (TEC) and dispersive delay. This exercises the
``pecsim`` backend that sits behind ``pylupnt.plasma`` (GCPM v2.4, IRI-2007, and
the Haselgrove ray-tracer).

Mirrors :doc:`../Python/ex4_plasmasphere`.

.. literalinclude:: ../../../cpp/examples/tutorials/ex4_plasmasphere.cc
   :language: cpp

.. code-block:: bash

   pixi run build-examples
   ./build-examples-pixi/ex4_plasmasphere

The higher-fidelity, closed-loop ray-trace example (with a full GNSS
constellation and occultation gating) is ``cpp/examples/environment/ex_plasma_raytrace.cc``.

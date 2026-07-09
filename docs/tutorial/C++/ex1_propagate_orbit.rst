.. _cpp_ex1_propagate_orbit:

Example 1: Propagating an ELFO Orbit
====================================

Initialise a lunar Elliptical Frozen Orbit (ELFO) from classical orbital
elements in the Moon-centred Orbital Plane (OP) frame, convert the state to the
Moon-Centred Inertial (MCI) frame, propagate it for 30 days with full N-body
dynamics and Solar Radiation Pressure, and inspect the long-term libration of the
osculating orbital elements. The OP frame tracks the Earth's instantaneous orbit
plane about the Moon, the natural frame for the frozen-orbit parameters of Ely
(2005).

Mirrors :doc:`../Python/ex1_propagate_orbit`.

.. literalinclude:: ../../../cpp/examples/tutorials/ex1_propagate_orbit.cc
   :language: cpp

.. code-block:: bash

   pixi run build-examples
   ./build-examples-pixi/ex1_propagate_orbit

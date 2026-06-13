.. _introduction:

Introducing LuPNT
############################

``LuPNT`` is an open-source C++/Python library for Lunar Positioning,
Navigation, and Timing (PNT) research. It provides astrodynamics, signal
propagation models, measurement models, and navigation algorithms tailored for
cislunar missions.

The current codebase is organized around a modern C++ core with Python
bindings:

* ``cpp/lupnt`` contains the C++ library source.
* ``python/pylupnt`` contains the Python package and high-level utilities.
* ``configs`` contains reusable simulation and application configurations.
* ``projects`` contains research workflows, examples, and notebooks.
* ``docs`` contains the Sphinx, Breathe, and Exhale documentation source.

For local development, install the Pixi environment and build the project from
the repository root:

.. code-block:: bash

   pixi install
   pixi run build
   pixi run build-py

See :doc:`development` for the full setup and workflow guide.

C++
===================================================================

These tutorials are the standalone example programs under
``cpp/examples/tutorials``. Each mirrors the correspondingly numbered Python
example notebook, computing the same results (the notebooks additionally render
3-D plots and maps, which these console programs summarise as statistics). They
are built by the ``lupnt_examples`` CMake project — each ``.cc`` becomes its own
executable. Build them all with:

.. code-block:: bash

   pixi run build-examples          # configures cpp/examples and builds all_examples
   ./build-examples-pixi/ex1_propagate_orbit

Each page below embeds the full source of one example. The GNSS trio (interface,
measurement simulation, and orbit determination) is data-driven (SP3 precise
ephemerides) and is documented separately under :doc:`GNSS <gnss>`. The
constellation-design (ex12), Cesium (ex13) and optical-navigation (ex14)
notebooks are Python-only — they build on the optimiser, Cesium and Blender
tooling that lives in the Python layer — so they have no C++ counterpart here.

.. toctree::
    :maxdepth: 1
    :caption: Examples

    ex1_propagate_orbit
    ex2_time_conversions
    ex4_plasmasphere
    ex7_groundstation_odts
    ex8_isl_odts
    ex9_ephemeris
    ex10_surface_rover
    ex11_lander_navigation
    ex15_marspnt
    ex16_leopnt

.. toctree::
    :maxdepth: 1
    :caption: GNSS

    gnss
    gnss_measurements

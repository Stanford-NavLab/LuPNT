.. _cpp_ex2_time_conversions:

Example 2: Relativistic Time Conversions
========================================

Reproduce the analytical secular drift rates between the IAU coordinate time
scales relevant to cislunar operations — Terrestrial Time (TT), Selenocentric
Coordinate Time (TCL), and Lunar Time (LT) — using the defining scale factors
(``L_G``, ``L_B``, ``L_L``, ``L_EM``) in ``lupnt/core/constants.h``. The TCG–TT
and LT–TCL relations are exact-linear; the LT–TT rate adds the mean tidal
correction of Fienga et al. (2023).

Mirrors :doc:`../Python/ex2_time_conversions`.

.. literalinclude:: ../../../cpp/examples/tutorials/ex2_time_conversions.cc
   :language: cpp

.. code-block:: bash

   pixi run build-examples
   ./build-examples-pixi/ex2_time_conversions

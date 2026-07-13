.. _cpp_ex7_groundstation_odts:

Example 7: Ground-Station Orbit Determination
=============================================

A lunar satellite is tracked by the three 70 m Deep Space Network antennas
(Goldstone, Canberra, Madrid) via two-way range and range-rate. The whole
scenario lives in ``configs/ground_station_odts.yaml`` and is run by the generic
agent-based engine: each antenna is a ``GroundStation`` agent running a
``GroundStationTrackingApp`` (a sensor that generates its own visibility-gated,
noisy observations of the target), and a ``GroundStationManager`` agent runs the
``GroundStationManagerApp`` that aggregates those observations and recovers the
orbit from a perturbed initial guess with an iterative batch (weighted
least-squares) filter plus a square-root information filter / smoother, whose
design matrix is built analytically from the autodiff state-transition matrix.
By default the truth target and the estimator share the one force model in the
``world:`` block, so the driver just builds the ``Simulation`` from the YAML and
calls ``Run()``. The manager app also accepts an optional ``filter_dynamics:``
block (or per-filter ``batch_dynamics:`` / ``sequential_dynamics:``) to run the
batch and/or sequential filter at a deliberately different fidelity than truth.

Each ``GroundStationTrackingApp`` can optionally add Earth signal-path delays
(troposphere, ionosphere, relativistic Shapiro) and a solid Earth tide station
displacement to its truth observables via the ``apply_*`` keys in the config;
they default off, and an estimator that models a pure geometric range sees them
as realistic tracking errors.  The closed forms and their magnitudes are given
in :doc:`../../math/measurements` (Ground-Station Signal-Path and
Station-Location Corrections).

Mirrors :doc:`../Python/ex7_groundstation_odts`.

.. literalinclude:: ../../../cpp/examples/tutorials/ex7_groundstation_odts.cc
   :language: cpp

.. code-block:: bash

   pixi run build-examples
   ./build-examples-pixi/ex7_groundstation_odts

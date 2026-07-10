.. _cpp_ex8_isl_odts:

Example 8: Distributed Inter-Satellite-Link ODTS
================================================

A five-satellite lunar relay/navigation constellation (NASA LCRNS Reference
Constellation 3.1), fully cross-linked, plus a small southern surface-station
beacon network. The scenario is built from PHYSICAL agents in
``configs/isl_odts_distributed.yaml``: each satellite is a ``Spacecraft`` (self-
propagating orbit+clock truth) hosting a ``SatelliteOdtsApp`` onboard Schmidt-EKF
— estimating its own ``[r, v, clock_bias, clock_drift]`` while carrying each
neighbour as a consider state, generating its own crosslink + station-aiding
measurements, and exchanging consider-state with its neighbours. Each
``SurfaceStation`` hosts a ``StationBeaconSensor`` feeding a
``SurfaceStationManager``'s ``GroundOdtsApp`` — the centralized ground filter
(no ISL) — for a head-to-head comparison. ``Simulation`` runs the event loop;
the per-epoch truth/estimate series are read off each app's accessors.

Mirrors :doc:`../Python/ex8_isl_odts`.

.. literalinclude:: ../../../cpp/examples/tutorials/ex8_isl_odts.cc
   :language: cpp

.. code-block:: bash

   pixi run build-examples
   ./build-examples-pixi/ex8_isl_odts

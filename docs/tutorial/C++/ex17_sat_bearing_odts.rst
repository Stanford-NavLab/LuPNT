.. _cpp_ex17_sat_bearing_odts:

Example 17: Satellite-to-Satellite Angles-Only OD
=================================================

The smallest agent-based orbit-determination scenario, and a template for authoring your
own. Two ``Spacecraft`` share the lunar ``world:`` force model. An **observer** measures only
the **bearing** (unit line-of-sight direction) to a **target** whose ephemeris is treated as
**known**, and an EKF on the observer estimates the observer's own 6-state orbit
:math:`[\mathbf r, \mathbf v]` (no clock). Measuring a direction to a known landmark is a
well-posed, observable problem that converges (here to the tens-of-meters level over a 6 h
arc from a 300 m a-priori with 2 arcsec bearings); estimating an unknown target's range from
bearings alone is the ill-posed angles-only problem, deliberately not attempted.

Only two new C++ classes are added — the observable ``SatBearingMeasurement`` (a
``MeasurementClone`` returning the unit line-of-sight with an autodiff Jacobian) and the
estimator ``AnglesOdtsApp`` (``REGISTER_FACTORY_CLASS(Application, ...)``); the agents,
dynamics, EKF and Monte-Carlo are all reused. This program loads the config, runs the
scenario, and prints the observer's converged position error.

Mirrors :doc:`../Python/ex17_sat_bearing_odts`.

.. literalinclude:: ../../../cpp/examples/tutorials/ex17_sat_bearing_odts.cc
   :language: cpp

.. code-block:: bash

   pixi run build-examples
   ./build-examples-pixi/ex17_sat_bearing_odts

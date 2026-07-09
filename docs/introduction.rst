.. _introduction:

Introducing LuPNT
############################

``LuPNT`` is an open-source C++/Python library for Lunar Positioning,
Navigation, and Timing (PNT) research. It provides high-fidelity astrodynamics,
signal-propagation models, measurement models, and navigation algorithms
tailored for cislunar missions. LuPNT is a product of the
`Stanford NAV Lab <https://navlab.stanford.edu/>`_.

The library is organized as a **config-driven, agent-based simulation
framework**. *Agents* (satellites, ground stations, rovers, landers, surface
stations, and constellations) host *Applications* — the mission and
navigation-filter logic — and run together on a shared, event-scheduled
*Simulation*. Agents, their devices, dynamics, and applications are all
assembled from YAML configuration, so new scenarios are described in
configuration rather than code. See :doc:`How to Create a New Simulation
<pages/new_simulation>` for the architecture and the workflow for building a
scenario.

Architecture at a glance
========================

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Module
     - Description
   * - **Agents & Applications**
     - Config-driven agents (satellite, ground station, rover, lander, surface
       station, constellation) that host navigation / mission **Applications** —
       ODTS, ISL, surface-rover and lander navigation, ephemeris / almanac
       generation — built from YAML via an asset factory.
   * - **Simulations**
     - Event-scheduled ``Simulation`` engine plus end-to-end scenario drivers
       (ground-station and inter-satellite ODTS, lunar GNSS ODTS, lander /
       surface navigation, ephemeris fitting).
   * - **Devices & sensors**
     - Clocks, IMUs, cameras, GNSS receivers, and communication devices
       attached to agents.
   * - **Dynamics**
     - Two-body, N-body, and high-fidelity force models (gravity, drag, SRP) for
       Earth, lunar, and arbitrary central-body orbits (see :doc:`math/dynamics`).
   * - **Environment**
     - Gravity fields, atmosphere, solar-system bodies, occultation, and plasma /
       ionosphere models.
   * - **Plasma**
     - GCPM v2.4 plasmasphere + IRI ionosphere electron density; GNSS
       ray-tracing with TEC and signal delay for cislunar links (integrated
       ``pecsim``).
   * - **GNSS**
     - GPS/GNSS signal generation, SP3/ANTEX-backed constellations,
       yaw-steering / attitude models, and space-user ranging.
   * - **Measurements**
     - Pseudorange, Doppler, and carrier-phase measurement models with
       light-time, Shapiro, and optional plasma corrections
       (see :doc:`math/gnss_measurements`).
   * - **Filters & States**
     - Extended Kalman Filter (EKF), Unscented Kalman Filter (UKF), square-root
       information filter (SRIF) / smoother, and batch least-squares estimators
       over composable state / joint-state abstractions (see :doc:`math/filters`).
   * - **Interfaces**
     - Data loaders and I/O: SP3, ANTEX, RINEX nav, TLE, SPICE kernels,
       EOP/TAI-UTC, LOLA DEM and crater data, plus Cesium and Matplotlib export.
   * - **Conversions**
     - Reference-frame transformations (ECI, ECEF, LVLH, Moon-centered, generic
       body-fixed / inertial) and time-system utilities
       (see :doc:`math/frame_conversions`).
   * - **Numerics**
     - Numerical integration (RK4, RK8, RKF45), a Nelder-Mead optimizer, and
       matrix utilities (see :doc:`math/integration`).
   * - **Visualization**
     - Matplotlib / Plotly plotting plus interactive 3-D
       `CesiumJS <https://cesium.com/platform/cesiumjs/>`_ scenes of
       constellations and surface assets (``pnt.plot.CesiumScene``).
   * - **Python bindings**
     - Full ``pylupnt`` Python package exposing the C++ library via pybind11.

Repository layout
=================

* ``cpp/lupnt`` — the C++ library source, organized by module (``agents/``,
  ``applications/``, ``simulations/``, ``devices/``, ``dynamics/``,
  ``environment/``, ``measurements/``, ``states/``, ``interfaces/``,
  ``conversions/``, ``core/``, ``numerics/``).
* ``cpp/examples`` — C++ example programs, including the
  :doc:`tutorial counterparts <tutorial/C++/index>` to the Python notebooks.
* ``python/pylupnt`` — the Python package (installed in place) and high-level
  utilities (plotting, interfaces, plasma).
* ``python/examples`` — the :doc:`tutorial notebooks <tutorial/Python/index>`
  (ex1–ex16).
* ``configs`` — reusable YAML scenario configuration (agents, applications,
  dynamics, environments, datasets).
* ``projects`` — research workflows, notebooks, and scripts.
* ``data/LuPNT_data`` — runtime data (ephemeris, GNSS products, plasma
  coefficients, TLEs), fetched automatically on first build.
* ``docs`` — the Sphinx, Breathe, and Exhale documentation source.

Getting started
===============

Install the `Pixi <https://pixi.sh>`_ environment (it manages the compiler
toolchain, Python, and every C++ dependency via conda-forge) and build the
project from the repository root:

.. code-block:: bash

   pixi install
   pixi run build       # build the C++ library
   pixi run build-py    # build and deploy the Python bindings

The runtime dataset ``data/LuPNT_data`` (~650 MB) is downloaded automatically on
the first build. See :doc:`development` for the full setup, the VS Code
debugging workflow, and the Windows/WSL instructions, and :doc:`builddocs` for
building this documentation. The Tutorial section walks through the library
hands-on in both :doc:`Python <tutorial/Python/index>` and
:doc:`C++ <tutorial/C++/index>`.

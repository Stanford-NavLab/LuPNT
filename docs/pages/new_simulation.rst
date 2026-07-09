.. _new_simulation:

How to Create a New Simulation
==============================

LuPNT is a **config-driven, agent-based** framework: a scenario is a set of
*Agents* that own *Dynamics*, *Devices*, and a *State*, and that host
*Applications* (mission / navigation logic). Applications build *Measurements*
and run *Filters* over composable *States*. Everything is scheduled by an
event-driven ``Simulation`` engine and assembled from YAML
configuration through a string-keyed asset factory.

This page maps those building blocks and their relationships, then gives the
minimal recipe for authoring a new scenario.

Two kinds of "simulation"
-------------------------

It is worth distinguishing two entry points up front:

* **The generic event-driven engine** — ``lupnt::Simulation``
  (``cpp/lupnt/simulations/simulation.{h,cc}``), built entirely from a YAML
  config and run by the loader ``cpp/examples/core/ex_sim.cc``. Agents,
  applications, devices, and dynamics are assembled by name via factories. This
  is the *"describe a scenario in config and register your class"* path
  (**Path A** below).
* **Monolithic scenario drivers** — e.g. ``GroundStationOdtsSimulation``,
  ``IslOdtsSimulation``, ``LunarGnssODTSSimulation``, ``EphemerisSimulation``,
  and the ``RunLanderNav`` / surface-nav drivers (under
  ``cpp/lupnt/simulations/<name>/``). These subclass ``Simulation`` (or are free
  functions), take a typed C++ ``Config`` struct, and wire agents/apps directly
  in code. This is the path for tightly-coupled estimators, error-state INS, and
  distributed filters (**Path B** below).

Both share the same building blocks.

The building blocks
-------------------

.. list-table::
   :header-rows: 1
   :widths: 16 30 54

   * - Block
     - Base class (header)
     - Role
   * - **Simulation**
     - ``Simulation`` (``simulations/simulation.h``)
     - Discrete-event engine: a time-ordered priority queue of ``Event`` s,
       ``Schedule`` / ``Run`` / pub-sub, and the ``agents_`` / ``channels_`` /
       ``constellations_`` registries built from config.
   * - **Agent**
     - ``Agent`` → ``AgentWithDynamics`` (``agents/agent.h``)
     - A platform (satellite, ground station, rover, lander, surface station).
       Owns a ``Dynamics``, ``Devices``, a ``State``, and hosts **one**
       ``Application``. Answers ``GetStateAt(t)`` for measurement geometry.
   * - **Application**
     - ``Application`` (``applications/application.h``)
     - The mission / navigation logic run per step. Builds ``Measurement`` s and
       runs a ``Filter``. Compose several via the sub-app pattern
       (``LunaNetSatApp`` + ``LunaNetSubApp``).
   * - **Device**
     - ``Device`` (``devices/device.h``)
     - A sensor / component on an agent (``Clock``, ``Imu``, ``Camera``,
       ``Transmitter`` / ``Receiver`` / ``Transponder``, ``GnssReceiver``).
       Steps on its own schedule; read by apps via ``GetDevice(name)``.
   * - **Measurement**
     - ``Measurement`` / ``ErrorStateMeasurement`` (``measurements/measurement.h``)
     - The observable model. ``Compute(x, H)`` returns the predicted value,
       covariance, and Jacobian :math:`H = \partial h/\partial x`;
       ``CreateFunction()`` wires it into a ``Filter``.
   * - **Dynamics**
     - ``Dynamics`` (``dynamics/dynamics.h``)
     - The propagation model. ``Propagate(x0,t0,tf[,u][,stm])`` advances a
       ``State`` and optionally its STM :math:`\Phi = \partial x_f/\partial x_0`.
   * - **State**
     - ``State`` and ``JointState`` (``states/state.h``, ``states/joint_state.h``)
     - A labeled Eigen vector with a type, per-element names/units, and a frame.
       ``JointState`` composes orbit + clock (+ parameters) into one filter state.

Simulation — the engine
~~~~~~~~~~~~~~~~~~~~~~~~~

``Simulation`` is a discrete-event engine backed by a
``std::priority_queue<Event>``. An ``Event`` (``core/event.h``) is a
``{time, frequency, priority, callback}`` tuple ordered **earliest-time first,
higher-priority first**. Priorities resolve ties at equal time so that, within a
step, ``DYNAMICS`` / ``AGENT`` callbacks run before ``DEVICE`` callbacks, which
run before ``APPLICATION`` / ``LOGGING`` callbacks:

.. code-block:: text

   DYNAMICS = AGENT (2)  >  DEVICE (1)  >  APPLICATION = LOGGING (0)

The main loop (``Simulation::Run``) repeatedly pops the earliest event, sets the
simulation clock, invokes its callback, and — if the event has a non-zero
``frequency`` — reschedules it at ``t + 1/frequency`` (this is how periodic
``Step`` callbacks recur). ``Schedule(...)`` adds events; ``Subscribe`` /
``Publish`` provide a topic-based message bus. ``GetAgent(name)`` /
``GetChannel(name)`` accept bare or ``<simname>/<name>`` keys.

Agent — a platform
~~~~~~~~~~~~~~~~~~~

``Agent`` owns ``name_``, a back-pointer to the ``Simulation``, a ``State``, a
map of ``Devices``, and a single ``Application``. Its pure-virtual
``Cart6 GetStateAt(Real t) const`` is the contract every agent must answer:
*"where am I at time t?"* — used by measurement models for light-time geometry.

Almost every physical agent derives from ``AgentWithDynamics``, which adds a
``Ptr<Dynamics>``, an ``AttitudeDynamics``, and a ``Cart6 state_``. Its
``Step(t)`` calls ``Propagate(t)`` (integrating the dynamics in place) then
``Log(t)``. Concrete agents (all registered with the factory except ``Lander``):
``Satellite``, ``GroundStation``, ``Rover``, ``SurfaceStation``, ``Lander``.
``Constellation`` / ``GnssConstellation`` are group generators built from the
``constellations:`` config block.

Application — the mission logic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Application`` holds a back-pointer to its owning ``Agent``. ``Setup()``
schedules a periodic ``Step(t)`` (at ``APPLICATION`` priority) and seeds the
filter state/covariance; the pure-virtual ``Step(Real t)`` runs one
predict/update cycle. An agent hosts **at most one** application; compose several
pieces of logic with the **sub-app pattern**: ``LunaNetSatApp`` holds a
``std::vector<Ptr<LunaNetSubApp>>`` and fans its ``Step`` out to each sub-app
(e.g. ``IslOdtsApp`` and ``EphemerisGenApp``).

Concrete apps: ``GroundStationOdtsApp`` (the config-driven, factory-registered
batch-OD app), ``SurfaceStationApp``, ``LunaNetSatApp`` (sub-app host),
``IslOdtsApp`` / ``EphemerisGenApp`` (sub-apps), ``LanderNavApp`` and
``SurfaceRoverNavApp`` (driven directly by their C++ simulations).

Device, Measurement, Dynamics, State
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Device** — created from the agent's ``devices:`` block, stored in
  ``Agent::devices_``, and stepped independently on its own ``DEVICE``-priority
  events, so its state is current when an application reads it with
  ``GetDevice(name)``.
* **Measurement** — ``MeasData Compute(const State& x, MatXd* H) const`` returns
  the predicted value :math:`z = h(x)`, its covariance :math:`R`, and (when
  ``H != nullptr``) the Jacobian :math:`H = \partial h/\partial x`. Derive from
  ``MeasurementClone<Derived>`` (which supplies ``Clone``) or, for INS/error-state
  observables, from ``ErrorStateMeasurement``. ``CreateFunction()`` adapts a
  measurement to ``Filter::SetMeasurementFunction``. See :doc:`../math/filters`
  and :doc:`../math/gnss_measurements`.
* **Dynamics** — the ``Propagate`` family is the contract (see
  :doc:`../math/dynamics` and :doc:`../math/integration`). The STM overload is
  used by filters for covariance propagation; if a subclass does not provide one,
  the base computes it by autodiff. Concrete: ``NBodyDynamics``,
  ``CartesianTwoBodyDynamics``, ``JToCartTwoBodyDynamics``,
  ``JointOrbitClockDynamics``, ``ClockDynamics``, ``ImuDynamics``, the analytical
  relative-motion models, and ``SurfaceDynamics2D``.
* **State** — a labeled vector (``Cart6``, ``ClassicalOE``, ``ClockState3``,
  ``ImuState``, ...). ``JointState`` concatenates ``(State, Dynamics, ParamState,
  est_types, tauq)`` triples and exposes the combined dynamics/STM/process-noise
  functions the ``Filter`` consumes — this is how orbit + clock (+ estimated or
  consider parameters) become one filter state.

How they relate
---------------

Containment and ownership:

.. code-block:: text

   Simulation
     ├── priority_queue<Event>              (the clock / scheduler)
     ├── agents_         : map<name, Agent>
     ├── channels_       : map<name, Channel>
     └── constellations_ : map<name, Constellation>

   Agent (AgentWithDynamics)
     ├── Dynamics  dynamics_     ── owns ParamState (force / clock parameters)
     ├── AttitudeDynamics
     ├── State state_, Attitude attitude_, State control_
     ├── devices_ : map<name, Device>   (Clock, Imu, Camera, Tx/Rx, GnssReceiver)
     └── Application application_        (one; compose many via LunaNetSatApp + sub-apps)

   Application
     ├── Agent* agent_                   (back-pointer, GetAgent())
     ├── builds Measurement model(s) + a Filter (EKF / UKF / SRIF / batch)
     └── JointState composes State + Dynamics (+ params) for the Filter

**Per-step data flow.** At a given time, the engine runs due callbacks in
priority order: ``AgentWithDynamics::Step`` first advances the *truth* state via
``Propagate``; then ``Device`` steps advance clocks/IMUs/comms; then the
``Application::Step`` reads other agents' truth via ``GetStateAt(t)`` and its own
devices via ``GetDevice(...)``, forms a ``Measurement`` (``z, H, R``), and runs a
filter predict (using the ``Dynamics`` / ``JointState`` STM) and update. Results
are logged through ``Application::Log`` ← ``Agent::Log``. Periodic events
reschedule themselves.

Configuration and the asset factory
-----------------------------------

Classes are instantiated by name through
``AssetFactory<Base, Config&>`` (``core/asset_factory.h``). A class opts in with

.. code-block:: cpp

   REGISTER_FACTORY_CLASS(Application, MyApp)   // at file scope in MyApp.cc

which registers a creator under the string ``"MyApp"``; the YAML ``class:`` value
must match that string exactly. Typed aliases exist for each base
(``AgentFactory``, ``DeviceFactory``, ``ApplicationFactory``,
``DynamicsFactory``, ``ChannelFactory``).

``Simulation::Setup`` reads the top-level ``name`` / ``epoch`` / ``duration`` /
``log_level`` / ``channels``, then for each ``agents:`` entry uses its ``class:``
to build the agent, whose constructor in turn builds its ``dynamics:``,
``initial_state:``, ``devices:``, and ``application:`` sub-blocks through the
matching factories. For example (abridged from
``configs/ground_station_odts.yaml``):

.. code-block:: yaml

   agents:
     sat:
       class: Satellite
       dynamics: { class: NBodyDynamics, bodies: [...], integrator: RK4 }
       initial_state: { class: ClassicalOE, frame: MOON_OP, a: ..., e: ..., i: ... }
     DSS14:
       class: GroundStation
       latitude_deg: 35.4
       longitude_deg: -116.9
       altitude_m: 1001
       application:
         class: GroundStationOdtsApp
         target: sat                 # resolved via GetSimulation()->GetAgent("sat")
         elevation_mask_deg: 10
         use_range: true
         range_sigma_m: 1.0

Run it with the generic loader:

.. code-block:: bash

   ./build/examples/ex_sim configs/ground_station_odts.yaml

Authoring a new simulation
--------------------------

Path A — config-driven (preferred)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use this when the scenario fits the agent model. **Template:**
``configs/ground_station_odts.yaml`` +
``cpp/examples/simulations/ex_ground_station_odts_app.cc``.

#. **Reuse existing classes.** Write a YAML in ``configs/`` listing ``agents:``
   with already-registered ``class:`` values (``Satellite``, ``GroundStation``,
   ``Rover``, ``SurfaceStation``), a ``dynamics:`` block (e.g. ``NBodyDynamics``),
   ``devices:``, and an ``application:`` block. Run with ``ex_sim``.
#. **New Application** (the most common extension) — subclass ``Application``,
   implement ``MyApp(Config&)``, ``Setup()`` (resolve target agents via
   ``GetAgent``, seed the filter; the base ``Setup`` schedules the steps), and
   ``Step(Real)`` (build the ``Measurement``, run the filter). Add
   ``REGISTER_FACTORY_CLASS(Application, MyApp)``. Reference:
   ``applications/ground_station_odts_app.{h,cc}``.
#. **New Measurement** (only for a new observable) — subclass
   ``MeasurementClone<MyMeas>`` (or ``ErrorStateMeasurement``) and implement
   ``Compute(const State&, MatXd* H)``. Wire it into a filter with
   ``CreateFunction()``. Reference: ``measurements/crosslink_measurement.{h,cc}``.
#. **New Dynamics** (only for new physics) — subclass ``NumericalDynamics``
   (implement ``ComputeRates`` + ``Propagate`` + ``GetStateType``) or ``Dynamics``
   directly; ``REGISTER_FACTORY_CLASS(Dynamics, MyDyn)``. The STM is provided by
   autodiff if you do not supply one.
#. **New Agent** (only for a new platform kind) — subclass ``AgentWithDynamics``,
   implement ``GetStateAt``, optionally override ``Propagate`` / ``Setup``;
   ``REGISTER_FACTORY_CLASS(Agent, MyAgent)``. Reference: ``agents/satellite.cc``,
   ``agents/rover.cc``.

Path B — monolithic C++ driver
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use this for tightly-coupled estimators, error-state INS, or distributed
filters. **Template:** ``simulations/lander_nav/`` (``LanderNavApp``, driven
directly) or ``simulations/isl_odts/`` (``LunaNetSatApp`` + ``IslOdtsApp``
sub-app).

#. Subclass ``Simulation`` and override ``Setup()`` / ``Run()``, or write a
   ``RunMySim(config)`` free function like ``RunLanderNav``.
#. In C++, construct the agents (``std::make_shared<Lander>(cfg)``), attach the
   application (``agent.SetApplication(app)``, or ``sat_app->AddSubApp(app)`` for
   the LunaNet sub-app composition), and drive ``Predict`` / ``Update*``
   explicitly rather than through the scheduler.

Reference files
---------------

``simulations/simulation.{h,cc}``, ``core/event.h``, ``core/asset_factory.h``,
``agents/agent.h``, ``applications/application.h``,
``applications/lunanet_sat_app.h``, ``devices/device.h``,
``measurements/measurement.h``, ``dynamics/dynamics.h``, ``states/state.h`` +
``states/joint_state.h``, ``configs/ground_station_odts.yaml``,
``examples/core/ex_sim.cc``.

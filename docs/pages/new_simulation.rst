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

Every scenario is built from a YAML config by the same generic event-driven
engine, ``lupnt::Simulation`` (``cpp/lupnt/simulations/simulation.{h,cc}``).
Agents, applications, devices, and dynamics are assembled by name through
factories, and a top-level ``world:`` block defines the shared read-only
environment (epoch, integration frame, force model) that both the truth agents
and the estimators draw from. The two things that vary are *how much logic lives
in config vs. code* and *how you launch the config*:

* **Path A — config-driven (preferred).** Describe the scenario entirely in
  YAML with already-registered ``class:`` values, and put any new mission /
  navigation logic in an ``Application`` that self-drives from the ``world:``
  block. This is how the ground-station ODTS, GNSS ODTS, ISL ODTS, ephemeris,
  surface-rover and lander scenarios all work now (see ``configs/*.yaml``).
* **Path B — in-code driver.** For a tightly-coupled estimator you may still
  wire the agents/apps together in a small C++ or Python driver, but it is run
  by the *same* engine — there is no separate monolithic ``Simulation``
  subclass with a hand-written ``Run()`` loop anymore.

Launch a config from Python with ``pnt.Simulation(yaml_or_dict)`` (the pattern
used by every ``python/examples/exN_run_*.py`` script) or from a small C++
driver such as the tutorials under ``cpp/examples/tutorials/``. Both share the
same building blocks.

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
       ``constellations_`` registries built from config. Owns the ``World``.
   * - **World**
     - ``World`` (``simulations/world.h``)
     - The shared, **read-only** physical environment built from the top-level
       ``world:`` block: epoch, integration ``frame``, and a force model.
       ``MakeDynamics()`` hands out a fresh ``NBodyDynamics`` from that one force
       model (so truth and estimator share the same physics), and it exposes a
       point-mass ``Gravity()`` plus an optional DEM/terrain service for surface
       scenarios. It never propagates agents — it only *provides* the environment.
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

World — the shared environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Simulation::Setup`` builds a single ``World`` from the top-level ``world:``
block *before* any agents, and every agent can reach it with
``Agent::GetWorld()``. The World holds the epoch, the integration ``frame``, and
one force-model ``Config``; ``MakeDynamics()`` returns a fresh ``NBodyDynamics``
built from that force model (autodiff on, so the estimator gets an analytic STM).
Because the truth agent and the estimator both derive their dynamics from the
*same* World force model, there is no truth/filter model mismatch. For surface
scenarios the ``world:`` block also accepts a point-mass ``gravity:`` body and a
``dem:`` terrain block (``World::GetElevation``, ``EnuToWorld``,
``SiteCenterWorld``). The World is **read-only**: it provides the environment and
a truth facade, but each ``AgentWithDynamics`` still self-propagates via its own
``dynamics_``.

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

Concrete apps (all factory-registered): the **sensor / estimator split** for
ground-station OD — ``GroundStationTrackingApp`` on each ``GroundStation``
(a sensor: visibility-gated range / range-rate, pushed to the manager) feeding a
``GroundStationManagerApp`` on a ``GroundStationManager`` agent (the centralized
batch + SRIF/smoother estimator); ``LunarGnssOdtsApp`` (GNSS ODTS),
``IslOdtsCoordinatorApp`` (centralized ISL) and the distributed ``SatelliteOdtsApp``
(per-satellite onboard filter) + ``GroundOdtsApp``; ``EphemerisApp`` /
``LunaNetSatApp`` (sub-app host) + ``IslOdtsApp`` / ``EphemerisGenApp`` (sub-apps);
``SurfaceStationApp``; and the self-driving error-state INS apps
``SurfaceRoverNavApp`` and ``LanderNavApp`` (a thin ``Rover`` / ``Lander`` agent
hosts them and they precompute the truth trajectory and synthesize their own
measurements from the shared ``World`` each ``Step``).

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

As pseudocode, the engine loop (``Simulation::Run``) and a periodic navigation
``Application::Step`` are:

.. code-block:: text

   # Simulation::Run  (simulations/simulation.cc)
   while queue not empty and queue.top().time <= duration:
       e = queue.pop()                       # earliest time, then highest priority
       time = e.time
       e.callback(time)                      # Agent/Device/Application Step, or a Publish
       if e.frequency > 0:                   # periodic -> reschedule
           queue.push(Event(time + 1/e.frequency, e.frequency, e.priority, e.callback))
   DataLogger::Flush()

   # AgentWithDynamics::Step(t)   (priority DYNAMICS) — advance the truth state
   state_, attitude_ = dynamics_.Propagate(state_, time_, t, control_)
   time_ = t;  Log(t)

   # A navigation Application::Step(t)   (priority APPLICATION)
   for tx in visible_transmitters(t):        # geometry via tx.GetStateAt(t)
       z, H, R = measurement(tx, receiver=agent_).Compute(x_hat, &H)   # h(x), Jacobian, noise
       stack z, H, R
   filter.Predict(t)                         # x_hat, P via Dynamics/JointState STM + process noise
   filter.Update(z, H, R)                    # innovation, gain, covariance update
   Log(t)

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
``log_level`` / ``channels``, builds the shared ``World`` from the ``world:``
block, then for each ``agents:`` entry uses its ``class:`` to build the agent,
whose constructor in turn builds its ``dynamics:``, ``initial_state:``,
``devices:``, and ``application:`` sub-blocks through the matching factories. For
example (abridged from ``configs/ground_station_odts.yaml``):

.. code-block:: yaml

   world:                            # shared, read-only environment
     frame: MOON_CI
     force_model: { integrator: RKF45, autodiff: true,
                    bodies: [ { MOON: { n: 2, m: 2 } }, { EARTH: {} }, { SUN: {} } ] }

   agents:
     sat:                            # truth target
       class: Satellite
       dynamics: { class: NBodyDynamics, frame: MOON_CI, bodies: [...] }
       initial_state: { class: ClassicalOE, frame: MOON_OP, a: ..., e: ..., i: ... }
     gs_manager:                     # centralized OD estimator
       class: GroundStationManager
       application: { class: GroundStationManagerApp, target: sat, run_srif: true }
     DSS14:                          # a tracking station (sensor)
       class: GroundStation
       latitude_deg: 35.4
       longitude_deg: 243.1
       altitude_m: 1001
       application:
         class: GroundStationTrackingApp
         target: sat                 # resolved via GetSimulation()->GetAgent("sat")
         manager: gs_manager         # where observations are pushed
         elevation_mask_deg: 10
         use_range: true
         range_sigma_m: 10.0

Run it from Python (the pattern in every ``exN_run_*.py`` script):

.. code-block:: python

   import yaml, pylupnt as pnt
   sim = pnt.Simulation(yaml.safe_load(open("configs/ground_station_odts.yaml")))
   sim.run()
   results = sim.get_agent("gs_manager").get_application().get_results()

Authoring a new simulation
--------------------------

Path A — config-driven (preferred)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use this when the scenario fits the agent model. **Template:**
``configs/ground_station_odts.yaml`` + ``python/examples/ex7_run_odts.py``.

#. **Add the shared environment.** Write a ``world:`` block (``frame`` +
   ``force_model``, or a ``gravity:`` / ``dem:`` block for surface scenarios);
   both truth and estimator dynamics derive from it.
#. **Reuse existing classes.** List ``agents:`` with already-registered
   ``class:`` values (``Satellite``, ``GroundStation``, ``GroundStationManager``,
   ``Rover``, ``SurfaceStation``, ``IslSatellite``, ...), each with a
   ``dynamics:`` block (e.g. ``NBodyDynamics``), ``devices:``, and an
   ``application:`` block. Launch with ``pnt.Simulation(config)``.
#. **New Application** (the most common extension) — subclass ``Application``,
   implement ``MyApp(Config&)``, ``Setup()`` (resolve target agents via
   ``GetAgent``, reach the environment via ``GetAgent()->GetWorld()``, seed the
   filter; the base ``Setup`` schedules the steps), and ``Step(Real)`` (build the
   ``Measurement``, run the filter). Add ``REGISTER_FACTORY_CLASS(Application,
   MyApp)``. Reference:
   ``applications/ground_station/ground_station_manager_app.{h,cc}``.
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

The skeleton of the two most common extensions — a new ``Application`` and a new
``Measurement`` — is:

.. code-block:: cpp

   // --- A new navigation application ---------------------------------------
   class MyOdtsApp : public Application {
    public:
     MyOdtsApp(Config& cfg) : Application(cfg) {
       target_ = cfg["target"].as<std::string>();      // agent name to track
       sigma_  = cfg["range_sigma_m"].as<double>();
     }
     void Setup() override {
       Application::Setup();                            // schedules Step at frequency_
       target_agent_ = GetAgent()->GetSimulation()->GetAgent(target_);
       // seed filter_ state x_hat_ and covariance P_ here
     }
     void Step(Real t) override {
       Cart6 x_tx = target_agent_->GetStateAt(t);       // truth geometry
       MatXd H; MyMeasurement meas(BuildConfig(x_tx));
       MeasData zHR = meas.Compute(x_hat_, &H);         // z = h(x), R, and H
       filter_.Predict(t);                              // STM + process noise
       filter_.Update(zHR.value, H, zHR.covariance);    // innovation, gain, update
       Log(t);
     }
   };
   REGISTER_FACTORY_CLASS(Application, MyOdtsApp)        // YAML: class: MyOdtsApp

   // --- A new observable ---------------------------------------------------
   class MyMeasurement : public MeasurementClone<MyMeasurement> {
    public:
     MeasData Compute(const State& x, MatXd* H = nullptr) const override {
       VecXd z(1);   z(0) = /* h(x): range, Doppler, ... */;
       if (H) { H->resize(1, x.size()); *H = /* dh/dx */; }
       return MeasData{timestamp_, z, R_};              // value, covariance
     }
   };

Path B — tightly-coupled logic in one Application
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use this for tightly-coupled estimators, error-state INS, or distributed
filters. Rather than a separate ``Simulation`` subclass, the whole per-epoch body
lives in a single factory-registered ``Application`` that is still driven by the
generic engine (``pnt.Simulation(config)``). Two shapes recur:

* **Centralized coordinator.** One manager agent hosts an app that runs every
  filter itself. **Template:** ``IslOdtsCoordinatorApp`` on an ``IslOdtsManager``
  agent (``applications/lunar_sat_odts/isl_odts_coordinator_app.{h,cc}``), driven
  by ``configs/isl_odts.yaml``.
* **Fully distributed.** Each platform is its own agent hosting its own app that
  generates its measurements, runs its own filter, and exchanges state over the
  ``Simulation::Publish`` / ``Subscribe`` bus. **Template:** ``IslSatellite`` +
  ``SatelliteOdtsApp`` and ``SurfaceStationManager`` + ``GroundOdtsApp``, driven
  by ``configs/isl_odts_distributed.yaml``.

Because ``Simulation(config)`` runs ``Setup()`` in its constructor, any state a
host must inject after construction (e.g. a lander reference trajectory) should be
applied through a lazy ``Initialize`` inside the app's first ``Step``.

Reference files
---------------

``simulations/simulation.{h,cc}``, ``simulations/world.{h,cc}``, ``core/event.h``,
``core/asset_factory.h``, ``agents/agent.h``, ``applications/application.h``,
``applications/lunanet_sat_app.h``, ``devices/device.h``,
``measurements/measurement.h``, ``dynamics/dynamics.h``, ``states/state.h`` +
``states/joint_state.h``, ``configs/ground_station_odts.yaml``,
``python/examples/ex7_run_odts.py``.

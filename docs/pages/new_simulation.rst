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

One engine, one path
--------------------

**Every** scenario in LuPNT is built from a YAML config by the same generic
event-driven engine, ``lupnt::Simulation`` (``cpp/lupnt/simulations/simulation.{h,cc}``).
Agents, applications, devices, and dynamics are assembled by name through
factories; a top-level ``world:`` block defines the shared read-only environment
(epoch, integration frame, force model) that both the truth agents and the
estimators draw from; and any new mission / navigation logic goes into a
factory-registered ``Application`` that self-drives from the ``world:`` block.
There is no longer a "monolithic" alternative — nothing subclasses ``Simulation``
or hand-writes a ``Run()`` loop; the ground-station ODTS, GNSS ODTS, ISL ODTS
(centralized *and* distributed), ephemeris, surface-rover and lander scenarios
all run this one way (see ``configs/*.yaml``).

Launch a config from Python with ``pnt.Simulation(yaml_or_dict)`` (the pattern
used by every ``python/examples/exN_run_*.py`` script) or from a small C++
driver such as the tutorials under ``cpp/examples/tutorials/`` (which just build a
``Simulation`` from the YAML and call ``Run()``).

The only real design choice is **how you decompose the per-epoch logic across
apps**, and even a heavyweight Monte-Carlo engine can live inside a single app
(``LunarGnssOdtsApp`` calls ``RunLunarGnssODTSMonteCarlo`` from its ``Step``):

* **Sensor / estimator split** — separate sensor apps push observations to a
  manager app that estimates (ground-station ODTS:
  ``GroundStationTrackingApp`` → ``GroundStationManagerApp``).
* **Single coordinator app** — one manager agent hosts an app that runs the
  whole estimator itself (centralized ground ODTS: ``SurfaceStationManager`` +
  ``GroundOdtsApp``, a single EKF over the surface-station beacons).
* **Fully distributed** — every platform is its own agent running its own app,
  exchanging state over the ``Publish`` / ``Subscribe`` bus (distributed ISL:
  ``Spacecraft`` + ``SatelliteOdtsApp``, plus ``SurfaceStation`` +
  ``StationBeaconSensor`` → ``SurfaceStationManager`` + ``GroundOdtsApp``).
* **Several apps on one platform** — one physical agent hosts multiple apps in an
  ``applications:`` list, run in order and sharing the agent's truth state (lander
  descent: a guidance ``LanderGncApp`` that owns the truth trajectory + a
  navigation ``LanderNavApp`` MEKF that reads it, resolved via
  ``GetApplicationByName``).

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
       model (so truth and estimator share the same physics by default), and it
       exposes a point-mass ``Gravity()``, an optional DEM/terrain service for
       surface scenarios, and an optional ``plasma:`` block (the shared
       ionosphere/plasmasphere signal-delay environment, read by the GNSS ODTS
       app). It never propagates agents — it only *provides* the environment.
   * - **Agent**
     - ``Agent`` → ``AgentWithDynamics`` (``agents/agent.h``)
     - A platform (satellite, ground station, rover, lander, surface station).
       Owns a ``Dynamics``, ``Devices``, a ``State``, and hosts **one or more**
       ``Application`` s (an ``application:`` map, or an ``applications:`` list run
       in order). Answers ``GetStateAt(t)`` for measurement geometry.
   * - **Application**
     - ``Application`` (``applications/application.h``)
     - The mission / navigation logic run per step. Builds ``Measurement`` s and
       runs a ``Filter``. Compose several on one agent by listing them under
       ``applications:`` (e.g. a lander's guidance ``LanderGncApp`` + navigation
       ``LanderNavApp``), or via the sub-app pattern (``LunaNetSatApp`` +
       ``LunaNetSubApp``).
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
An app that reuses ``MakeDynamics()`` for *both* the truth grid and the filter has
no truth/filter model mismatch by construction; when you want a mismatch you set
the two sides separately (see `Different dynamics for truth and filter`_ below).
For surface scenarios the ``world:`` block also accepts a point-mass ``gravity:``
body and a ``dem:`` terrain block (``World::GetElevation``, ``EnuToWorld``,
``SiteCenterWorld``). The World is **read-only**: it provides the environment and
a truth facade, but each ``AgentWithDynamics`` still self-propagates via its own
``dynamics_``.

Different dynamics for truth and filter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Truth and filter dynamics are **independent knobs**, so you can (and for a
realistic OD study usually should) run the estimator against a deliberately
mismodeled truth:

* **Truth** is whatever propagates the physical agent. A self-propagating
  ``AgentWithDynamics`` (``Satellite``, ``Spacecraft``, ...) integrates its own
  ``dynamics:`` block; an app reads that truth through ``GetStateAt(t)`` /
  ``World::GetStateAt(name, t)`` (which just forwards to ``agent->GetStateAt(t)``).
* **Filter** dynamics are built *inside* the estimator ``Application`` — either by
  reusing the shared ``World::MakeDynamics()`` (the ``world: force_model:`` block)
  or from app-specific config. Whatever the filter uses also supplies its STM, so
  keep ``autodiff: true`` on that model.

*Same model (no mismatch).* Give the target agent a ``dynamics:`` block that
matches ``world: force_model:`` and have the app build its filter from
``World::MakeDynamics()`` — this is exactly what ex7 does
(``configs/ground_station_odts.yaml``).

*One force-model schema everywhere.* Every place that specifies an orbit force model —
the shared ``world: force_model:``, a physical agent's ``dynamics:``, and an estimator's
filter force model — uses the **same** ``force_model: { integrator, autodiff, bodies: [...] }``
block: a ``bodies:`` list of ``BODY: { n, m }`` (omit ``n``/``m`` or use ``{}`` for point-mass),
plus optional ``relativity`` and SRP (``mass``/``area``/``CR``). There is no separate
``moon_gravity_degree_*`` / ``include_earth`` dialect — filters and agents read the identical
schema (parsed by ``ParseForceModelSpec``).

*Different models (mismatch).* Configure the two sides separately: the **truth**
force model lives on the physical agent (or the shared ``world:`` block), while the
**filter** fidelity is an app knob. In the ISL ODTS scenario the truth is propagated
from ``world: force_model:`` — inherited by every ``Spacecraft`` — at a high gravity
resolution, and each ``SatelliteOdtsApp`` runs its EKF at a lower one
(``configs/isl_odts_distributed.yaml``):

.. code-block:: yaml

   world:
     force_model:                            # truth propagated at 16x16 lunar gravity
       bodies:
         - MOON: { n: 16, m: 16 }
         - EARTH: {}
         - SUN: {}
   # ... every Spacecraft inherits that truth force model; the onboard app filters coarser
   # using the SAME force_model schema:
   application:
     class: SatelliteOdtsApp
     force_model:                            # EKF runs at 8x8 -> intentional mismatch
       bodies:
         - MOON: { n: 8, m: 8 }
         - EARTH: {}
         - SUN: {}
       relativity: true

More generally: give the truth agent a high-fidelity ``dynamics:`` block (more
third bodies, higher gravity degree/order, SRP, drag) and point the filter at a
coarser model — a lower-fidelity ``world: force_model:`` consumed by
``World::MakeDynamics()``, or a separate app-level filter config as above. In a
custom app, simply build two ``Dynamics`` objects (one for the truth grid, one for
``filter_``) instead of sharing ``dynamics_``.

An estimator app can expose the filter model as config, too. ``GroundStationManagerApp``
(ex7) accepts an optional ``filter_dynamics:`` block (an ``NBodyDynamics``-style force
model applied to **both** its batch and its sequential SRIF/smoother filters), or
``batch_dynamics:`` / ``sequential_dynamics:`` to set each estimator independently.
Absent, both inherit the shared ``world: force_model:`` (truth == filter, the default).

Agent — a platform
~~~~~~~~~~~~~~~~~~~

``Agent`` owns ``name_``, a back-pointer to the ``Simulation``, a ``State``, a
map of ``Devices``, and one or more ``Application`` s (an ``application:`` map or an
``applications:`` list; retrieve with ``GetApplication()`` /
``GetApplications()`` / ``GetApplicationByName("...")``). Its pure-virtual
``Cart6 GetStateAt(Real t) const`` is the contract every agent must answer:
*"where am I at time t?"* — used by measurement models for light-time geometry.

Almost every physical agent derives from ``AgentWithDynamics``, which adds a
``Ptr<Dynamics>``, an ``AttitudeDynamics``, and a ``Cart6 state_``. Its
``Step(t)`` calls ``Propagate(t)`` (integrating the dynamics in place) then
``Log(t)``. Concrete agents (all registered with the factory except ``Lander``):
``Satellite``, ``Spacecraft`` (orbit + clock 8-state truth), ``GroundStation``,
``Rover``, ``SurfaceStation``, ``Lander``. Group generators:
``Constellation`` / ``GnssConstellation`` (from the ``constellations:`` block) and
``LunarNavConstellation`` — a single ``agents:`` entry that expands into **N child
``Spacecraft``** so a scenario need not declare one agent per navigation satellite,
from either an explicit ``satellites:`` list or a symmetric ``walker:`` spec
(``n_planes``, ``sats_per_plane``, frozen ELFO ``a/e/i/omega``); the surface-rover
example (ex10) sources its relay truth from it via ``GetSatelliteStateAt(j, t)``.

Application — the mission logic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Application`` holds a back-pointer to its owning ``Agent``. ``Setup()``
schedules a periodic ``Step(t)`` (at ``APPLICATION`` priority) and seeds the
filter state/covariance; the pure-virtual ``Step(Real t)`` runs one
predict/update cycle. An agent can host **several** applications (an
``applications:`` list, run in insertion order and resolved by
``GetApplicationByName``); you can also compose logic *within* one app via the
**sub-app pattern**: ``LunaNetSatApp`` holds a ``std::vector<Ptr<LunaNetSubApp>>``
and fans its ``Step`` out to each sub-app (e.g. ``IslOdtsApp`` and
``EphemerisGenApp``).

Concrete apps (all factory-registered): the **sensor / estimator split** for
ground-station OD — ``GroundStationTrackingApp`` on each ``GroundStation``
(a sensor: visibility-gated range / range-rate, pushed to the manager) feeding a
``GroundStationManagerApp`` on a ``GroundStationManager`` agent (the centralized
batch + SRIF/smoother estimator, with configurable filter dynamics); ``LunarGnssOdtsApp``
(GNSS ODTS); the distributed ``SatelliteOdtsApp`` (per-satellite onboard filter) with
``StationBeaconSensor`` → ``GroundOdtsApp`` (the centralized station-only ground filter);
``EphemerisApp`` / ``LunaNetSatApp`` (sub-app host) + ``IslOdtsApp`` / ``EphemerisGenApp``
(sub-apps); ``SurfaceStationApp``; the surface-rover INS ``SurfaceRoverNavApp`` (which
sources its relay truth from a ``LunarNavConstellation`` agent); and the lander, whose
descent is split across **two apps on one ``Lander`` agent** — a guidance ``LanderGncApp``
that owns the truth trajectory and writes the lander state, and a navigation
``LanderNavApp`` (error-state INS MEKF) that reads that truth to synthesize its
IMU / altimeter / crater / LunaNet measurements each ``Step``.

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

**Template:** ``configs/ground_station_odts.yaml`` +
``python/examples/ex7_run_odts.py``.

#. **Add the shared environment.** Write a ``world:`` block (``frame`` +
   ``force_model``, or a ``gravity:`` / ``dem:`` block for surface scenarios);
   both truth and estimator dynamics derive from it.
#. **Reuse existing classes.** List ``agents:`` with already-registered
   ``class:`` values (``Satellite``, ``GroundStation``, ``GroundStationManager``,
   ``Rover``, ``SurfaceStation``, ``Spacecraft``, ...), each with a
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

Working in Python
-----------------

The whole simulation lifecycle is driven from Python — you **compose a scenario and
run it**, reusing the registered C++ agent/app/measurement/dynamics classes. A scenario
is just a ``dict`` (or a loaded YAML), so it can be built, mutated, and swept
programmatically:

.. code-block:: python

   import copy, yaml, pylupnt as pnt

   cfg = yaml.safe_load(open("configs/ground_station_odts.yaml"))
   cfg["agents"]["gs_manager"]["application"]["filter"]["process_accel_sigma_mps2"] = 1e-8

   sim = pnt.Simulation(cfg)                 # instantiates World + agents + apps from `class:` keys
   sim.run()                                 # runs the event loop
   app = sim.get_agent("gs_manager").get_application()
   results = app.get_results()               # each app binds its own result accessors

Monte-Carlo trials run at the simulation level via ``pnt.run_monte_carlo`` (parallel across
processes; each trial deep-copies the config with a distinct top-level ``seed``), or — for a
scenario with an expensive shared precompute — via the engine's in-config
``monte_carlo_runs`` (one precompute, ``seed = base_seed + i`` per run). See
``python/examples/ex6_gnss_odts.ipynb`` §7d for both.

**Authoring a new Application in Python.** You can subclass ``pnt.Application`` in pure Python
and have the C++ simulation drive it — no C++ build required. Register the class with
``pnt.register_application(name, cls)`` and reference it by ``class: name`` in a config; the
engine constructs it as ``cls(config_dict)`` (the ``application:`` block, as a dict) and calls
its ``step(t)`` (and ``setup()`` / ``log(t)``) each epoch. The app reads truth through its host
agent, builds a filter dynamics model, and stores results on ``self`` (``sim.get_agent(...)
.get_application()`` hands the same Python object back):

.. code-block:: python

   import numpy as np, yaml, pylupnt as pnt

   class MyOdtsApp(pnt.Application):
       def __init__(self, config):
           pnt.Application.__init__(self)
           self.target = config["target"]
           self.set_frequency(config.get("frequency", 1.0 / 60))
           self.est = []
       def setup(self):
           pnt.Application.setup(self)          # base schedules step() at get_frequency()
           self.dyn = pnt.NBodyDynamics()       # a filter force model, built in Python
           for b in (pnt.Body.Moon(8, 8), pnt.Body.Earth(), pnt.Body.Sun()):
               self.dyn.add_body(b)
           self.dyn.set_frame(pnt.Frame.MOON_CI); self.dyn.set_autodiff(True)
       def step(self, t):
           world = self.get_agent().get_world()
           r_obs = self.get_agent().get_state_at(t)             # this agent's truth [r; v]
           r_tgt = world.get_state_at(self.target, t)           # any agent's truth by name
           xf, F = self.dyn.propagate_stm(self.x, t0, t)        # predict with STM (numpy in/out)
           # ... EKF update from a measurement computed in numpy; store on self ...
           self.est.append(self.x.copy())

   pnt.register_application("MyOdtsApp", MyOdtsApp)
   cfg = yaml.safe_load(open("configs/sat_bearing_odts.yaml"))
   cfg["agents"]["observer"]["application"] = {"class": "MyOdtsApp", "target": "target"}
   sim = pnt.Simulation(cfg); sim.run()
   est = sim.get_agent("observer").get_application().est     # results read straight off self

A complete, converging example — an angles-only OD EKF authored this way — is
``python/examples/py_authored_app_demo.py``.

**Authoring a new Agent in Python.** A physical platform (whose truth trajectory you define in
Python — an analytic ephemeris, a scripted path) works the same way: subclass ``pnt.Agent``,
implement ``get_state_at(self, t)`` returning a numpy ``[r; v]`` 6-vector, and register it with
``pnt.register_agent(name, cls)``:

.. code-block:: python

   class PyEphemeris(pnt.Agent):
       def __init__(self, config):
           pnt.Agent.__init__(self)
           self.a = config["radius_m"]; self.w = (pnt.GM_MOON / self.a**3) ** 0.5
       def get_state_at(self, t):
           c, s = np.cos(self.w * t), np.sin(self.w * t)
           return np.array([self.a*c, self.a*s, 0, -self.a*self.w*s, self.a*self.w*c, 0])

   pnt.register_agent("PyEphemeris", PyEphemeris)   # config: {class: PyEphemeris, radius_m: 2.0e6}

**Authoring a new Measurement in Python.** An observable model is a ``pnt.Measurement`` subclass
implementing ``compute(self, x) -> (z, H, R)`` (the predicted observable, its Jacobian
``H = dh/dx``, and the noise ``R``, all numpy). ``measurement.evaluate(x)`` runs it through the
C++ base — the same path a filter uses — so the *same* Python class both **generates** an
observation (apply it to a truth state from the target agent/device, add noise) and **predicts**
one (apply it to the filter state):

.. code-block:: python

   class BearingMeasurement(pnt.Measurement):
       def __init__(self, sigma_rad):
           pnt.Measurement.__init__(self); self.sigma = sigma_rad; self.target = np.zeros(3)
       def compute(self, x):
           d = self.target - x[:3]; rn = np.linalg.norm(d); u = d / rn
           H = np.zeros((3, 6)); H[:, :3] = -(np.eye(3) - np.outer(u, u)) / rn
           return u, H, self.sigma**2 * np.eye(3)

   m = BearingMeasurement(np.radians(2 / 3600)); m.target = tgt_truth[:3]
   z = m.evaluate(obs_truth)[0] + rng.normal(0, m.sigma, 3)   # generate from the target's truth
   u, Hx, R = m.evaluate(self.x)                              # predict at the filter state

The angles-only demo (``python/examples/py_authored_app_demo.py``) uses exactly this. So new
**Applications, Agents, and Measurements** are all authorable in pure Python; new ``Dynamics``
classes remain C++ (reuse ``pnt.NBodyDynamics`` from Python).

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

Decomposing tightly-coupled logic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a tightly-coupled estimator, error-state INS, or distributed filter, the same
config-driven path still applies — the whole per-epoch body simply lives in one or
more factory-registered ``Application`` s driven by ``pnt.Simulation(config)``
(there is no ``Simulation`` subclass). How you spread the work across apps is the
only choice, and three shapes recur:

* **Centralized coordinator.** One manager agent hosts an app that runs the whole
  estimator itself. **Template:** ``SurfaceStationManager`` + ``GroundOdtsApp``
  (``applications/lunar_sat_odts/ground_odts_app.{h,cc}``) — a single EKF over the
  surface-station beacons — or ``GroundStationManager`` + ``GroundStationManagerApp``
  (ex7), driven by ``configs/isl_odts_distributed.yaml`` /
  ``configs/ground_station_odts.yaml``.
* **Fully distributed.** Each platform is its own agent hosting its own app that
  generates its measurements, runs its own filter, and exchanges state over the
  ``Simulation::Publish`` / ``Subscribe`` bus. **Template:** ``Spacecraft`` +
  ``SatelliteOdtsApp`` with ``SurfaceStation`` + ``StationBeaconSensor``, driven
  by ``configs/isl_odts_distributed.yaml``.
* **Several apps on one platform.** One physical agent hosts multiple
  ``applications:``, run in order and sharing the agent's truth state. **Template:**
  a ``Lander`` hosting a guidance ``LanderGncApp`` (owns the descent truth) + a
  navigation ``LanderNavApp`` MEKF that resolves the guidance app via
  ``GetApplicationByName``, driven by ``configs/lander_nav.yaml``.
* **Heavyweight engine inside one app.** An existing batch/Monte-Carlo engine can
  be called straight from a single app's ``Step`` — e.g. ``LunarGnssOdtsApp``
  invokes ``RunLunarGnssODTSMonteCarlo`` — so no rewrite of the numerics is needed.

Because ``Simulation(config)`` runs ``Setup()`` in its constructor, any state a
host must inject after construction (e.g. the lander guidance app's reference
trajectory) should be applied through a lazy ``Initialize`` inside the app's first
``Step``.

Reference files
---------------

``simulations/simulation.{h,cc}``, ``simulations/world.{h,cc}``, ``core/event.h``,
``core/asset_factory.h``, ``agents/agent.h``, ``applications/application.h``,
``applications/lunanet_sat_app.h``, ``devices/device.h``,
``measurements/measurement.h``, ``dynamics/dynamics.h``, ``states/state.h`` +
``states/joint_state.h``, ``configs/ground_station_odts.yaml``,
``python/examples/ex7_run_odts.py``.

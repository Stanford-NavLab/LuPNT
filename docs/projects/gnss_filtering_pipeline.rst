GNSS Filtering Pipeline
===================================================================

The ``projects/GNSS_Filtering`` scenario configures the core
``LunarGnssODTSSimulation`` implemented under
``cpp/lupnt/simulations/LunarGnssODTS``. The project folder contains scenario
assets only: YAML configs, Python delay precompute, plotting scripts, and the
notebook. The workflow is organized as a staged pipeline so that expensive
environment simulations can be cached and inspected separately from the Monte
Carlo filter.

The default case is a 3 minute ELFO receiver simulation around the Moon using
GPS SP3/ANTEX data. The receiver application is scheduled in the receiver's
local clock at ``receiver_app.rate_hz`` and each app tick is mapped to the TDB
coordinate timeline used for propagation and measurement generation.

Pipeline Stages
-------------------------------------------------------------------

1. ``LunarGnssODTSSimulation::Precompute`` in C++
   Builds the GNSS constellations from SP3/ANTEX data, propagates the nominal
   lunar receiver trajectory, solves receiver-dependent light time, applies
   transmitter PCO corrections through the constellation setup, and writes one
   row per epoch, transmitter, and frequency. Ionosphere/plasma columns are set
   to zero in this stage.

2. ``precompute_delays.py`` in Python
   Reads the link CSV and runs GCPM/IRI ray tracing for each saved link. This
   step is intentionally Python multiprocessing based because the raytrace
   stack is process-parallel, not OpenMP/thread-parallel.

3. ``LunarGnssODTSSimulation::Run`` in C++
   Loads the precomputed links and, when enabled, merges the delay CSV into the
   truth channels. It then generates noisy pseudorange, Doppler, and optional
   TDCP measurements and runs the UDU EKF. When TDCP is enabled, the C++ stage
   uses UDU stochastic cloning with state ``[x_k, x_{k-1}]`` so the
   time-differenced carrier range can depend on both receiver epochs. Truth and
   filter propagation use ``NBodyDynamics``. The truth clock is propagated with
   ``JointOrbitClockDynamics`` and OCXO process noise; relativistic clock drift
   is evaluated from the receiver orbit. SRP can be modeled as a fixed force or
   estimated as an additional filter state.

4. Post-processing
   The notebook and plotting script visualize precomputed plasma delays,
   tracked satellite counts, RTN position/velocity errors, clock errors, SRP
   coefficient error when enabled, and covariance-derived 3-sigma bounds.

Run Commands
-------------------------------------------------------------------

The full pipeline can be launched with:

.. code-block:: bash

   pixi run run-gnss-pipeline

The GCPM delay batch can be slow. Skip it only when the configured delay table
already exists, or when the config has ``plasma.simulate_truth: false``:

.. code-block:: bash

   pixi run run-gnss-pipeline --skip-delays

The wrapper accepts:

.. code-block:: bash

   --config PATH
   --workers N
   --serial-delays
   --overwrite-delays
   --skip-delays
   --no-plot

The stages can also be run one at a time:

.. code-block:: bash

   pixi run precompute-gnss-links
   pixi run precompute-gnss-delays
   pixi run run-gnss-monte-carlo
   pixi run plot-gnss-filtering

Inputs
-------------------------------------------------------------------

The main configuration is:

.. code-block:: text

   projects/GNSS_Filtering/gnss_filtering_config.yaml

Receiver hardware, tracking-loop, link-budget, measurement-noise, and
filter-mismatch defaults are selected from:

.. code-block:: text

   projects/GNSS_Filtering/gnss_designs.yaml

Measurement selection is controlled in the ``measurements`` block. The default
uses pseudorange and Doppler. Set ``measurements.use_tdcp: true`` to add TDCP
rows after the first receiver app tick; ``carrier_phase_sigma_m`` and
``tdcp_sigma_m`` set the carrier/TDCP noise levels in meters.

The constellation configuration can auto-select SP3 files from
``constellation.sp3_directory`` for the configured start epoch and simulation
duration. ANTEX PCO data are read from ``constellation.antex_file``. The
default scenario uses every GPS PRN available in the selected SP3 file; set
``constellation.use_all_gps: false`` and list ``constellation.gps_prns`` for a
small debugging subset. Galileo can be enabled with
``constellation.include_galileo: true`` when the selected SP3/ANTEX products
contain Galileo satellites.

Intermediate Files
-------------------------------------------------------------------

By default, staged products are written under ``output/gnss_filtering``:

.. list-table::
   :header-rows: 1

   * - File
     - Producer
     - Purpose
   * - ``precomputed_links.csv``
     - ``LunarGnssODTSSimulation::Precompute``
     - Receiver-specific light-time links with zero plasma delay columns.
   * - ``precomputed_delays.csv``
     - ``precompute_delays.py``
     - GCPM/IRI delay values keyed by epoch, constellation, PRN, and frequency.
   * - ``trajectory_mc<N>.csv``
     - ``LunarGnssODTSSimulation::Run``
     - Per-run truth, estimate, RTN errors, covariance bounds, clock errors,
       optional SRP error, and tracked satellite count.
   * - ``summary.csv``
     - ``LunarGnssODTSSimulation::Run``
     - Final and RMS Monte Carlo statistics.

Truth and Estimated Measurement Split
-------------------------------------------------------------------

The truth measurement channels are built from interpolated SP3/ANTEX
transmitter data, receiver-dependent light-time iteration, Sun-only Shapiro
delay, and optional precomputed GCPM/IRI ionosphere/plasma delay. The filter
measurement uses the ephemeris metadata carried in each ``GnssChannel`` and
re-solves light time from the current estimated receiver state. This keeps the
truth data product and the estimated measurement function close to the way a
receiver consumes a navigation message.

The filter can either ignore plasma delay and inflate measurement noise with
``plasma.filter_pseudorange_noise_inflation_m`` /
``plasma.filter_doppler_noise_inflation_hz``, or model the same delay carried
by the channel when ``plasma.model_in_filter`` is enabled.

Time and Unit Conventions
-------------------------------------------------------------------

The receiver trajectory and app schedule are propagated on a TDB coordinate
timeline. GNSS constellation ephemeris tables are stored and interpolated in
TAI seconds. The measurement API converts receiver receive epochs from
``options.receive_time_scale`` to ``options.ephemeris_time_scale`` before
interpolating transmitter states.

Receiver states use SI position and velocity units. Clock bias and drift can
be represented in seconds, meters, or kilometers through ``ClockBiasUnit``; the
filter example uses seconds in its configuration, while the library supports
range-like clock units for better numerical scaling in other filters.

Outputs and Visualization
-------------------------------------------------------------------

Open the notebook:

.. code-block:: text

   projects/GNSS_Filtering/plot_gnss_filtering_results.ipynb

or regenerate PNGs from the command line:

.. code-block:: bash

   pixi run plot-gnss-filtering

The command-line script writes:

.. code-block:: text

   output/gnss_filtering/gnss_filtering_rtn_tracking.png
   output/gnss_filtering/gnss_filtering_plasma_delays.png

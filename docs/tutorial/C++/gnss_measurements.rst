GNSS Measurement Simulation
===================================================================

This page outlines the receiver-side GNSS simulation flow in C++.
`GnssConstellation` stores GNSS satellite ephemerides and transmitter
metadata. `GNSSMeasurements` owns the receiver-dependent work: light-time
iteration, visibility, relativistic range corrections, plasma delays, channel
selection, and observable generation.

Time and State Conventions
-------------------------------------------------------------------

Epochs passed to ``GNSSMeasurements::Compute`` and
``GNSSMeasurements::Precompute`` are receiver signal-reception epochs. The
scale is explicit. A cislunar simulation will often use TDB for the receiver
coordinate timeline while keeping GNSS ephemeris tables in TAI:

.. code-block:: cpp

   GnssMeasurementOptions options;
   options.receive_time_scale = Time::TDB;     // receive epochs you pass in
   options.ephemeris_time_scale = Time::TAI;   // GnssConstellation tables
   options.frame = Frame::MOON_CI;             // common inertial state frame

The receiver state must contain position, velocity, clock bias, and clock
drift at the configured indices:

.. code-block:: cpp

   options.indices.position = 0;
   options.indices.velocity = 3;
   options.indices.clock_bias = 6;
   options.indices.clock_drift = 7;
   options.clock_bias_unit = ClockBiasUnit::SECONDS;

``ClockBiasUnit::METERS`` and ``ClockBiasUnit::KILOMETERS`` are also
supported when a filter state should use range-like clock units for numerical
scaling.

Constellation Setup
-------------------------------------------------------------------

The constellation is a satellite catalog. It does not decide which satellites
are visible to a receiver.

.. code-block:: cpp

   auto gps = MakePtr<GnssConstellation>(GnssConst::GPS);
   gps->SetSatelliteStates(prns, t_tai, rv_eci_by_prn);
   gps->SetupTransmitters();

``SetSatelliteStates`` stores the tabulated ephemerides in TAI seconds and
fits a piecewise-Chebyshev state model. Runtime satellite interpolation uses
the Chebyshev model.

Online Flow
-------------------------------------------------------------------

Online computation is for filters and closed-loop simulations where the
receiver state is known one epoch at a time.

.. code-block:: cpp

   GNSSMeasurements gnss(gps);
   gnss.SetFrequency(GnssFreq::L1);
   gnss.SetOptions(options);

   GnssOccludingBody earth;
   earth.radius_m = R_EARTH;
   earth.position_m = Vec3::Zero();
   gnss.SetOccludingBodies({earth});

   gnss.SetSunPositionProvider([](Real t_rx) {
     return Vec3(0.0, AU, 0.0);
   });
   gnss.SetBoresightTargetProvider([](Real t_rx) {
     return Vec3::Zero();
   });

   MatXd H;
   GNSSMeasurementsEpoch epoch = gnss.Compute(t_rx, receiver_state, &H);

For each PRN, ``Compute``:

1. Converts the receiver receive epoch to the ephemeris scale.
2. Iterates transmit time with geometric light time.
3. Builds a ``GnssChannel`` containing transmit state, receive/transmit epochs,
   ephemeris, C/N0, and correction terms.
4. Applies receiver-dependent visibility and C/N0 thresholding.
5. Computes pseudorange, Doppler, carrier phase, covariance, and optional
   Jacobian.

Relativistic and Plasma Corrections
-------------------------------------------------------------------

The measurement model includes transmitter clock relativity, Shapiro delay,
and optional ionosphere/plasma delays.

.. code-block:: cpp

   options.apply_transmitter_relativity = true;

   options.apply_shapiro_delay = true;
   options.shapiro_mu = GM_EARTH;
   options.shapiro_body_position = Vec3::Zero();

   options.apply_ionosphere_plasma_delay = true;
   options.default_ionosphere_plasma_delay_m = 0.0;

   GnssIonospherePlasmaRayTraceOptions plasma;
   plasma.config.kp = kp_index;
   plasma.config.step_size = step_size_km;
   plasma.config.cutoff_r = cutoff_radius_km;
   plasma.config.gradn_dx = gradient_step_km;
   plasma.config.integ_method = integration_method;
   plasma.config.correction_method = correction_method;
   plasma.config.correction = apply_ray_correction;
   plasma.config.fine_correction = apply_fine_correction;
   plasma.config.corr_tol = correction_tolerance_m;
   plasma.config.use_fortran_gcpm = use_fortran_gcpm;
   plasma.config.compute_higher_order = compute_higher_order;
   plasma.config.use_adaptive_step = use_adaptive_step;
   plasma.config.straight_ray = use_straight_ray;
   plasma.raytrace_frame = Frame::ECEF;
   plasma.raytrace_epoch_scale = Time::UTC;
   plasma.iri_model = pecsim::IRIModel::IRI_2020;
   gnss.SetIonospherePlasmaRayTraceOptions(plasma);

Shapiro delay is non-dispersive and is added to both pseudorange and carrier
phase range. Ionosphere/plasma delay is represented as a positive code delay
in meters: pseudorange uses ``+delay`` and carrier phase uses ``-delay``.
The GNSS raytrace hook calls ``pecsim::trace_ray`` and uses the returned TEC
delay. The caller owns the plasma simulation settings; the measurement class
does not pick Kp, step size, cutoff radius, integrator, or correction policy.

For cislunar GNSS filtering runs, the recommended high-throughput workflow is
to precompute all receiver-specific links first, then run the expensive
GCPM/IRI delay calculation as a separate Python multiprocessing step (the
Precompute Flow below).

Precompute Flow
-------------------------------------------------------------------

Precompute is for generating a measurement sequence before filtering. It uses
the same channel-building logic as online computation, but it can call a batch
ionosphere/plasma model after all epochs and links are known.

.. code-block:: cpp

   std::vector<Real> receive_times = {t0, t1, t2};
   std::vector<State> receiver_states = {x0, x1, x2};

   gnss.SetBatchCustomIonospherePlasmaDelayModel(
       [](const std::vector<Real>& times,
          const std::vector<State>& states,
          const std::vector<std::vector<GnssChannel>>& channels) {
         std::vector<std::vector<Real>> delays(times.size());
         for (size_t i = 0; i < times.size(); i++) {
           delays[i].resize(channels[i].size());
           for (size_t j = 0; j < channels[i].size(); j++) {
             delays[i][j] = 0.0;  // batch model result in meters
           }
         }
         return delays;
       });

   std::vector<GNSSMeasurementsEpoch> epochs =
       gnss.Precompute(receive_times, receiver_states, true);

The batch provider receives:

* receiver receive epochs in ``options.receive_time_scale``
* receiver states in ``options.frame``
* receiver-specific channels with light-time transmit epochs already solved

This is the right hook for ionosphere/plasma models that need the full
simulation interval, endpoint geometry, cached environmental grids, or shared
line-of-sight integrations. When the batch provider is installed,
``Precompute`` skips the online plasma hook during channel construction and
fills each channel delay from the batch result.

Using Precomputed Data in a Filter
-------------------------------------------------------------------

Each ``GNSSMeasurementsEpoch`` contains the channels and the stacked
observables:

.. code-block:: cpp

   for (const GNSSMeasurementsEpoch& epoch : epochs) {
     const VecXd& y = epoch.values;
     const MatXd& R = epoch.covariance;
     const MatXd& H = epoch.jacobian;

     for (const GnssChannel& channel : epoch.channels) {
       int prn = channel.prn;
       Real t_tx = channel.transmit_time;
       Real shapiro_m = channel.shapiro_delay_m;
       Real plasma_m = channel.ionosphere_plasma_delay_m;
     }

Ionosphere, Plasmasphere, and Ray Tracing
===================================================================

Purpose
-------------------------------------------------------------------

This specification defines the mathematical contract for the ionospheric and
plasmaspheric group-delay model applied to terrestrial GNSS signals received
by lunar (cislunar) users.  Unlike near-Earth links, lunar GNSS paths traverse
long, near-horizontal trajectories through the ionosphere and plasmasphere;
the delay is therefore computed by numerically ray-tracing the bent signal
path and integrating the electron density along it.

The primary modeling reference is Keidai Iiyama and Grace Gao, *"Ionospheric
and Plasmaspheric Delay Characterization for Lunar Terrestrial GNSS Receivers
with Global Core Plasma Model,"* submitted to NAVIGATION: Journal of the
Institute of Navigation, arXiv:2510.10059 (2025), and the author's PhD thesis
Chapter 5.  Equation numbers ``(5.x)`` below refer to that chapter.

The ray tracer lives in the ``pecsim`` namespace under
``cpp/lupnt/environment/plasma/tec/raytrace.{h,cc}`` (``pecsim::trace_ray``,
``RayTraceConfig``, ``PathProfile``).  Electron density comes from GCPM v2.4
(``cpp/lupnt/environment/plasma/gcpm/gcpm_interface.{h,cc}``) wrapping the IRI
ionosphere (``gcpm/iri_interface.{h,cc}``); the magnetic field comes from
IGRF-14 (``cpp/lupnt/environment/plasma/igrf/igrf_interface.{h,cc}``).  The
delay is consumed by the GNSS measurement model through
``GNSSMeasurements::ComputeIonospherePlasmaRayTraceDelay``
(``cpp/lupnt/measurements/gnss_measurement.cc``), specified as
:math:`\Delta\rho_\mathrm{plasma}` in :doc:`gnss_measurements`.

.. note::

   The entire plasma module is **km-based**.  Positions ``x1``, ``x2`` passed
   to ``trace_ray`` are in km, the step size and cutoff radius are in km, the
   Earth radius is ``RE = 6371.0`` km, and the speed of light is
   ``C = 299792.458`` km/s (``cpp/lupnt/environment/plasma/core/constants.h``).
   Electron density is evaluated internally in cm\ :sup:`-3` (GCPM) and
   converted to m\ :sup:`-3` before entering the delay integrals; one TEC unit
   is ``TECU = 1e16`` electrons/m\ :sup:`2`.  Final delays are returned in
   meters.  Callers such as ``GNSSMeasurements`` divide their SI (meter)
   positions by 1000 before calling ``trace_ray``.

Coordinate, Unit, and Time Contract
-------------------------------------------------------------------

The ray tracer works in a geocentric Earth-fixed Cartesian frame.  The GNSS
coupling converts receiver/transmitter states from ``options.frame`` into
``raytrace_frame`` (default ``Frame::ECEF``) and passes an epoch in
``raytrace_epoch_scale`` (default ``Time::UTC``), expressed as seconds from
J2000:

.. math::

   x_1 = r_S^\mathrm{ecef}\,[\mathrm{km}], \qquad
   x_2 = r_R^\mathrm{ecef}\,[\mathrm{km}], \qquad
   t = t_\mathrm{rx}^\mathrm{UTC}\ (\text{s since J2000}).

Electron-density evaluation converts the geocentric position into **solar
magnetic (SM)** coordinates required by GCPM
(``cpp/lupnt/environment/plasma/tec/raytrace.cc :: compute_ne``):

.. code-block:: cpp

   std::array<double,3> pos_gei_arr = {pos_geo[0]/RE, pos_geo[1]/RE, pos_geo[2]/RE};
   geo_to_sm(itime, pos_gei_arr, pos_sm);       // geographic -> solar magnetic
   cart_to_pol(pos_sm, alatr, along, r);        // r in RE, alatr in rad
   double amlt = 12.0 + along / AMLTRAD;         // magnetic local time [hours]
   if (amlt > 24.0) amlt -= 24.0;
   out = gcpm_v24_fortran(datetime, r, amlt, alatr, kp);   // out[0] = ne [cm^-3]

The SM frame (thesis, following the plasma-physics convention) has
:math:`X_\mathrm{SM}` in the Sun-Earth plane, :math:`Z_\mathrm{SM}` aligned
with the geomagnetic pole, and :math:`Y_\mathrm{SM}` completing the right-hand
system.

Electron Density Model
-------------------------------------------------------------------

The electron density is supplied by the **Global Core Plasma Model (GCPM
v2.4)**, an empirical model of the inner magnetosphere that stitches together
region-specific submodels for the ionosphere (based on the International
Reference Ionosphere, IRI), the plasmasphere, the plasmapause, the
magnetospheric trough, and the polar cap, with hyperbolic-tangent transitions.
As a mathematical contract, the density along the path is

.. math::

   n_e = n_e\big(r,\ \phi_\mathrm{SM},\ \mathrm{MLT};\ K_p,\ R_{12},\ t\big),

where :math:`r` is geocentric radius, :math:`\phi_\mathrm{SM}` is SM latitude,
and :math:`\mathrm{MLT}` is magnetic local time.  The inputs are:

* :math:`K_p` -- planetary geomagnetic index (``config.kp``).  It sets the
  plasmapause radius: higher :math:`K_p` erodes the plasmasphere, contracting
  the plasmapause from ~4-6 :math:`R_E` (quiet) to ~2-3 :math:`R_E` (active).
* :math:`R_{12}` / :math:`F_{10.7}` -- solar activity (``config.rz12``).  The
  IRI submodel derives :math:`F_{10.7}` and the ionospheric index
  :math:`IG_{12}` from :math:`R_{12}` via the PHaRLAP relations (thesis
  Eqs. 5.19-5.20):

  .. math::

     F_{10.7} = 63.75 + 0.728\,R_{12} + 0.00089\,R_{12}^2,

  .. math::

     IG_{12} = -12.349154 + 1.4683266\,R_{12} - 0.00267690893\,R_{12}^2.

* :math:`t` -- epoch (date/time), which sets IRI solar-activity parameters and
  the SM frame orientation.

GCPM returns total electron density in cm\ :sup:`-3`;
``compute_ne`` converts to m\ :sup:`-3`:

.. math::

   n_e\,[\mathrm{m}^{-3}] = 10^{6}\, n_e\,[\mathrm{cm}^{-3}].

.. note::

   ``config.kp`` must be set to a non-negative value; ``compute_ne`` throws if
   ``kp < 0``, and an IRI model must be selected (``set_iri_model``) or it
   throws.  ``config.rz12`` follows IRI's convention: ``> 0`` uses the value
   directly (:math:`0 < R_{12} \le 200`), ``-1`` uses IRI's
   historical/projected :math:`R_{12}` with the storm model, and ``-2`` the
   same without the storm model.  ``trace_ray`` applies ``config.rz12`` to the
   active IRI option (``IRI2007Option``/``IRI2020Option``) before tracing.
   GCPM/IRI model only the inner-magnetospheric cold-plasma environment;
   contributions outside the cutoff radius (default :math:`4\,R_E`) are
   neglected, where densities are much lower.

.. note::

   GCPM cannot be evaluated for arbitrary future epochs (it requires
   historical solar/geomagnetic parameters).  The thesis baseline uses
   2025-01-01 12:00 UTC as a proxy epoch for a 2027 scenario; the
   characterization spans :math:`K_p \in \{1,3,5,6.667,9\}` and
   :math:`R_{12} \in \{10,50,100,150,200\}`.

Magnetic Field Model
-------------------------------------------------------------------

The magnetic field, needed only for the second- and third-order terms, is
evaluated with IGRF-14 in geographic spherical coordinates
(``raytrace.cc :: compute_B``):

.. code-block:: cpp

   cart_to_pol(pos_geo_arr, lat_rad, lon_rad, r);
   std::vector<double> B = igrf14(lat_rad, lon_rad, r, decimal_year);   // nT

The projection of the field on the ray direction enters the delay through
:math:`\cos\theta`, the angle between the ray path and the field line:

.. math::

   B = \lVert \mathbf{B} \rVert, \qquad
   \cos\theta = \frac{\mathbf{B}\cdot\hat{s}}{B}.

.. note::

   The delay coefficients in the code assume :math:`B` in **nanotesla**.  This
   is the source of the numeric difference from the thesis symbolic
   coefficients, which are written for :math:`B` in tesla (see the note under
   *Refractive Index and Group Delay*).

Refractive Index and Group Delay
-------------------------------------------------------------------

For a right-hand circularly polarized GNSS signal the phase refractive index
:math:`n` and group refractive index :math:`n_\mathrm{gr}` are (thesis
Eqs. 5.1-5.2)

.. math::

   n = 1
     - \frac{f_p^2}{2f^2}
     - \frac{f_p^2 f_g \cos\theta}{2f^3}
     - \frac{f_p^2}{4f^4}\!\left[\frac{f_p^2}{2} + f_g^2(1+\cos^2\theta)\right],

.. math::

   n_\mathrm{gr} = 1
     + \frac{f_p^2}{2f^2}
     + \frac{f_p^2 f_g \cos\theta}{f^3}
     + \frac{3f_p^2}{4f^4}\!\left[\frac{f_p^2}{2} + f_g^2(1+\cos^2\theta)\right],

with plasma frequency :math:`f_p = n_e e^2/(4\pi\epsilon_0 m)` and gyro
frequency :math:`f_g = eB/(2\pi m)`.  The delay is dispersive (:math:`\propto
1/f^2` at leading order).  The ionospheric/plasmaspheric group delay along the
path (thesis Eqs. 5.4-5.7) is

.. math::

   d_{I_\mathrm{gr}} = \int_\Gamma (n_\mathrm{gr}-1)\,ds
       = \frac{p}{f^2} + \frac{q}{f^3} + \frac{u}{f^4},

.. math::

   p = 40.3 \int_\Gamma n_e\, ds = 40.3\,\mathrm{TEC}
     = 40.3\,(\mathrm{TEC}_\mathrm{LOS} + \Delta\mathrm{TEC}_\mathrm{bend}),

.. math::

   q = -2.2566\times10^{12} \int_\Gamma n_e\, B \cos\theta\, ds,

.. math::

   u = 2437 \int_\Gamma n_e^2\, ds
     + 4.74\times10^{22} \int_\Gamma n_e\, B^2 (1+\cos^2\theta)\, ds .

The integral of :math:`n_e` along the path is the **TEC**.  Because the path
:math:`\Gamma` is curved, the TEC splits into a line-of-sight part
:math:`\mathrm{TEC}_\mathrm{LOS}` and a bending increment
:math:`\Delta\mathrm{TEC}_\mathrm{bend}` (the difference between the bent-path
and straight-line TEC).  The first-order term dominates; the second and third
orders are below ~1% of the first for terrestrial and LEO users.

The phase index used to bend the ray is
``raytrace.cc :: refractive_index_neB``, and the accumulators for the three
group-delay orders are set in ``raytrace.cc :: ray_derivative``:

.. code-block:: cpp

   // refractive_index_neB (phase index that bends the ray)
   double p = p_coeff * ne_m3;                          // 40.3 * ne
   double q = q_coeff * ne_m3 * B * cos_theta;          // 2.2566e3 * ne*B*cos (B in nT)
   double u = u1_coeff*ne_m3*ne_m3
            + u2_coeff*ne_m3*B*B*(1 + cos_theta*cos_theta);
   double n = 1 - (p/f2 + q/(2*f3) + u/(3*f4));

   // ray_derivative: per-arc group-delay accumulators (s in km -> *1000)
   dz(7) = ne_m3 * 1000;                                // d(TEC)/ds
   dz(8) = q_coeff * ne_m3 * B * cos_theta * 1000;      // 2nd-order integrand
   dz(9) = 1000*(u1_coeff*ne_m3*ne_m3
              + u2_coeff*ne_m3*B*B*(1 + cos_theta*cos_theta));   // 3rd-order integrand

where the coefficients are ``p_coeff = 40.3``, ``q_coeff = 2.2566e3``,
``u1_coeff = 2437``, ``u2_coeff = 4.74e4``.

.. note::

   The code coefficients are scaled for :math:`B` in **nanotesla**, so
   ``q_coeff = 2.2566e3`` and ``u2_coeff = 4.74e4`` correspond to the thesis
   SI-tesla values :math:`2.2566\times10^{12}` and :math:`4.74\times10^{22}`
   (factors :math:`10^{-9}` and :math:`10^{-18}`).  The code stores
   ``q`` with a positive coefficient and applies the sign of the delay through
   the phase index (:math:`-q/(2f^3)` in ``refractive_index_neB``) and through
   the sign of :math:`\cos\theta`; the thesis absorbs the sign into the
   symbolic definition of :math:`q` in Eq. (5.6).  Both the phase index
   (:math:`q/(2f^3)`, :math:`u/(3f^4)`) and the group-delay accumulators
   (:math:`q/f^3`, :math:`u/f^4`) match Eqs. (5.1)-(5.4) -- the group form omits
   the :math:`\tfrac12` and :math:`\tfrac13` phase factors.

Ray-Tracing Algorithm
-------------------------------------------------------------------

Signal propagation is governed by the vector eikonal equation (thesis
Eqs. 5.11-5.13):

.. math::

   \frac{dt}{ds} = \frac{n}{c}, \qquad
   \frac{d\mathbf{r}}{ds} = \hat{s}, \qquad
   \frac{d\hat{s}}{ds} = \frac{\nabla n - \hat{s}(\hat{s}\cdot\nabla n)}{n},

subject to the transmitter/receiver boundary conditions
:math:`\mathbf{r}(t_\mathrm{tx})=\mathbf{r}_\mathrm{tx}`,
:math:`\mathbf{r}(t_\mathrm{rx})=\mathbf{r}_\mathrm{rx}`.  Given the reception
epoch/position, the transmit epoch is solved by the standard light-time
iteration (thesis Eq. 5.16), consistent with :doc:`gnss_measurements`.

The state vector integrated along arc length :math:`s` is
:math:`z = [\,\mathbf{r},\ \hat{s},\ t,\ \mathrm{TEC},\ I_2,\ I_3\,]` (10
components), and the refractive-index gradient is a central finite difference
with step ``config.gradn_dx`` (default 1 km).  The eikonal right-hand side is
``raytrace.cc :: ray_derivative``:

.. code-block:: cpp

   dz.segment<3>(0) = s_hat;                                   // dr/ds = s_hat
   dz.segment<3>(3) = (grad_n - s_hat*grad_n.dot(s_hat)) / n;  // ds_hat/ds  (Eq. 5.13)
   dz(6) = 1.0 / C;                                            // dt/ds ~ 1/c

.. note::

   The geometric propagation time is accumulated as ``dz(6) = 1/C``
   (vacuum) rather than ``n/c``; the ionospheric delay is then added
   separately as ``(tec_delay + second_delay + third_delay)/c``.  This is a
   deliberate machine-precision choice (the code comment notes that the
   per-step TEC time increment is :math:`\sim n_e\cdot 10^{-23}`\ s, too small
   to co-accumulate with the epoch), and is numerically equivalent to
   integrating :math:`dt/ds = n/c`.

**Shooting method with cutoff split.**  The path is divided into a
plasmasphere section and a vacuum section at the cutoff radius
``config.cutoff_r`` (default :math:`4\,R_E \approx 25484` km): the ray is
assumed straight once it exits the cutoff.  Integration marches from the
transmitter with the fourth-order Runge-Kutta or Euler integrator
(``config.integ_method``) until :math:`r \ge` cutoff, then shoots straight to
the receiver (thesis Eqs. 5.17-5.18):

.. math::

   s_f = s_\mathrm{exit} + \Delta s, \qquad
   \Delta s = (\mathbf{x}_\mathrm{rx} - \mathbf{x}_\mathrm{exit})\cdot
              \hat{s}_\mathrm{exit},

.. code-block:: cpp

   // propagate_ray: exit at cutoff, project remaining distance onto s_hat_exit
   Vec3d delta_x2 = x2 - pos_now;
   double delta_s = delta_x2.dot(dir_now);
   sf = s + delta_s;
   Vec3d final_pos = pos_now + dir_now * (sf - s);

An **adaptive step size** is used inside the plasmasphere
(``raytrace.cc :: adjust_stepsize``), keyed on altitude
:math:`\mathrm{alt}=r-R_E`:

.. code-block:: cpp

   if (alt < 1000)      step_size = 10.0;   // km
   else if (alt < 4000) step_size = 20.0;
   else                 step_size = 100.0;

**Initial-direction correction (the iterative bent-path solve).**  The unknown
initial ray direction :math:`\hat{s}_0` is adjusted so the shot lands on the
receiver.  ``trace_ray`` runs up to 20 *coarse* iterations that (a) refresh the
propagation time from the accumulated delay
(``correct_proptime``) and (b) re-optimize :math:`\hat{s}_0` over the
2-D azimuth/elevation of the initial direction with a Nelder-Mead simplex
minimizing the terminal miss :math:`\lVert\bar{\mathbf{x}}_\mathrm{rx} -
\mathbf{x}_\mathrm{rx}\rVert` (``correct_dir_neldermead``).  To keep cost low,
the coarse pass reuses the stored per-step refractive index and its gradient
(``store_mat``) instead of re-querying GCPM at every point.  An optional finer
pass (``config.fine_correction``) refines with Nelder-Mead or a
finite-difference Newton solve (``correct_dir_newton``) to the tolerance
``config.corr_tol`` (default 1 m).  Without correction the terminal error is
10-100 km due to bending; with correction it is driven to sub-meter-100 m.

.. code-block:: cpp

   // trace_ray: coarse loop -- proptime update then direction re-solve
   correct_proptime(prop_time, dir, sf, x1, x2, epoch_utc_rx, config, zf);
   correct_dir_neldermead(prop_time, dir, sf, x1, x2, epoch_utc_rx, config, zf,
                          store_mat, /*coarse_mode=*/true, debug_corr);

.. note::

   ``config.correction_method`` must be ``"neldermead"`` or ``"newton"``;
   ``trace_ray`` throws otherwise.  The struct default is ``"grid"`` (which
   would throw), so a caller enabling correction must set this field.  The
   thrown message text (``"Supported methods are 'naive' and 'newton'"``) is
   stale relative to the actual check.

Delay Decomposition (``PathProfile``)
-------------------------------------------------------------------

After the corrected path is integrated, ``compute_path_profile`` assembles the
delay terms in meters from the accumulators (``z(7..9)``) and the geometry:

.. code-block:: cpp

   double tec_delay_m    = p_coeff * z(7) / freq2;          // 40.3 * TEC_bent / f^2
   double second_delay_m = z(8) / freq3;
   double third_delay_m  = z(9) / freq4;
   pp.dist_bend_m        = (sf - dist) * 1000;              // excess path length
   pp.tec_delay_bend_m   = tec_delay_m - tec_delay_straight_m;

These map onto the thesis decomposition (Eqs. 5.8-5.10).  The excess path
length from bending is

.. math::

   d_\mathrm{len} = \int_\Gamma ds - \rho, \qquad
   \rho = \lVert \mathbf{r}_\mathrm{tx} - \mathbf{r}_\mathrm{rx} \rVert,

and the total group delay relative to the LOS distance is

.. math::

   d_\mathrm{total}
     = d_{I_\mathrm{gr}} + d_\mathrm{len}
     = \underbrace{40.3\,\mathrm{TEC}_\mathrm{LOS}}_{d_{I1,\mathrm{LOS}}}
     + \underbrace{q/f^3}_{d_{I2}}
     + \underbrace{u/f^4}_{d_{I3}}
     + \underbrace{40.3\,\Delta\mathrm{TEC}_\mathrm{bend}}_{d_{I1,B}}
     + d_\mathrm{len}.

The ``PathProfile`` fields carry each term: ``tec_delay_m``
(:math:`40.3\,\mathrm{TEC}` along the **bent** path, i.e.
:math:`d_{I1,\mathrm{LOS}}+d_{I1,B}`), ``tec_delay_bend_m`` (:math:`d_{I1,B}`),
``second_delay_m`` (:math:`d_{I2}`), ``third_delay_m`` (:math:`d_{I3}`),
``dist_bend_m`` (:math:`d_\mathrm{len}`), and ``total_delay_m``
(:math:`d_\mathrm{total}`).

Configuration Knobs (``RayTraceConfig``)
-------------------------------------------------------------------

The complete set of knobs (``raytrace.h``), with defaults:

* ``freq_Hz`` -- carrier frequency [Hz]; L1/L2/L5 constants are provided
  (``freq_L1 = 1575.42e6`` etc.).  All delays scale as :math:`1/f^2` to leading
  order, so L5 (1176.45 MHz) sees ~1.8x the L1 delay.
* ``step_size = 10.0`` -- nominal integrator step [km].
* ``use_adaptive_step = true`` -- enable the 10/20/100 km altitude schedule.
* ``integ_method = "Euler"`` -- ``"Euler"`` or ``"RK4"``.
* ``correction = true`` -- run the bent-path direction correction.
* ``fine_correction = false`` -- run the additional fine pass.
* ``correction_method = "grid"`` -- must be set to ``"neldermead"`` or
  ``"newton"`` when correction is on (see note above).
* ``corr_tol = 1.0`` -- terminal miss tolerance [m].
* ``cutoff_r = 4*RE`` -- plasmasphere/vacuum boundary [km].
* ``gradn_dx = 1.0`` -- central-difference step for :math:`\nabla n` [km].
* ``kp = -1`` -- geomagnetic :math:`K_p` (must be set :math:`\ge 0`).
* ``rz12 = -1.0`` -- solar :math:`R_{12}` (see IRI convention above).
* ``compute_higher_order = true`` -- include second/third orders (needs IGRF).
* ``straight_ray = false`` -- assume a straight path (no bending; skips coarse
  correction, zeros :math:`\nabla n`), for the straight-line TEC baseline.
* ``use_fortran_gcpm = true`` -- use the Fortran GCPM v2.4 backend.

Coupling into the GNSS Measurement Model
-------------------------------------------------------------------

When ``options.apply_ionosphere_plasma_delay`` is true and ray-trace options
have been registered via ``SetIonospherePlasmaRayTraceOptions``, the GNSS
model calls the tracer per link
(``gnss_measurement.cc :: ComputeIonospherePlasmaRayTraceDelay``):

.. code-block:: cpp

   pecsim::RayTraceConfig config = ray_opts.config;
   if (ray_opts.use_channel_frequency) config.freq_Hz = GnssFrequencyHz(freq);
   if (ray_opts.iri_model != pecsim::IRIModel::NONE) set_iri_model(ray_opts.iri_model);
   Real t_raytrace = ConvertGnssEpoch(receive_time, options_.receive_time_scale,
                                      ray_opts.raytrace_epoch_scale);
   // r_rx, r_tx converted options_.frame -> ray_opts.raytrace_frame, then to km
   pecsim::PathProfile profile = pecsim::trace_ray(
       t_raytrace.val(), ToVec3d(r_tx_ray)/1000.0, ToVec3d(r_rx_ray)/1000.0, config);
   Real delay_m = profile.tec_delay_m;
   if (ray_opts.delay_mode == GnssIonospherePlasmaRayTraceDelayMode::TEC_PLUS_HIGHER_ORDER)
     delay_m += profile.second_delay_m + profile.third_delay_m;

The returned :math:`\Delta\rho_\mathrm{plasma}` is the positive code delay of
:doc:`gnss_measurements`, entering the pseudorange as
:math:`+\Delta\rho_\mathrm{plasma}` and the carrier phase as
:math:`-\Delta\rho_\mathrm{plasma}`.  The ``GnssIonospherePlasmaRayTraceOptions``
fields are: ``config`` (the tracer knobs; ``GNSSMeasurements`` never fills
:math:`K_p`, integrator, or step sizes), ``raytrace_frame`` (default
``ECEF``), ``raytrace_epoch_scale`` (default ``UTC``), ``use_channel_frequency``
(default true -- overwrite ``config.freq_Hz`` with the channel frequency),
``iri_model``, and ``delay_mode``.

.. note::

   ``delay_mode`` defaults to ``TEC_ONLY``, so the returned delay is
   ``profile.tec_delay_m`` (the first-order bent-path TEC delay,
   :math:`d_{I1,\mathrm{LOS}}+d_{I1,B}`) alone; ``TEC_PLUS_HIGHER_ORDER`` adds
   :math:`d_{I2}+d_{I3}`.  The geometric excess-path term
   :math:`d_\mathrm{len}` (``profile.dist_bend_m``) is computed but **not**
   added to the pseudorange delay; its contribution to timing is already the
   :math:`\Delta\mathrm{TEC}_\mathrm{bend}` part of ``tec_delay_m``, and the
   sub-meter geometric term is left out of the range equation.

Batch vs Online Precompute
-------------------------------------------------------------------

Ray tracing is expensive (a GCPM call per step per iteration), so two paths
exist.  The **online** path calls ``trace_ray`` link-by-link inside
``Compute``.  The **batch/precompute** path decouples the tracer from filtering:
in the staged Lunar GNSS ODTS scenario, ``Precompute`` writes all
light-time-corrected links with zero plasma delay, a Python driver
(``precompute_delays.py``) fills each epoch/channel delay with GCPM/IRI ray
tracing, and ``Run`` merges the delay file back before forming observables (see
the *Online and Precompute Equivalence* section of :doc:`gnss_measurements`).
The batch tracer is exposed through
``cpp/lupnt/environment/plasma/tec/raytrace_wrapper.{h,cc}``
(``pecsim::run_raytrace_batch``, ``load_raytrace_data`` /
``RaytraceData``), which stores per-ray ``total_delays``, ``tec_delays``,
``second_delays``, ``third_delays``, ``dist_bend_m``, ``tec_delay_bend_m``,
``tecus``, and ``min_alts``.  The :math:`K_p` time series consumed by these
runs is fetched by ``python/pylupnt/plasma/kp_loader.py`` from GFZ Potsdam and
cached under the plasma data directory.

Model Boundaries and Notes
-------------------------------------------------------------------

* **Km-unit convention.**  All tracer inputs/outputs are km / km-s based; SI
  callers scale by 1000 at the boundary.
* **Magnetic-field units.**  Delay coefficients assume :math:`B` in nT; the
  thesis symbolic coefficients are the SI-tesla equivalents.
* **Second/third order are small.**  For terrestrial and LEO users they are
  below ~1% of the first-order delay; the thesis lunar characterization finds
  them several orders of magnitude below :math:`d_{I1}` (Table 5.2), and the
  GNSS default (``TEC_ONLY``) omits them.
* **Cutoff at** :math:`4\,R_E`.  Plasma outside the inner magnetosphere is
  neglected and the path is straight beyond the cutoff.
* **Time integration.**  ``dt/ds`` uses the vacuum :math:`1/c` with the delay
  added separately, equivalent to :math:`n/c` but numerically robust.
* **Excess path length** :math:`d_\mathrm{len}` is characterized but not added
  to the measurement delay.
* **Validity ranges.**  GCPM/IRI require historical solar/geomagnetic inputs;
  the characterization covers :math:`K_p \in [1,9]` and
  :math:`R_{12}\in[10,200]`, and delays are strongly altitude-dependent
  (mean first-order delay ~50 m for tangential altitudes below 500 km, falling
  below ~0.2 m above 10,000 km).

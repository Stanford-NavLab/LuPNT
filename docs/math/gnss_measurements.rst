GNSS Measurement Model
===================================================================

Purpose
-------------------------------------------------------------------

This specification defines the GNSS receiver measurement model implemented by
``GnssMeasurement`` and ``GNSSMeasurements``.  The model generates
pseudorange, Doppler, and carrier phase observations from receiver state,
GNSS transmitter ephemerides, receiver/transmitter clock corrections, and
optional propagation delays.

State, Time, and Frame Contract
-------------------------------------------------------------------

All vectors in this specification are expressed in ``options.frame``.  The
current C++ implementation requires receiver and transmitter states to share a
defined common inertial frame; the lunar filtering project uses
``Frame::MOON_CI`` for receiver states and GNSS channel states after
constellation setup.

The receiver state vector is interpreted through
``GnssMeasurementOptions::indices``:

.. math::

   x =
   \begin{bmatrix}
     \cdots &
     r_R^\mathsf{T} &
     \cdots &
     v_R^\mathsf{T} &
     \cdots &
     b_R &
     d_R &
     \cdots
   \end{bmatrix}^\mathsf{T},

where

.. math::

   r_R \in \mathbb{R}^3 \quad [\mathrm{m}], \qquad
   v_R \in \mathbb{R}^3 \quad [\mathrm{m/s}].

The receiver clock bias ``b_R`` and drift ``d_R`` are represented in
``options.clock_bias_unit``.  They are converted to SI seconds before entering
the range equations:

.. math::

   \Delta t_R = \operatorname{seconds}(b_R), \qquad
   \dot{\Delta t}_R = \operatorname{seconds}(d_R).

Epochs passed to ``Compute`` and ``Precompute`` are receiver signal-reception
epochs:

.. math::

   t_R \equiv t_\mathrm{receive}.

The scale of ``t_R`` is ``options.receive_time_scale``.  GNSS constellation
ephemerides are tabulated in TAI seconds and are evaluated in
``options.ephemeris_time_scale``; the current implementation enforces TAI for
the ephemeris scale. Cislunar simulations commonly set the receiver receive
epoch scale to TDB and let the measurement model convert receive epochs to TAI
for GNSS ephemeris interpolation.

GNSS Transmitter State Interpolation
-------------------------------------------------------------------

For each PRN, ``GnssConstellation::SetSatelliteStates`` receives tabulated
states

.. math::

   \{(t_k,\; r_S(t_k),\; v_S(t_k))\}_{k=0}^{N-1}.

The constellation fits a piecewise Chebyshev model for the six-dimensional
state:

.. math::

   y_S(t) =
   \begin{bmatrix}
     r_S(t) \\
     v_S(t)
   \end{bmatrix}
   \approx
   \sum_{n=0}^{p} c_n T_n(\tau),

where

.. math::

   \tau = \frac{t - t_\mathrm{mid}}{t_\mathrm{radius}},
   \qquad
   -1 \le \tau \le 1.

Runtime satellite state evaluation uses this Chebyshev model.  The original
state table may still be carried in each ``GnssChannel`` as ephemeris-message
metadata.

Light-Time Iteration
-------------------------------------------------------------------

For a receiver position ``r_R(t_R)`` and transmitter PRN ``s``, the transmit
epoch is solved iteratively in the ephemeris time scale:

.. math::

   t_S^{(0)} = t_R,

.. math::

   t_S^{(i+1)}
   =
   t_R
   -
   \frac{
     \left\lVert r_R(t_R) - r_S(t_S^{(i)}) \right\rVert
   }{c}.

The iteration stops when

.. math::

   \left| t_S^{(i+1)} - t_S^{(i)} \right|
   <
   \epsilon_\mathrm{lt},

or when ``options.light_time_max_iterations`` is reached.  The convergence
tolerance is ``options.light_time_tolerance_s``.

Clock terms and group delays are not used to solve the physical light-time
epoch; they enter the measurement equations below.

Geometry
-------------------------------------------------------------------

After light-time convergence, define

.. math::

   r_S = r_S(t_S), \qquad
   v_S = v_S(t_S),

.. math::

   \Delta r = r_R - r_S, \qquad
   \Delta v = v_R - v_S,

.. math::

   \rho = \lVert \Delta r \rVert,
   \qquad
   \hat{u} = \frac{\Delta r}{\rho},

.. math::

   \dot{\rho}
   =
   \frac{\Delta r^\mathsf{T}\Delta v}{\rho}.

Transmitter Clock Corrections
-------------------------------------------------------------------

The transmitter effective clock bias is

.. math::

   \Delta t_S^\mathrm{eff}
   =
   \Delta t_S
   +
   \Delta t_S^\mathrm{rel}
   -
   \Delta t_S^\mathrm{grp},

where ``tx_clock_bias_s`` is :math:`\Delta t_S`,
``relativistic_correction_s`` is :math:`\Delta t_S^\mathrm{rel}`, and
``group_delay_s`` is :math:`\Delta t_S^\mathrm{grp}`.

The current transmitter relativistic correction is the state-based broadcast
GNSS form

.. math::

   \Delta t_S^\mathrm{rel}
   =
   -\frac{2 r_S^\mathsf{T} v_S}{c^2}.

The transmitter effective clock drift is

.. math::

   \dot{\Delta t}_S^\mathrm{eff}
   =
   \dot{\Delta t}_S.

Propagation Delays
-------------------------------------------------------------------

Shapiro Delay
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``options.apply_shapiro_delay`` is true, the one-way gravitational delay
from the Sun is represented as a range correction:

.. math::

   \Delta \rho_\mathrm{Shapiro}
   =
   \frac{2\mu_\odot}{c^2}
   \ln
   \left(
     \frac{
       r_R^\odot + r_S^\odot + \rho
     }{
       r_R^\odot + r_S^\odot - \rho
     }
   \right),

where

.. math::

   r_R^\odot = \lVert r_R - r_\odot(t_R) \rVert, \qquad
   r_S^\odot = \lVert r_S - r_\odot(t_R) \rVert.

The Sun position :math:`r_\odot(t_R)` is supplied by
``SetSunPositionProvider`` or by the default Sun provider.

Ionosphere / Plasma Delay
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ionosphere/plasma correction is represented as a positive code delay in
meters:

.. math::

   \Delta \rho_\mathrm{plasma} \ge 0.

It is disabled by default.  When
``options.apply_ionosphere_plasma_delay`` is true, the delay is supplied by one
of ``SetCustomIonospherePlasmaDelayModel`` (custom online link-by-link),
``SetIonospherePlasmaRayTraceOptions`` (built-in plasma ray tracer in
``environment/plasma/tec/raytrace``), ``SetBatchCustomIonospherePlasmaDelayModel``
(precompute over all epochs/channels), or the static
``options.default_ionosphere_plasma_delay_m``.  The ray-traced TEC delay uses

.. math::

   \Delta \rho_\mathrm{plasma}
      = \frac{40.3}{f^2}\int_\Gamma n_e \, ds ,

and the full electron-density / ray-bending model is specified in
:doc:`ionosphere_plasmasphere`.  The sign convention is dispersive:

.. math::

   \text{code range uses } +\Delta \rho_\mathrm{plasma},
   \qquad
   \text{carrier phase range uses } -\Delta \rho_\mathrm{plasma}.

Observable Equations
-------------------------------------------------------------------

Let

.. math::

   \lambda = \frac{c}{f}.

The receiver clock range and range-rate terms are

.. math::

   \Delta \rho_R = c \Delta t_R, \qquad
   \Delta \dot{\rho}_R = c \dot{\Delta t}_R.

The transmitter clock range and range-rate terms are

.. math::

   \Delta \rho_S = c \Delta t_S^\mathrm{eff}, \qquad
   \Delta \dot{\rho}_S = c \dot{\Delta t}_S^\mathrm{eff}.

Pseudorange
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The pseudorange observable is

.. math::

   P
   =
   \rho
   +
   \Delta \rho_\mathrm{Shapiro}
   +
   \Delta \rho_\mathrm{plasma}
   +
   \Delta \rho_R
   -
   \Delta \rho_S.

Doppler
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Doppler observable in Hz is

.. math::

   D
   =
   -
   \frac{
     \dot{\rho}
     +
     \Delta \dot{\rho}_R
     -
     \Delta \dot{\rho}_S
   }{\lambda}.

The current implementation does not include time derivatives of Shapiro or
plasma delay in Doppler.

Carrier Phase
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The non-integer carrier phase observable in cycles is

.. math::

   \Phi
   =
   \frac{
     \rho
     +
     \Delta \rho_\mathrm{Shapiro}
     -
     \Delta \rho_\mathrm{plasma}
     +
     \Delta \rho_R
     -
     \Delta \rho_S
   }{\lambda}
   +
   \phi_0,

where :math:`\phi_0` is ``phase_bias_cycles``.

If the measurement state includes an integer ambiguity :math:`N`, the emitted
carrier observable is

.. math::

   \Phi_\mathrm{meas} = \Phi + N.

Visibility and Channel Selection
-------------------------------------------------------------------

Visibility is receiver-dependent and is evaluated in ``GNSSMeasurements``.
For each occluding body, the line of sight is retained only if the spherical
body does not block the receiver-transmitter segment.

For near-surface points, visibility is evaluated using an elevation threshold

.. math::

   e > -10^\circ.

For space-to-space links, visibility is evaluated using a tangent-line
occlusion test.  C/N0 thresholding is then optionally applied:

.. math::

   C/N_0 > (C/N_0)_\mathrm{min}.

The geometric visibility test and its two regimes are specified in
:doc:`link_budget` (``ComputeVisibility``).

Link Budget and C/N0
-------------------------------------------------------------------

For each surviving channel, ``GNSSMeasurements::ComputeCN0``
(``cpp/lupnt/measurements/gnss_measurement.cc :: ComputeCN0``) evaluates the
received carrier-to-noise density in dB-Hz.  The link-budget equation and the
free-space path-loss, EIRP, gain-pattern, and tracking-loop-noise details are
specified in full in :doc:`link_budget`; this section covers only
how the *variables that feed the budget* are computed inside the GNSS model.

.. code-block:: cpp

   GnssAttitude::Compute(r_tx_gcrf, v_tx_gcrf, r_sun_gcrf, ex, ey, ez);
   Real phi_tx   = safe_acos(u_tx2rx_gcrf.dot(ez));
   Real theta_tx = atan2(u_tx2rx_gcrf.dot(ey), u_tx2rx_gcrf.dot(ex));
   Real G_tx = constellation->GetTransmitterAntenna(prn, freq)
                   .ComputeGain(theta_tx, phi_tx);
   Real G_rx = rx_antenna_.ComputeGain(0.0, phi_rx);
   Real P_tx = constellation->GetTransmitPowerDbw(prn, freq);
   return LinkBudget(P_tx, G_tx, G_rx, range, freq, rx_params_);

In closed form (thesis Eq. 6.59),

.. math::

   \left(C/N_0\right)_\mathrm{dB\text{-}Hz}
   =
   P_\mathrm{tx}
   + G_\mathrm{tx}(\theta_\mathrm{tx}, \varphi_\mathrm{tx})
   + G_\mathrm{rx}(\theta_\mathrm{rx})
   - L_\mathrm{fs}
   - L_\mathrm{atm}
   - L_\mathrm{ad}
   - L_\mathrm{pol}
   - 10\log_{10}(k_B)
   - 10\log_{10}(T_\mathrm{eff}).

Transmitter Attitude and Yaw Steering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The transmit gain :math:`G_\mathrm{tx}` is looked up in the GNSS satellite's
yaw-steering body frame, computed by ``GnssAttitude::Compute``
(``cpp/lupnt/attitude/gnss_attitude.cc``).  ``ComputeCN0`` first transforms the
transmitter state, receiver position, and Sun position from ``options.frame``
into ``Frame::GCRF`` before building the attitude, so the gain pattern is
evaluated in the same Earth-centered inertial frame the ANTEX data assume.
The orthonormal triad is

.. math::

   e_z = \widehat{-r_S}\ (\text{nadir/boresight}),
   \quad
   e_y = \widehat{e_z \times \widehat{(r_\odot - r_S)}}\ (\text{Sun side}),
   \quad
   e_x = \widehat{e_y \times e_z},

and the boresight angles of the transmitter-to-receiver line of sight
:math:`\hat{u}` are

.. math::

   \varphi_\mathrm{tx} = \cos^{-1}\!\big(\hat{u}^\mathsf{T} e_z\big),
   \qquad
   \theta_\mathrm{tx}
   = \operatorname{atan2}\!\big(\hat{u}^\mathsf{T} e_y,\ \hat{u}^\mathsf{T} e_x\big).

The velocity-aware overload derives this same frame from the documented
**nominal yaw-steering law** (``cpp/lupnt/attitude/gnss_yaw_steering.cc ::
NominalYawAngle``, Cheng et al. 2025, Eq. 1):

.. code-block:: cpp

   Real GnssYawSteering::NominalYawAngle(Real beta, Real mu) {
     return atan2(-tan(beta), sin(mu));            // Eq. (1)
   }

with :math:`\beta` the Sun elevation above the orbital plane
(``BetaAngle``) and :math:`\mu` the orbit angle from the midnight point
(``OrbitAngle``).  ``GnssYawSteering`` additionally provides the block-specific
maneuver laws (GPS IIF/IIR/III, Galileo IOV/FOC, BDS-3 CAST/SECM; Eqs. 3-16)
as stateless building blocks, but ``ComputeCN0`` uses only the nominal frame,
matching the thesis assumption (Chapter 6.4.3) of nominal yaw steering with
boresight along the spacecraft-Earth direction and eclipse steering left as
future work.

Antenna Gain Patterns (ANTEX / ACE)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The gains are table lookups into measured patterns,
``Antenna::ComputeGain`` (``cpp/lupnt/measurements/antenna.cc``), over the
off-boresight angle :math:`\varphi` and azimuth :math:`\theta`:

.. code-block:: cpp

   Real gain = LinearInterp2d(phi_, theta_, gain_, phi, theta);   // 2D pattern
   if (phi > phi_max_ || phi < -phi_max_) return NAN;             // sidelobe cutoff

Patterns are read by ``Antenna::LoadAntennaPattern`` from named data files and
normalized to :math:`\varphi\in[-180,180]^\circ`,
:math:`\theta\in[0,360]^\circ` by ``FormatAntennaPattern``.  Per the thesis
(Chapter 6.4.3), the transmit patterns are populated from the NASA antenna
characterization experiment (ACE) study for GPS Block II-F, Lockheed-Martin
data for Block IIR/IIR-M, the U.S. Coast Guard Navigation Center for
Block III, the EU Publications Office for Galileo, and the Cabinet Office of
Japan for QZSS; the lunar receive antenna has a peak gain of 14 dBi with a
:math:`12.2^\circ` half-power beamwidth.  The transmit power values
:math:`P_\mathrm{tx}` (thesis Table 6.3) are supplied separately by
``GnssConstellation::GetTransmitPowerDbw``.  An empty antenna name gives an
omni pattern (:math:`G\equiv 0`), and off-pattern directions return ``NaN``,
which acts as an implicit sidelobe cutoff.

IGS ANTEX (``.atx``) files are parsed by
``cpp/lupnt/interfaces/antex_loader.{h,cc}`` (``AntexLoader``, exposed to Python
as ``pnt.AntexLoader``), which exposes the antenna phase-center offsets (see
below) rather than the gain pattern.

Precise (SP3) and Broadcast (BRDC) Ephemeris
-------------------------------------------------------------------

``GnssConstellation`` consumes *precomputed* transmitter ephemerides (ECI
position/velocity history per PRN) fit with a piecewise Chebyshev model.
Those ephemerides are produced by one of two loaders, which differ in how the
transmitter position, velocity, and clock are obtained.

Precise SP3 ephemeris
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``cpp/lupnt/interfaces/sp3_loader.{h,cc}`` (``Sp3Loader``, exposed to Python as
``pnt.Sp3Loader``) parses IGS SP3 precise-ephemeris files.  SP3 tabulates satellite **center-of-mass** ECEF
positions and clock bias at a fixed cadence; ``Sp3Loader::GetPosVelClock``
returns an interpolated position/velocity/clock from a per-satellite
piecewise Chebyshev fit, with velocity taken as the analytic derivative of the
fitted position polynomial:

.. code-block:: cpp

   void Sp3Loader::GetPosVelClock(const std::string& sat_id, Real t_tai,
                                  Vec6& rv_ecef, Real& clock_bias_s) const;

Because SP3 positions are referred to the center of mass, the antenna
**phase-center offset** (PCO) must be added to place the transmit point at the
antenna.  ``cpp/lupnt/interfaces/antex_loader.{h,cc}`` (``AntexLoader``)
supplies the PCO in the satellite North/East/Up frame and applies it via the
satellite "IJK" rotation:

.. math::

   r_\mathrm{ant}^\mathrm{ecef}
   =
   r_\mathrm{CoM}^\mathrm{ecef}
   +
   C_{ijk}(t, r_\mathrm{CoM}^\mathrm{ecef})\,\Delta_\mathrm{PCO}^\mathrm{NEU},

.. code-block:: cpp

   static Vec3d AntexLoader::ApplyPcoCorrectionEcef(
       Real t_tai, const Vec3d& pos_sp3_ecef, const Vec3d& pco_neu_m);

.. note::

   The PCO "IJK" triad ``AntexLoader::ComputeIjkToEcefRotation``
   (``jvec = normalize(r_sun - r_sat)``, ``kvec = -normalize(r_sat)``,
   ``ivec = jvec x kvec``) is intentionally the *raw, un-orthogonalized*
   Sun-pointing frame used to generate the precomputed PCO-corrected
   ephemerides.  It is deliberately distinct from the orthonormal
   ``GnssAttitude`` body frame used for gain-pattern lookups; the two must not
   be conflated.

Broadcast BRDC ephemeris
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``cpp/lupnt/interfaces/rinex_nav_loader.{h,cc}`` (``RinexNavLoader``, exposed to
Python as ``pnt.RinexNavLoader``) parses RINEX V3 navigation ("BRDC") files and
evaluates the transmitter state by **Keplerian propagation of the broadcast
navigation message** nearest the requested epoch (GPS / Galileo / BeiDou /
QZSS; GLONASS, which uses tabulated state vectors, is intentionally
unsupported).  Broadcast ephemerides carry the real-time signal-in-space error
(position and clock) relative to the precise product, and are used to inject
realistic ephemeris/clock modeling errors into the measurement truth.  The
LunaNet-style broadcast-message generation on the transmit side lives in
``cpp/lupnt/applications/ephemeris/lunanet_ephemeris.*`` and
``ephemeris_gen_app.*`` (specified in :doc:`ephemeris_almanac`).

.. note::

   The thesis (Chapter 6.4.2 / 6.5.2) treats the antenna-phase-center-corrected
   IGS precise (SP3) product as truth and derives per-satellite ephemeris
   errors by differencing the broadcast (BRDC) position/clock against it; a
   per-constellation median clock-bias offset (thesis Eq. 6.58) is removed to
   absorb the broadcast/precise time-reference difference.  In LuPNT, the
   ``RinexNavLoader`` omits the Galileo-specific GST/GPST ("GAGP")
   system-time-correction term; this affects only the satellite *clock*
   (sub-100 ns, i.e. sub-30 m range-equivalent) and has no effect on the
   broadcast position or velocity.

Almanac constellation source (future-epoch scenarios)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The precise/broadcast paths above both require IGS products for the run epoch,
which do not exist for a *future* mission epoch.  For those scenarios the Lunar
GNSS ODTS driver (``constellation.source == "almanac"``,
``cpp/lupnt/simulations/lunar_gnss_odts/lunar_gnss_odts_simulation.cc``)
synthesizes each PRN's transmitter ephemeris instead of loading it:

#. **Seed elements.** Each PRN's coarse Keplerian set is read from a YUMA GPS
   almanac (``RinexNavLoader::LoadYumaFile``,
   ``cpp/lupnt/interfaces/rinex_nav_loader.{h,cc}``) or, absent a YUMA file, from
   the broadcast (BRDC) navigation message.  The seed Cartesian state is
   evaluated at the almanac reference epoch (:math:`t_k=0`) in ECEF and rotated
   to ECI.

#. **GPS week-number rollover.** A YUMA almanac carries a mod-1024 week number,
   so its raw reference epoch can land ~20 years from the intended date.  The
   seed epoch is snapped to the 1024-week era nearest the run epoch,

   .. math::

      t_\mathrm{seed}
      \leftarrow
      t_\mathrm{seed}
      +
      \Delta_\mathrm{1024}\,
      \operatorname{round}\!\left(
        \frac{t_\mathrm{run} - t_\mathrm{seed}}{\Delta_\mathrm{1024}}
      \right),
      \qquad
      \Delta_\mathrm{1024} = 1024 \times 7 \times 86400\ \text{s},

   so the propagation span is the intended months, not decades.  A full-week
   BRDC seed makes this a no-op.

#. **Numerical propagation.** From the seed state each PRN is propagated over the
   ephemeris grid under Earth gravity (a :math:`4\times4` field) plus Sun and
   Moon third-body point masses (``CreateGnssEarthDynamics``), and the resulting
   ECI states are handed to ``GnssConstellation::SetSatelliteStates`` exactly like
   a precise/broadcast fit.

Because no precise reference exists to difference against, the realistic
broadcast **signal-in-space error** (SISE) is *modeled* rather than measured: the
``SyntheticSISE`` transmitter-error model draws a per-PRN orbit error in the
local orbital (radial / along-track / cross-track) frame and a clock error,

.. math::

   \delta r_\mathrm{ECI}
   = \delta_R\,\hat r + \delta_A\,\hat t + \delta_C\,\hat c,
   \qquad
   \delta_{R,A,C} \sim \mathcal{N}\!\big(0,\ \sigma_{R,A,C}^2\big),
   \qquad
   \delta t \sim \mathcal{N}\!\big(0,\ \sigma_\mathrm{clk}^2\big),

once per PRN from a seeded RNG, held fixed over the (short) arc — a first-order
stand-in for the near-constant systematic error a real receiver sees.  The
per-PRN magnitudes :math:`\sigma_R,\sigma_A,\sigma_C,\sigma_\mathrm{clk}` come
from ``constellation.synthetic_sise_{radial,along,cross,clock}_m``, and the delta
is injected through the same filter-transmitter code path
(``BroadcastEphemerisError``) that the measured broadcast error uses.

Jacobian Convention
-------------------------------------------------------------------

The analytic Jacobian currently includes first-order derivatives of geometric
range, range rate, receiver clock bias, receiver clock drift, and optional
carrier integer ambiguity with respect to the receiver state.

For pseudorange:

.. math::

   \frac{\partial P}{\partial r_R}
   =
   \hat{u}^\mathsf{T},
   \qquad
   \frac{\partial P}{\partial b_R}
   =
   c \frac{\partial \Delta t_R}{\partial b_R}.

For Doppler:

.. math::

   \frac{\partial D}{\partial r_R}
   =
   -\frac{1}{\lambda}
   \left(
     \frac{\Delta v}{\rho}
     -
     \frac{\dot{\rho}\Delta r}{\rho^2}
   \right)^\mathsf{T},

.. math::

   \frac{\partial D}{\partial v_R}
   =
   -\frac{1}{\lambda}
   \hat{u}^\mathsf{T},
   \qquad
   \frac{\partial D}{\partial d_R}
   =
   -\frac{c}{\lambda}
   \frac{\partial \dot{\Delta t}_R}{\partial d_R}.

For carrier phase:

.. math::

   \frac{\partial \Phi}{\partial r_R}
   =
   \frac{1}{\lambda}\hat{u}^\mathsf{T},
   \qquad
   \frac{\partial \Phi}{\partial b_R}
   =
   \frac{c}{\lambda}
   \frac{\partial \Delta t_R}{\partial b_R},
   \qquad
   \frac{\partial \Phi_\mathrm{meas}}{\partial N}
   =
   1.

Derivatives of transmitter ephemeris interpolation, light-time coupling,
Shapiro delay, plasma delay, and C/N0-based noise with respect to receiver
state are not included in the current analytic Jacobian.

Online and Precompute Equivalence
-------------------------------------------------------------------

The online path computes :math:`\mathcal{Y}(t_i, x_i) =
\operatorname{Compute}(t_i, x_i)` one epoch at a time.  The precompute path
computes the same measurement epochs for vectors ``receive_times`` and
``user_states``.  When a batch ionosphere/plasma provider is configured,
precompute first builds all light-time-corrected channels without running the
online plasma model, then calls the batch provider to fill
:math:`\Delta \rho_{\mathrm{plasma},ij}` for each epoch/channel before
evaluating the observables.  This keeps visibility, light-time, transmitter
clock terms, and batch plasma simulation ordered explicitly, and avoids an
online raytrace that would be overwritten by the batch result.

In the staged Lunar GNSS ODTS scenario, this contract is implemented with an
explicit file boundary: ``LunarGnssOdtsApp::Precompute`` writes all
light-time-corrected links with zero plasma delay, ``precompute_delays.py``
fills the delay terms with GCPM/IRI ray tracing (see
:doc:`ionosphere_plasmasphere`), and ``LunarGnssOdtsApp::Step`` merges the
delay file back into the truth channels before generating pseudorange,
Doppler, and optional TDCP measurements. TDCP is represented as a carrier-range
difference in meters and is processed by the UDU stochastic-cloning filter
because its measurement model depends on both the current and previous receiver
state.

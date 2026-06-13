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
``SetSunPositionProvider`` or by the default Sun provider.  The legacy
``options.shapiro_mu`` and ``options.shapiro_body_position`` fields are not
used by the current implementation.

Ionosphere / Plasma Delay
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ionosphere/plasma correction is represented as a positive code delay in
meters:

.. math::

   \Delta \rho_\mathrm{plasma} \ge 0.

It is disabled by default.  When
``options.apply_ionosphere_plasma_delay`` is true, the delay is supplied by one
of:

* ``SetCustomIonospherePlasmaDelayModel`` for a custom online link-by-link
  model.
* ``SetIonospherePlasmaRayTraceOptions`` for built-in link-by-link evaluation
  with the plasma ray tracer in ``environment/plasma/tec/raytrace``.
* ``SetBatchCustomIonospherePlasmaDelayModel`` for a custom precompute model
  over all epochs and all channels.
* ``options.default_ionosphere_plasma_delay_m`` when no provider is installed.

The built-in raytrace option calls ``pecsim::trace_ray`` directly.  The caller
must supply the ``pecsim::RayTraceConfig`` and, when desired, the IRI model
selection; ``GNSSMeasurements`` does not choose Kp, step size, integrator,
correction method, cutoff radius, or GCPM/IRI settings.  The GNSS channel
frequency is copied into ``RayTraceConfig::freq_Hz`` by default, but this can
be disabled through ``use_channel_frequency``.  Receiver and transmitter
positions are converted from ``options.frame`` to
``GnssIonospherePlasmaRayTraceOptions::raytrace_frame`` before tracing, and
the receive epoch is converted from ``options.receive_time_scale`` to
``raytrace_epoch_scale``.

The raytracer returns a path profile.  The GNSS measurement uses
``PathProfile::tec_delay_m`` by default:

.. math::

   \Delta \rho_\mathrm{plasma}
      = \frac{40.3}{f^2}\int_\Gamma n_e \, ds .

When ``delay_mode`` is ``TEC_PLUS_HIGHER_ORDER``, the second- and third-order
terms reported by the raytracer are added to this code delay.

The sign convention is dispersive:

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

Noise and Covariance
-------------------------------------------------------------------

For each selected channel, the measurement covariance is diagonal in the
current implementation:

.. math::

   R =
   \operatorname{diag}
   \left(
     \sigma_P^2,
     \sigma_D^2,
     \sigma_\Phi^2
   \right),

using the subset and order requested by ``options.observables``.

When C/N0 is finite, the tracking-loop noise models provide:

.. math::

   \sigma_P \quad [\mathrm{m}],
   \qquad
   \sigma_D \quad [\mathrm{Hz}],
   \qquad
   \sigma_\Phi \quad [\mathrm{cycles}].

If a channel does not have a finite positive standard deviation for an
observable, the filter wrapper uses a fallback standard deviation of ``1`` in
that observable's native unit.

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

The online path computes

.. math::

   \mathcal{Y}(t_i, x_i)
   =
   \operatorname{Compute}(t_i, x_i)

one epoch at a time.

The precompute path computes the same measurement epochs for vectors
``receive_times`` and ``user_states``:

.. math::

   \left\{
     \mathcal{Y}(t_i, x_i)
   \right\}_{i=0}^{M-1}.

When a batch ionosphere/plasma provider is configured, precompute first builds
all light-time-corrected channels without running the online plasma model:

.. math::

   \mathcal{C}_{ij}
   =
   \operatorname{BuildChannel}(t_i, x_i, \mathrm{PRN}_j),

then calls the batch provider to fill
:math:`\Delta \rho_{\mathrm{plasma},ij}` for each epoch/channel before
evaluating the observables.  This keeps visibility, light-time, transmitter
clock terms, and batch plasma simulation ordered explicitly, and it avoids
doing an online raytrace that would be overwritten by the batch result.

In the staged Lunar GNSS ODTS scenario, this contract is implemented with an
explicit file boundary: ``LunarGnssODTSSimulation::Precompute`` writes all
light-time-corrected links with zero plasma delay, ``precompute_delays.py``
fills the delay terms with GCPM/IRI ray tracing, and
``LunarGnssODTSSimulation::Run`` merges the delay file back into the truth
channels before generating pseudorange, Doppler, and optional TDCP
measurements. TDCP is represented as a carrier-range difference in meters and
is processed by the UDU stochastic-cloning filter because its measurement model
depends on both the current and previous receiver state.

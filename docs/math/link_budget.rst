Link Budget and Visibility
===================================================================

Purpose
-------------------------------------------------------------------

This specification defines the mathematical contract for LuPNT's GNSS
line-of-sight visibility test and its carrier-to-noise-density
(:math:`C/N_0`) link budget.  Together these decide *which* GNSS channels a
lunar receiver can track and *how noisy* each observable is.  The equations
here feed the C/N0 threshold and the tracking-loop noise
(:math:`\sigma_P, \sigma_D, \sigma_\Phi`) referenced by
``docs/math/gnss_measurements.rst``.  The concrete code lives in
``cpp/lupnt/measurements/comms_utils.{h,cc}`` (path loss, tracking-loop noise,
visibility), ``cpp/lupnt/measurements/antenna.{h,cc}`` (gain patterns),
``cpp/lupnt/measurements/gnss_measurement.cc`` (the ``LinkBudget`` /
``ComputeCN0`` glue), and ``cpp/lupnt/agents/gnss_attitude.*`` /
``gnss_yaw_steering.*`` (transmitter boresight).  The reference model is the
author's PhD thesis, Chapter 6.4.3 (link budget) and Chapter 3.3 (EIRP); the
corresponding equations are cited inline.

Visibility Geometry
-------------------------------------------------------------------

Line-of-sight visibility between two points :math:`r_1` and :math:`r_2` in a
common frame, in the presence of one spherical occulting body of radius
:math:`R_B` centered at :math:`r_B`, is evaluated by
``ComputeVisibility`` (``cpp/lupnt/measurements/comms_utils.cc ::
ComputeVisibility``).  It is invoked per transmitter channel from
``GNSSMeasurements::BuildChannels`` when ``options.apply_visibility`` is set,
once per configured ``GnssOccludingBody`` (typically the Moon and the Earth).

.. code-block:: cpp

   bool ComputeVisibility(const Vec3& r1, const Vec3& r2, Real R_body,
                          const Vec3& r_body = Vec3::Zero(),
                          const Real min_alt = 10e3,
                          const Real min_elev_deg = 5.0);

The test distinguishes two regimes based on whether either endpoint is within
:math:`R_B + h_\mathrm{min}` of the body center
(:math:`h_\mathrm{min} =` ``min_alt``).

**Elevation-mask regime (near-surface endpoint).**  If :math:`r_1` is a
surface point, the local elevation of the link toward :math:`r_2` is

.. math::

   e_1
   =
   \cos^{-1}
   \!\left(
     \frac{(r_2 - r_1)^\mathsf{T}(r_B - r_1)}
          {\lVert r_2 - r_1 \rVert \, \lVert r_B - r_1 \rVert}
   \right)
   - \frac{\pi}{2},

and the link is rejected when :math:`e_1 < e_\mathrm{min}` with
:math:`e_\mathrm{min} =` ``min_elev_deg`` (default :math:`5^\circ`).  The same
test is applied symmetrically if :math:`r_2` is a surface point.

**Horizon-occlusion regime (both endpoints elevated).**  For a space-to-space
link the body is treated as a hard sphere.  With :math:`r = r_2 - r_1`,
:math:`r_{1B} = r_1 - r_B`,

.. math::

   \theta_1 = \cos^{-1}\!\left(\frac{-r_{1B}^\mathsf{T} r}
                                    {\lVert r \rVert\,\lVert r_{1B}\rVert}\right),
   \qquad
   \theta_2 = \sin^{-1}\!\left(\frac{R_B}{\lVert r_{1B}\rVert}\right),
   \qquad
   d_\mathrm{hor} = \sqrt{\lVert r_{1B}\rVert^2 - R_B^2}.

The line of sight is blocked (returns ``false``) when the transmitter lies
behind the body's limb *and* beyond the tangent point:

.. math::

   \theta_1 < \theta_2
   \quad\text{and}\quad
   \lVert r \rVert > d_\mathrm{hor}.

.. note::

   ``ComputeVisibility`` handles a **single** occluding body per call;
   ``BuildChannels`` loops over all configured ``GnssOccludingBody`` entries
   and requires every one to pass.  The GNSS spec's near-surface threshold of
   :math:`e > -10^\circ` is the scenario-level ``min_elev_deg`` passed by the
   caller; the code default is :math:`+5^\circ`.

C/N0 Link-Budget Equation
-------------------------------------------------------------------

The received carrier-to-noise-density ratio for one channel is assembled by
the anonymous-namespace helper ``LinkBudget`` and driven by
``GNSSMeasurements::ComputeCN0``
(``cpp/lupnt/measurements/gnss_measurement.cc :: LinkBudget`` /
``ComputeCN0``):

.. code-block:: cpp

   Real LinkBudget(Real P_tx_dbw, Real G_tx_db, Real G_rx_db, Real range_m,
                   GnssFreq freq, const GnssReceiverParams& rx_params) {
     constexpr double k_B = 1.38e-23;                 // [J/K] Boltzmann
     Real L_fs = FreeSpacePathLoss(range_m, freq_Hz);
     double kb_db = 10.0 * std::log10(k_B);
     return P_tx_dbw + G_tx_db + G_rx_db
            - rx_params.L_atm - L_fs - rx_params.L_ad
            - rx_params.L_pol - kb_db - 10.0 * log10(rx_params.T_eff);
   }

In closed form, with all terms in decibels,

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
   - 10\log_{10}(T_\mathrm{eff}),

where :math:`P_\mathrm{tx}` is the transmit power [dBW],
:math:`G_\mathrm{tx}`/:math:`G_\mathrm{rx}` are the transmit/receive antenna
gains [dBi], :math:`L_\mathrm{fs}` is free-space path loss,
:math:`L_\mathrm{atm}` is atmospheric loss (``0`` for the vacuum lunar link),
:math:`L_\mathrm{ad}` is A/D-converter / implementation loss,
:math:`L_\mathrm{pol}` is polarization loss, :math:`k_B` is Boltzmann's
constant, and :math:`T_\mathrm{eff}` is the effective system noise
temperature.  These live on ``GnssReceiverParams``
(``cpp/lupnt/agents/gnss_constellation.h``) with defaults
:math:`L_\mathrm{ad}=0.6`, :math:`L_\mathrm{pol}=1.0`,
:math:`L_\mathrm{atm}=0.0` dB, :math:`T_\mathrm{eff}=167.98` K.

.. note::

   This is thesis Eq. (6.59),
   :math:`(C/N_0) = P_\mathrm{tx} + G_\mathrm{tx} + G_\mathrm{rx}
   - 20\log_{10}(4\pi d f / c) - L_\mathrm{pol}
   - 10\log_{10}(k_B T_\mathrm{sys}) - R_\mathrm{loss}`,
   with two bookkeeping differences: the code splits the single thesis noise
   term :math:`10\log_{10}(k_B T_\mathrm{sys})` into
   :math:`10\log_{10}(k_B) + 10\log_{10}(T_\mathrm{eff})` (identical value),
   and it carries an extra :math:`L_\mathrm{atm}` term (zero on the lunar
   link) and names the implementation loss :math:`L_\mathrm{ad}` rather than
   :math:`R_\mathrm{loss}`.  Thesis Table 6.3 gives the reference parameter
   values (GPS :math:`P_\mathrm{tx}\!\approx\!16\text{--}19` dBW,
   :math:`G_\mathrm{rx}=14` dBi, :math:`T_\mathrm{sys}=162` K,
   :math:`L_\mathrm{pol}=1.0` dB, :math:`R_\mathrm{loss}=0.9` dB).

Free-Space Path Loss
-------------------------------------------------------------------

Free-space path loss is ``FreeSpacePathLoss``
(``cpp/lupnt/measurements/comms_utils.cc``):

.. code-block:: cpp

   Real FreeSpacePathLoss(Real dist, Real freq) {
     return 20.0 * log10((4.0 * PI * dist) / (C / freq));
   }

i.e.

.. math::

   L_\mathrm{fs}
   =
   20\log_{10}\!\left(\frac{4\pi d f}{c}\right)
   =
   20\log_{10}\!\left(\frac{4\pi d}{\lambda}\right)
   \quad[\mathrm{dB}],

with :math:`d` the transmitter-receiver slant range, :math:`f` the carrier
frequency (from ``GNSS_FREQ_MAP``), and :math:`\lambda = c/f`.  This is thesis
Eq. (3.18) / the loss term of Eq. (6.59).

EIRP
-------------------------------------------------------------------

LuPNT carries the transmit power and transmit antenna gain **separately**:
``GnssConstellation::GetTransmitPowerDbw`` supplies :math:`P_\mathrm{tx}` and
``GnssConstellation::GetTransmitterAntenna`` supplies the pattern that yields
:math:`G_\mathrm{tx}(\theta_\mathrm{tx}, \varphi_\mathrm{tx})`.  The
equivalent isotropically radiated power is the sum

.. math::

   \mathrm{EIRP}(\theta_\mathrm{tx}, \varphi_\mathrm{tx})
   =
   P_\mathrm{tx}
   + G_\mathrm{tx}(\theta_\mathrm{tx}, \varphi_\mathrm{tx}).

For constellations that publish EIRP directly rather than
power-plus-pattern (e.g. Galileo), :math:`P_\mathrm{tx}` is loaded as the EIRP
and the transmit pattern is set to its peak-referenced shape.

The thesis instead poses the *required* EIRP as a design output
(Chapter 3.3, Eq. (3.17)) for a surface user at a fixed received-power floor
:math:`P_\mathrm{rx}`:

.. math::

   \mathrm{EIRP}
   =
   P_\mathrm{rx}
   + L_\mathrm{fs}
   - G_\mathrm{rx}
   + L_\mathrm{tx}
   + L_\mathrm{rx},

with :math:`L_\mathrm{tx}, L_\mathrm{rx}` the transmitter/receiver cable
losses (thesis Table 3.6, :math:`f = 2492.028` MHz,
:math:`P_\mathrm{rx}=-160/-147` dBW, :math:`G_\mathrm{rx}=3.0` dBi).  This is
the inverse of the link-budget code above: LuPNT propagates a *known* EIRP
forward to :math:`C/N_0`; the thesis solves for the EIRP that meets a
:math:`C/N_0` (or received-power) requirement.

Antenna Gain Patterns
-------------------------------------------------------------------

The gains :math:`G_\mathrm{tx}` and :math:`G_\mathrm{rx}` are table lookups
into measured patterns, ``Antenna::ComputeGain``
(``cpp/lupnt/measurements/antenna.cc :: ComputeGain``).  A pattern is a
function of the off-boresight (elevation) angle :math:`\varphi` and azimuth
:math:`\theta` about the boresight:

.. math::

   G(\theta, \varphi)
   =
   \operatorname{interp}\big(\varphi, \theta;\ \texttt{gain\_}\big)
   \quad[\mathrm{dB}].

.. code-block:: cpp

   Real Antenna::ComputeGain(Real theta, Real phi) const {
     if (n_dim_ == 0) return 0.0;                  // omni
     theta = WrapToTwoPi(theta) * DEG;
     phi   = WrapToPi(phi) * DEG;
     if (phi > phi_max_ || phi < -phi_max_) return NAN;   // outside coverage
     if (phi < 0 && phi_(0) >= 0) phi = -phi;             // mirror symmetry
     if (n_dim_ == 2) gain = LinearInterp2d(phi_, theta_, gain_, phi, theta);
     else             gain = LinearInterp1d(phi_, gain_, phi);
     return gain;
   }

Patterns are read from ANTEX/ACE data files by ``Antenna::LoadAntennaPattern``
(1D :math:`G(\varphi)` or 2D :math:`G(\varphi,\theta)` tables), normalized to
:math:`\varphi \in [-180,180]^\circ`, :math:`\theta \in [0,360]^\circ` by
``FormatAntennaPattern``.  An empty name gives an omni antenna
(:math:`G \equiv 0`), and directions beyond the pattern's ``phi_max_`` return
``NaN`` (used as an implicit sidelobe cutoff).  A closed-form parabolic
alternative is available for quick studies,
``cpp/lupnt/measurements/comms_utils.cc :: ParabolicAntennaGain``:

.. math::

   G(\varphi)
   =
   G_\mathrm{max}
   - 12\left(\frac{\varphi}{\varphi_\mathrm{3dB}}\right)^2.

Transmitter Attitude and Boresight Angles
-------------------------------------------------------------------

Evaluating :math:`G_\mathrm{tx}` requires the transmitter body frame, because
the off-boresight angles :math:`(\theta_\mathrm{tx}, \varphi_\mathrm{tx})` are
defined in the GNSS satellite's yaw-steering attitude.
``GNSSMeasurements::ComputeCN0`` builds that frame with
``GnssAttitude::Compute`` (``cpp/lupnt/agents/gnss_attitude.cc``) at the
transmit epoch, then projects the line of sight onto it:

.. code-block:: cpp

   GnssAttitude::Compute(r_tx_gcrf, v_tx_gcrf, r_sun_gcrf, ex, ey, ez);
   Real phi_tx   = safe_acos(u_tx2rx_gcrf.dot(ez));               // off-boresight
   Real theta_tx = atan2(u_tx2rx_gcrf.dot(ey), u_tx2rx_gcrf.dot(ex));
   Real phi_rx   = safe_acos(u_rx2body.dot(u_rx2tx));             // receiver off-boresight
   Real G_tx = constellation->GetTransmitterAntenna(prn, freq)
                   .ComputeGain(theta_tx, phi_tx);
   Real G_rx = rx_antenna_.ComputeGain(0.0, phi_rx);

The orthonormal transmitter triad :math:`(e_x, e_y, e_z)` is the nominal
"yaw-steering" body frame:

.. math::

   e_z = \widehat{-r_\mathrm{tx}}
   \quad(\text{nadir / boresight}),
   \qquad
   e_y = \widehat{e_z \times \widehat{(r_\odot - r_\mathrm{tx})}}
   \quad(\text{Sun side}),
   \qquad
   e_x = \widehat{e_y \times e_z}.

The boresight/off-boresight and azimuth angles of the transmitter-to-receiver
unit vector :math:`\hat{u}` are

.. math::

   \varphi_\mathrm{tx} = \cos^{-1}\!\big(\hat{u}^\mathsf{T} e_z\big),
   \qquad
   \theta_\mathrm{tx}
   = \operatorname{atan2}\!\big(\hat{u}^\mathsf{T} e_y,\ \hat{u}^\mathsf{T} e_x\big).

For the receiver, ``BoresightTarget(t)`` supplies the receive-antenna
boresight aim point (default the frame origin; a provider can point it at the
Earth, per the thesis Earth-pointing receiver), and
:math:`\varphi_\mathrm{rx}` is the angle between the boresight direction and
the receiver-to-transmitter line of sight.

Yaw-Steering Attitude Law
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The velocity-aware overload ``GnssAttitude::Compute(r, v, r_sun, ...)`` builds
the same Sun-pointing frame through the documented **nominal yaw-steering
law** (Kouba 2009; Cheng et al. 2025, Eq. 1), implemented in
``cpp/lupnt/agents/gnss_yaw_steering.cc``.  The satellite-Sun-Earth geometry
is parameterized by the Sun elevation above the orbital plane :math:`\beta`
and the orbit angle :math:`\mu` from the midnight point:

.. code-block:: cpp

   Real GnssYawSteering::BetaAngle(...) {
     Vec3 n = r_sat.cross(v_sat).normalized();     // orbit normal
     return safe_asin(n.dot(r_sun.normalized()));
   }
   Real GnssYawSteering::NominalYawAngle(Real beta, Real mu) {
     return atan2(-tan(beta), sin(mu));            // Eq. (1)
   }

.. math::

   \beta = \sin^{-1}\!\big(\hat{n}^\mathsf{T}\hat{s}\big),
   \quad
   \hat{n} = \frac{r_\mathrm{tx}\times v_\mathrm{tx}}
                  {\lVert r_\mathrm{tx}\times v_\mathrm{tx}\rVert},
   \qquad
   \psi(\beta,\mu) = \operatorname{atan2}\!\big(-\tan\beta,\ \sin\mu\big).

``ComputeFromYawAngle`` then rotates the orbital reference axes
(:math:`e_x^\mathrm{ref}=\hat{v}`, :math:`e_y^\mathrm{ref}=e_z\times e_x^\mathrm{ref}`)
about the nadir axis by :math:`\psi`,

.. math::

   e_x = \cos\psi\, e_x^\mathrm{ref} + \sin\psi\, e_y^\mathrm{ref},
   \qquad
   e_y = -\sin\psi\, e_x^\mathrm{ref} + \cos\psi\, e_y^\mathrm{ref},

which is numerically identical to the direct Sun-vector construction for the
nominal case.  The class also provides block-specific *modeled* yaw laws for
eclipse/noon/midnight turns (GPS IIF/IIR/III, Galileo IOV/FOC, BDS-3
CAST/SECM; Eqs. 3-16), used to replace :math:`\psi` during maneuver windows.

.. note::

   The thesis (Chapter 6.4.3) assumes **nominal yaw steering** with the
   antenna boresight aligned to the spacecraft-Earth direction and the solar
   panel axis orthogonal to the Earth-satellite-Sun plane; off-nominal
   eclipse-season steering is noted as future work.  LuPNT implements the
   nominal law used by ``ComputeCN0`` and additionally exposes the maneuver
   laws as building blocks, but no maneuver-window state machine is wired into
   the C/N0 path yet.

Receiver Tracking-Loop Noise from C/N0
-------------------------------------------------------------------

Once :math:`C/N_0` is known, the per-observable measurement standard
deviations follow from the tracking-loop thermal-noise models.  Let
:math:`\gamma = 10^{(C/N_0)_\mathrm{dB\text{-}Hz}/10}` be the linear
carrier-to-noise ratio (``DecibelToDecimal``).

**Pseudorange (DLL).**  ``GNSSMeasurements::ComputeSigmaRange`` scales the
dimensionless code-tracking jitter ``SigmaDll``
(``cpp/lupnt/measurements/comms_utils.cc``) by the chip length in meters
:math:`c\,T_c`:

.. code-block:: cpp

   return SigmaDll(params, CN0_w) * (Real(C) * Tc);   // Tc = 1 / Rc

For the wide-band case (:math:`d \le 1/(T_c B_\mathrm{fe})`) the DLL variance
is

.. math::

   \sigma_P^2
   =
   (c\,T_c)^2
   \frac{B_n}{2\gamma}
   \frac{1}{B_\mathrm{fe} T_c}
   \left(1 + \frac{1}{T_n\,\gamma}\right),

matching thesis Eq. (6.60) (with :math:`T_n` the DLL integration time,
:math:`B_n` the DLL bandwidth, :math:`d` the early-late spacing,
:math:`B_\mathrm{fe} = b\,R_c`).  ``SigmaDll`` also carries the intermediate-
and wide-correlator cases.

**Carrier phase (PLL).**  ``GNSSMeasurements::ComputeSigmaCarrierPhase``
converts the PLL phase jitter to length via :math:`\lambda/2\pi`:

.. math::

   \sigma_\Phi
   =
   \frac{\lambda}{2\pi}\,\sigma_\mathrm{PLL},
   \qquad
   \sigma_\mathrm{PLL}^2
   =
   \frac{B_p}{\gamma}
   \left(1 + \frac{1}{2\,T_p\,\gamma}\right).

.. note::

   The thesis Eq. (6.61) writes
   :math:`\sigma_\phi^2 = (\lambda/2\pi)^2\,[\,B_p/(2\gamma)\,(1 + 1/(2 T_p
   \gamma))\,]`.  The LuPNT ``SigmaPll`` uses :math:`B_p/\gamma` (no factor of
   two in the denominator), so the code's carrier-phase variance is larger
   than the thesis expression by a factor of two.  The DLL term matches the
   thesis (both use :math:`B_n/(2\gamma)`).

**Doppler / range-rate (FLL).**
``GNSSMeasurements::ComputeSigmaRangeRate`` uses a frequency-lock-loop model:

.. math::

   \sigma_D
   =
   \frac{\lambda}{2\pi\,T}
   \sqrt{\frac{4 B_f}{\gamma}\left(1 + \frac{1}{T\,\gamma}\right)},

with :math:`B_f` the FLL bandwidth and :math:`T` the integration time.

Model Boundaries
-------------------------------------------------------------------

* Visibility uses spherical occulting bodies and a per-body pass/fail; no
  atmospheric refraction or terrain horizon beyond the elevation mask.
* The C/N0 budget is per-channel and static in loss terms
  (:math:`L_\mathrm{atm}, L_\mathrm{ad}, L_\mathrm{pol}, T_\mathrm{eff}`);
  only :math:`L_\mathrm{fs}` and the two antenna gains vary with geometry.
* The transmit gain uses the **nominal** yaw-steering attitude; maneuver-law
  attitudes are available in ``GnssYawSteering`` but not yet applied in
  ``ComputeCN0``.
* The PLL noise coefficient differs from the thesis by a factor of two (see
  note above); the DLL and FLL forms match their reference expressions.

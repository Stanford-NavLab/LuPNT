Ephemeris and Almanac Design
===================================================================

Purpose
-------------------------------------------------------------------

This specification defines the mathematical contract for LuPNT's broadcast
**navigation-message** models for lunar satellites: the precise,
short-validity ``CartesianEphemeris`` and the coarse, long-validity
``Almanac``, both in ``cpp/lupnt/applications/ephemeris/``.  It fixes the state
representation, the fitting objective, the Chebyshev/Fourier bases and their
evaluation, the almanac element set with its argument-of-latitude
correction, and the bit-budget/quantization used to size a message against
the LunaNet budget.  Algorithms follow the author's PhD thesis Chapters 7
(ephemeris) and 8 (almanac) and Iiyama & Gao, "Ephemeris and Almanac Design
for Lunar Navigation Satellites"; equations are cited inline, and each is
tied to the implementing function.

State, Frame, and Unit Contract
-------------------------------------------------------------------

Both models represent a sampled Cartesian trajectory
:math:`rv \in \mathbb{R}^{N\times 6}`, row :math:`i` being
:math:`[r\ (\text{m});\ v\ (\text{m/s})]` at epoch :math:`t_s(i)` (seconds,
strictly increasing, arbitrary fixed origin).  ``Fit`` returns a flat
parameter vector; ``Eval`` reconstructs :math:`rv` at query epochs;
``EvalError`` returns RTN-decomposed RMS/95th-percentile statistics
(``EphemerisFitErrorStats``, ``applications/ephemeris/ephemeris_basis.h``).

Three frames are used (thesis section 7.1.1).  The **MCI** frame is inertial
(MOON_CI).  The **PA** frame (MOON_PA) rotates with the Moon's principal
axes at rate :math:`\omega_b = \Omega_\mathrm{Moon}`.  The **PAI**
(Principal-Axis Inertial) frame is the instantaneous PA orientation frozen at
each epoch -- a quasi-inertial frame used only as an intermediate for
osculating-element fitting.  When ``options.frame == Frame::MOON_PA`` the
code fits in PAI and returns in PA: the caller supplies MOON_PA states, the
fit adds the :math:`\omega\times r` offset to recover PAI velocity, and
``Eval`` subtracts it back.

.. code-block:: cpp

   // lunanet_ephemeris.cc :: FrameSpinRate / SpinCross
   double FrameSpinRate(Frame f) { return f == Frame::MOON_PA ? OMEGA_MOON : 0.0; }
   Vec3d  SpinCross(double spin, const Vec3d& r) { return spin * Vec3d(-r(1), r(0), 0.0); }

.. note::

   The frame handling is what lets a modest polynomial stay accurate: in
   MOON_PA the ascending node drifts at :math:`-\omega_b`, so the fitted
   two-body baseline's node tracks the Moon's rotation
   (``EphemerisGenApp`` broadcasts in ``output_frame = MOON_PA``,
   ``ephemeris_gen_app.cc``).

Ephemeris Representation
-------------------------------------------------------------------

``CartesianEphemeris`` (thesis section 7.2, Algorithm 3) models position as an
osculating two-body Kepler baseline plus a Chebyshev residual and an
optional Fourier residual, with velocity from the analytic time
derivatives.  The eight base parameters are the reference epoch, fit span,
and the osculating elements fixed at the window midpoint
:math:`t_0`:

.. math::

   \big[t_\mathrm{ref},\ t_\mathrm{fit},\ a,\ e,\ i,\ \Omega,\ \omega,\ M_\mathrm{ref}\big].

The **Chebyshev basis** :math:`T_j` (first kind, constant term first) is
evaluated in normalized time :math:`z = 2t_k/t_\mathrm{fit}` by the standard
recurrence (thesis Eqs. 7.1-7.2), with an analytic derivative basis for
velocity (``ephemeris_basis.h``):

.. math::

   T_0 = 1,\quad T_1 = z,\quad T_{j} = 2zT_{j-1} - T_{j-2},
   \qquad z = 2t_k/t_\mathrm{fit}\in[-1,1].

.. code-block:: cpp

   // ephemeris_basis.h :: ChebyshevBasis
   T.col(1) = z;
   for (int j = 2; j <= order; ++j)
     T.col(j) = 2.0 * z.array() * T.col(j-1).array() - T.col(j-2).array();

The **baseline** (``EvalKeplerianBaseline``) propagates the osculating orbit
via ``ClassicalToCart`` and, for a rotating frame, rigidly rotates the PAI
state by :math:`R_z(-\omega_b t_k)` and adds the corresponding
:math:`\dot R_z\, r` velocity term (Algorithm 3 steps 5-7).  The **Fourier
residual** uses harmonics of the argument of latitude :math:`u`
(``FourierBasis``):

.. math::

   \Delta\xi^{\mathrm{res}}(t) = \sum_{j=0}^{n} C^\xi_{F,j}\,T_j\!\Big(\tfrac{2t_k}{t_\mathrm{fit}}\Big)
     + \sum_{h=1}^{M}\big[C^\xi_{c,h}\cos(h u) + C^\xi_{s,h}\sin(h u)\big],
   \quad \xi\in\{x,y,z\},

and ``Eval`` reconstructs
:math:`\xi^{\mathrm{PA}} = \xi^{\mathrm{PA}}_{\mathrm{OE}} + \Delta\xi^{\mathrm{res}}`,
:math:`\dot\xi^{\mathrm{PA}} = \dot\xi^{\mathrm{PA}}_{\mathrm{OE}} + \Delta\dot\xi^{\mathrm{res}}`
(thesis Algorithm 3 steps 8-9).  The full parameter layout (``ParamNames``)
is the 8 base parameters, then per axis :math:`x,y,z` the
:math:`(\text{order}+1)` Chebyshev coefficients, then per axis the
:math:`2M` Fourier coefficients.

.. note::

   **Implementation vs. thesis.**  Thesis Eq. 7.3 shows a single
   second-harmonic Fourier term :math:`C_c\cos 2u + C_s\sin 2u`.  The code
   uses fundamental-and-multiples :math:`\cos(hu),\sin(hu)`,
   :math:`h=1..M`; with the default ``num_fourier_terms = 0`` the model is
   pure Chebyshev + Kepler.  Fourier terms require
   ``use_keplerian_baseline`` (the argument of latitude comes from the
   baseline orbit).

Ephemeris Least-Squares Fitting
-------------------------------------------------------------------

``CartesianEphemeris::Fit`` (thesis section 7.3.1, Eqs. 7.5-7.7) fixes the
osculating elements from the midpoint state (in PAI), subtracts the baseline
to form the position residual :math:`\Delta r`, and solves one joint linear
least-squares per axis over the stacked
:math:`[\,\text{Chebyshev}\mid\text{Fourier}\,]` design matrix:

.. math::

   \min_{\alpha}\ \sum_{m=1}^{N}\big\lVert \Delta r(t_m) - T_m\,\alpha \big\rVert_2^2,
   \qquad
   T = \big[\,\text{Chebyshev}\ \big|\ \text{Fourier}\,\big].

.. code-block:: cpp

   // lunanet_ephemeris.cc :: Fit  (per-axis column-pivoted Householder QR)
   const auto qr = design.colPivHouseholderQr();
   for (int d = 0; d < 3; ++d) {
     const VecXd sol = qr.solve(residual_pos.col(d));
     cheb_coeffs.col(d) = sol.head(cheb_len);
     if (four_len > 0) four_coeffs.col(d) = sol.tail(four_len);
   }

.. note::

   Thesis Eq. 7.4 samples at Chebyshev-Lobatto nodes to suppress Runge
   oscillation; the LuPNT model fits whatever epochs the caller supplies.
   ``EphemerisGenApp::GenerateFromAgent`` samples the predicted arc
   **uniformly** across the validity window
   (``ephemeris_fit_samples = 121`` by default), so node placement is a
   caller responsibility, not enforced by ``Fit``.  The osculating elements
   and the polynomial coefficients are fit **sequentially** (elements first,
   then a linear solve for the residual), matching the thesis's fast/robust
   scheme rather than a joint nonlinear solve.

Almanac Representation
-------------------------------------------------------------------

``Almanac`` (thesis section 8.1, Eq. 8.1; Iiyama & Gao Algorithm 2) fits each
osculating element directly as a low-order polynomial plus an
element-specific Fourier term over a multi-day window.  For element
:math:`\xi(t)` in normalized time :math:`s = t_k/t_\mathrm{fit}`:

.. math::

   \xi(t) = \sum_{p=0}^{P}\beta^\xi_p\,s^p
     + \sum_{h=1}^{M}\big[C^\xi_{c,h}\cos(h\omega_\xi t_k) + C^\xi_{s,h}\sin(h\omega_\xi t_k)\big],

with the element-specific base frequency (``ElementFreq``): the orbital mean
motion for the semi-major axis, the doubled sidereal harmonic for the rest
(thesis Eq. 8.1):

.. math::

   \omega_a = \sqrt{GM/a_\mathrm{ref}^3} = \tfrac{2\pi}{T_\mathrm{orb}},
   \qquad
   \omega_{e,i,\Omega,M,u} = \tfrac{4\pi}{T_\mathrm{sid}},
   \quad T_\mathrm{sid} = 27.321661\ \text{day}.

.. code-block:: cpp

   // lunanet_almanac.cc :: ElementFreq
   if (d == 0) return std::sqrt(gm / (a_ref*a_ref*a_ref));   // a: mean motion
   return 4.0 * PI / sidereal_period_s;                      // e,i,node,M,u

Six elements are represented: :math:`a, e, i, \Omega`, the **mean-anomaly
residual** :math:`M`, and the **argument-of-latitude correction**
:math:`u`.  The mean anomaly is reconstructed as the nominal two-body drift
integrated from the fitted semi-major axis plus the fitted residual
(thesis Algorithm 4 steps 4-5):

.. math::

   M(t) = \int_0^{t_k}\! n(\tau)\,d\tau + \Delta M(t),
   \qquad
   n(\tau) = \sqrt{GM/a(\tau)^3},

computed by cumulative-trapezoid integration (``CumTrapz``,
``ephemeris_basis.h``).  The argument of periapsis is **replaced** by
:math:`u = \nu + u_\mathrm{corr}` so the reconstruction stays
well-conditioned near periapsis of eccentric orbits.

.. code-block:: cpp

   // lunanet_almanac.cc :: EvalSeries
   const VecXd n_series = (opt.gm * s.a.array().pow(-3)).sqrt();
   s.m = CumTrapz(n_series, t_k) + m_res;   // nominal drift + residual

The parameter layout (``ParamNames``) is
:math:`[t_\mathrm{ref}, t_\mathrm{fit}, a_\mathrm{ref}]` followed, for each of
:math:`a,e,i,\Omega,M,u`, by the :math:`(P+1)` polynomial and :math:`2M`
Fourier coefficients.  With the defaults ``poly_order = 1``,
``num_fourier_terms = 1`` each element has 4 coefficients
:math:`[\beta_0,\beta_1,C_c,C_s]` -- the 24 fitted parameters of thesis
section 8.1 (plus the 3 bookkeeping entries), versus 13 for a GPS almanac.

Almanac Least-Squares Fitting
-------------------------------------------------------------------

``Almanac::Fit`` (thesis section 8.2, Algorithm 4 fitting; Eqs. 8.2-8.6)
computes unwrapped osculating elements from the (PAI) states, forms one
design matrix per frequency family (``ElementBasis``), and solves each
element by QR least squares:

.. math::

   \bar a = \tfrac{1}{N}\sum_m a(t_m),
   \qquad
   \min_{\alpha_\xi}\ \sum_m \big\lVert \xi(t_m) - T_{\xi,m}\,\alpha_\xi\big\rVert_2^2 .

The node :math:`\Omega` is fit directly (it drifts secularly in PA).  The
mean-anomaly residual is fit against the nominal drift
:math:`M_\mathrm{nom}=\int n\,dt`, and the :math:`u` correction is fit to
:math:`\omega + \mathrm{wrap}(\nu_\mathrm{osc}-\nu_\mathrm{fit})` so the
target has no :math:`2\pi` jumps:

.. code-block:: cpp

   // lunanet_almanac.cc :: Fit  (mean-anomaly residual, then u-correction)
   const VecXd m_coeffs = qr_s.solve(coe.col(5) - m_nom);
   // ...
   u_target(i) = coe(i,4) + std::atan2(std::sin(nu_o - nu_f), std::cos(nu_o - nu_f));
   const VecXd u_coeffs = qr_s.solve(u_target);

``Eval`` reconstructs the Cartesian state by ``ClassicalToCart`` with
:math:`\mathrm{argp} = u_\mathrm{corr}` and the wrapped mean anomaly, then
removes the :math:`\omega\times r` offset for a rotating output frame
(thesis Algorithm 4 steps 6-10).

.. note::

   **Implementation vs. thesis.**  Thesis Algorithm 4 step 11 computes
   velocity by finite differencing the reconstructed position.  The code
   instead returns the **analytic two-body velocity** from
   ``ClassicalToCart`` (minus the frame :math:`\omega\times r` term),
   avoiding the finite-difference step size.

Message Size and Bit Budget
-------------------------------------------------------------------

The broadcast cost is sized per parameter so that quantization keeps the
position error under a tolerance (thesis section 7.3.2, Eqs. 7.8-7.12).  The
implementation lives in ``applications/ephemeris/ephemeris_app.cc`` (the
``EphemerisApp`` bit-budget helpers).
For each parameter, ``SearchParamBits`` binary-searches the minimum number
of fractional bits :math:`k` such that a :math:`2^{-k}` step changes the
**worst-case** position over the window by less than the tolerance
:math:`\delta_x`:

.. math::

   k_j = \min\Big\{k : \max_{m}\big\lVert
     f(\alpha + 2^{-k}e_j,\,t_m) - f(\alpha,\,t_m)\big\rVert_2 \le \delta_x\Big\}.

.. code-block:: cpp

   // ephemeris_app.cc :: SearchParamBits (worst-case over the window)
   const double eps = std::pow(2.0, -k);
   p(idx) += eps;  const MatXd plus  = eph.Eval(t_s, p);
   p(idx) -= 2*eps; const MatXd minus = eph.Eval(t_s, p);
   return std::max(err_plus, err_minus);   // max position offset over t_s

The total is the sum of per-parameter integer + fractional + sign + margin
bits (thesis Eqs. 7.11-7.12), ``RangeBits``:

.. math::

   b_j = \max\!\big(\lceil\log_2(\alpha_{j,\max}-\alpha_{j,\min})\rceil, 1\big)
     + k_j + b_s + b_m,
   \qquad
   \text{Message Size} = \sum_j b_j .

Angle parameters (:math:`i,\Omega,\omega,M,u`) use range :math:`2\pi`, signed,
no margin; eccentricity uses range 1, unsigned; the semi-major axis and
Chebyshev/Fourier coefficients use the observed magnitude range with a margin
bit (``RangeBits`` / ``IsAngleParam`` / ``IsEccentricityParam``).

.. code-block:: cpp

   // ephemeris_app.cc :: RangeBits
   const int bits_range = range > 0.0
       ? std::max(1, (int)std::ceil(std::log2(range))) : 1;
   return bits_range + std::max(bits_frac, 0)
          + (signed_range ? 1 : 0) + (margin_bit ? 1 : 0);

.. note::

   Thesis Eqs. 7.9-7.10 size bits from a closed-form
   :math:`k_j = -\lceil\log_2(\delta_x/|\partial f/\partial\alpha_j|)\rceil`
   plus :math:`\pm 1` refinement; LuPNT reaches the same target with a
   direct binary search over the realized worst-case quantized error.  The
   default tolerance is :math:`\delta_x = 0.01` m.  The LunaNet
   Interoperability Specification allots roughly **900 bits** to ephemeris
   data per broadcast frame (thesis section 1.2.3); the results (thesis
   Table 7.2) meet the 3-sigma SISE thresholds (13.43 m / 1.2 mm/s for LCRNS)
   within that budget for ELFO fits up to ~240 min, with Chebyshev +
   osculating elements (+ Fourier for :math:`\ge 240` min) as the best
   accuracy-per-bit representation.

Broadcast-Message Generation
-------------------------------------------------------------------

``EphemerisGenApp`` (``applications/ephemeris/ephemeris_gen_app.{h,cc}``) is the
``LunaNetSubApp`` that produces both message types from a satellite's own
predicted arc.  On each ``Step(t)`` it refreshes the precise ephemeris
(default 2-hour validity, ``ephemeris_window_s``) and the coarse almanac
(default 15-day validity, ``almanac_window_s``), samples the owning agent's
predicted state uniformly over the window (``Agent::GetStateAt``), converts
to ``output_frame`` (typically MOON_PA), fits the corresponding model, and
appends a ``BroadcastMessage`` holding the fitted parameter vector.
``LatestEphemeris(t)`` / ``LatestAlmanac(t)`` mirror a receiver selecting the
currently valid page.

Model Boundaries
-------------------------------------------------------------------

* Inputs to ``Fit``/``EvalError`` must already be in ``options.frame``;
  MOON_PA broadcasting requires the caller to ``ConvertFrame`` from MCI
  first.
* ``CartesianEphemeris`` Fourier terms require the Keplerian baseline;
  the almanac always carries the baseline (osculating elements are the
  representation).
* Node placement (Lobatto vs. uniform) is a caller responsibility; the
  models fit the supplied epochs directly.
* Almanac velocity is analytic two-body (not finite-differenced); ephemeris
  velocity is the analytic derivative of the Chebyshev/Fourier bases.
* The bit-budget/quantization analysis lives in ``EphemerisApp``, not
  in the model classes; ``EphemerisGenApp`` stores unquantized
  double-precision parameters.

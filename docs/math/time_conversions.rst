Time Systems and Conversions
===================================================================

Purpose
-------------------------------------------------------------------

This specification defines the mathematical contract for LuPNT's time
systems and the conversions between them.  It fixes the epoch origin, the
scale definitions, the exact offset/rate/periodic equations, the body and
reference-system each scale is proper to, and the leap-second and
Earth-orientation inputs.  Every conversion is tied to the function that
implements it in ``cpp/lupnt/conversions/time_conversions.{h,cc}``; the
Python surface is ``pnt.convert_time`` and friends from
``python/bindings/py_time_converter.cc``.  The relativistic content follows
the author's PhD thesis, Chapter 2.3 (Iiyama), and the Turyshev 2026 (ApJ
997:97) lunar-timescale formulation; equations are cited inline.

Epoch, Argument, and Unit Contract
-------------------------------------------------------------------

Every scalar time in the API is a **coordinate offset in seconds from the
J2000 epoch** of the corresponding scale, not a calendar date.  The J2000
origin is

.. math::

   \mathrm{MJD}_\mathrm{J2000} = 51544.5\ \text{(TT days)},
   \qquad
   \mathrm{JD}_\mathrm{J2000} = 2451545.0\ \text{(TT days)}.

The Julian/Modified-Julian bridges are exact linear maps
(``MjdToTime``/``TimeToMjd``/``JdToTime``/``TimeToJd``):

.. math::

   t = (\mathrm{MJD} - 51544.5)\cdot 86400,
   \qquad
   \mathrm{MJD} = t/86400 + 51544.5 .

.. code-block:: cpp

   // time_conversions.cc :: MjdToTime / TimeToMjd
   Real MjdToTime(Real mjd) { return (mjd - MJD_J2000_TT) * SECS_DAY; }
   Real TimeToMjd(Real t)   { return t / SECS_DAY + MJD_J2000_TT; }

The ``Time`` enum (``cpp/lupnt/core/constants.h``) enumerates the supported
scales:

.. math::

   \texttt{UT1},\ \texttt{UTC},\ \texttt{TAI},\ \texttt{TDB},\ \texttt{TT},\
   \texttt{TCG},\ \texttt{TCB},\ \texttt{GPS},\ \texttt{JD\_TT},\
   \texttt{JD\_TDB},\ \texttt{TCL},\ \texttt{LT}.

Conversion Graph
-------------------------------------------------------------------

``ConvertTime(t, from, to)`` (``time_conversions.cc :: ConvertTime``) routes
any pair through the canonical chain, with TAI/TT as the two hubs.  The
vectorized overload ``ConvertTime(VecX, from, to)`` maps element-wise,
with special handling for the position-dependent lunar scales.  The
Python bindings expose exactly this:

.. code-block:: python

   import pylupnt as pnt
   t_tt = pnt.convert_time(t_tai, pnt.Time.TAI, pnt.Time.TT)   # scalar
   t    = pnt.convert_time(t_vec, pnt.Time.GPS, pnt.Time.TDB)  # VecX

Proper vs. Coordinate Time
-------------------------------------------------------------------

Following the thesis (Ch. 2.3.1), a clock measures **proper time**
:math:`\tau` along its worldline, while a **coordinate time** :math:`t`
labels events in a 4-D relativistic reference system (BCRS, GCRS, LCRS) and
cannot be read off any single physical clock.  ``IsCoordinateTimeScale``
partitions the enum: TT, TCG, TDB, TCB, TCL, and LT are coordinate scales;
UT1, UTC, TAI, GPS are not.

.. code-block:: cpp

   // time_conversions.cc :: IsCoordinateTimeScale
   case Time::TT: case Time::TCG: case Time::TDB:
   case Time::TCB: case Time::TCL: case Time::LT: return true;

``ConvertCoordinateTime`` is the coordinate-only entry point; its
position-aware overloads take a BCRS position :math:`x_\mathrm{bcrs}` so
that the position-dependent terms of the TDB/TCB/TCL transforms
(below) can be evaluated at the event location rather than at the
body center.

Atomic and Terrestrial Scales
-------------------------------------------------------------------

**TAI -> TT** is a defining constant offset (thesis Eq. 2.38):

.. math::

   \mathrm{TT} = \mathrm{TAI} + 32.184\ \text{s} .

.. code-block:: cpp

   // time_conversions.cc :: TaiToTt / TtToTai   (TT_TAI_OFFSET = 32.184 s)
   Real TaiToTt(Real t_tai) { return t_tai + TT_TAI_OFFSET; }
   Real TtToTai(Real t_tt)  { return t_tt  - TT_TAI_OFFSET; }

**TAI -> GPS** is the fixed 19-second offset that removes the leap seconds
accumulated at the GPS epoch (thesis Eq. 2.39); GPS time runs without leap
seconds thereafter (``TaiToGps``/``GpsToTai``):

.. math::

   \mathrm{GPS} = \mathrm{TAI} - 19\ \text{s} .

.. note::

   Only GPS is implemented among the GNSS scales.  The thesis also lists
   GLONASS :math:`= \mathrm{TAI}-34` s (leap-tracking, offset by UTC),
   Galileo :math:`= \mathrm{TAI}-19` s, and BeiDou :math:`= \mathrm{TAI}-33`
   s (Eqs. 2.40-2.42); these are not in the ``Time`` enum.

Earth-Rotation and Leap-Second Scales (UTC, UT1)
-------------------------------------------------------------------

**TAI <-> UTC** applies the integer leap-second table
:math:`\Delta_\mathrm{AT}(\mathrm{MJD})` from
``interfaces/tai_utc.h`` (thesis Eq. 2.43):

.. math::

   \mathrm{UTC} = \mathrm{TAI} - \Delta_\mathrm{AT},
   \qquad
   \Delta_\mathrm{AT} = \texttt{GetTaiUtcDifference}(\mathrm{MJD}) .

**UTC <-> UT1** adds the observed Earth-orientation difference
:math:`(\mathrm{UT1}-\mathrm{UTC})` looked up from the EOP table
(``interfaces/eop.h``):

.. math::

   \mathrm{UT1} = \mathrm{UTC} + (\mathrm{UT1}-\mathrm{UTC}),
   \qquad
   (\mathrm{UT1}-\mathrm{UTC}) = \texttt{GetUt1UtcDifference}(\mathrm{MJD}) .

.. code-block:: cpp

   // time_conversions.cc :: UtcToUt1
   Real mjd_utc = t_utc / SECS_DAY + MJD_J2000_TT;
   Real ut1_utc = GetUt1UtcDifference(mjd_utc);
   return t_utc + ut1_utc;

UT1 in turn drives Earth rotation: ``EarthRotationAngle`` gives the ERA
:math:`\theta = 2\pi(0.7790572732640 + 1.00273781191135448\, t_\mathrm{UT1}/86400)`,
and ``GreenwichMeanSiderealTime``/``GreenwichApparentSiderealTime`` give GMST/GAST
(the latter adding the equation of the equinoxes).  These are the only
rotation quantities keyed on UT1.

Barycentric and Geocentric Coordinate Scales
-------------------------------------------------------------------

The four Einstein-scaled coordinate times share the reference epoch
:math:`T_0 =` 1977-01-01 00:00:00 TAI and the defining rates

.. math::

   L_B = 1.550519768\times10^{-8},\quad
   L_G = 6.969290134\times10^{-10},\quad
   \mathrm{TDB}_0 = -65.5\ \mu\text{s},

with :math:`\mathrm{MJD}_{T_0} = \texttt{MJD\_COORDINATE\_TT\_TCG\_TCB}`.

**TCG <-> TT** (geocentric rate, thesis Eq. 2.37) removes the geoid potential
secular rate:

.. math::

   \mathrm{TT} = \mathrm{TCG} - L_G\,(\mathrm{JD}_\mathrm{TCG}-\mathrm{JD}_0)\cdot 86400 .

.. code-block:: cpp

   // time_conversions.cc :: TtToTcg  (inverse of Eq. 2.37 in offset form)
   Real tt_tcg = -L_G / (1.0 - L_G) * (mjd_tt - MJD_COORDINATE_TT_TCG_TCB) * SECS_DAY;
   return t_tt - tt_tcg;

**TCB <-> TDB** (barycentric rate + offset, thesis Eq. 2.32):

.. math::

   \mathrm{TDB} = \mathrm{TCB}
     - L_B\,(\mathrm{JD}_\mathrm{TCB}-\mathrm{JD}_0)\cdot 86400 + \mathrm{TDB}_0 .

.. code-block:: cpp

   // time_conversions.cc :: TcbToTdb
   Real t_tdb = t_tcb - L_B * (mjd_tcb - MJD_COORDINATE_TT_TCG_TCB) * SECS_DAY + TDB_0;

**TT <-> TDB** has two implementations that intentionally differ:

* Forward ``TtToTdb`` uses the analytic Fairhead/Moyer series (thesis
  Eqs. 2.33-2.34), :math:`M_E = (357.53 + 0.9856003\,d_\mathrm{J2000})^\circ`,

  .. math::

     \mathrm{TDB} = \mathrm{TT} + 0.001658\sin M_E + 0.000014\sin 2M_E\ \text{s} .

* Inverse ``TDBToTt`` (no position) uses the NAIF single-term expansion
  :math:`\mathrm{TT} = \mathrm{TDB} - k\sin E`, with :math:`k=1.657\times10^{-3}`,
  eccentric anomaly :math:`E = M + e_b\sin M`,
  :math:`M = 6.239996 + 1.99096871\times10^{-7}\,t_\mathrm{TDB}`.

.. note::

   The scalar TT<->TDB pair is a periodic approximation good to
   :math:`\lesssim 2` ms (thesis: "deviates by a periodic oscillation of
   less than 2 milliseconds").  The **position-aware** overloads
   ``TtToTdb(t, x_bcrs)`` / ``TDBToTt(t, x_bcrs)`` instead evaluate the full
   GCRS<->BCRS 4-D transform (thesis Eq. 2.35, Turyshev Eq. 21/22): a
   trapezoidal integral of :math:`c^{-2}` and :math:`c^{-4}` integrands built
   from the Earth's barycentric velocity and the external Solar-System
   potential :math:`w_\mathrm{ext}=\sum_{B\neq E} GM_B/r_{EB}`, plus a
   position term :math:`-v_E\!\cdot\!r_E/c^2`.  ``TtToTdb(t, x_bcrs)`` inverts
   this by 10-step Newton iteration.  See
   ``TdbToTtMinusTdbEq21`` / ``IntegrateTdbToTtTerms``.

Lunicentric Coordinate Scales (TCL, LT)
-------------------------------------------------------------------

The Moon's analogues follow Turyshev 2026, with LCRS as the lunar 4-D
system.  **TCB -> TCL** is the barycentric-to-lunicentric 4-D transform
(thesis Eqs. 2.44-2.45; Turyshev Eq. 23/25), structurally identical to the
Earth case with the Moon substituted:

.. math::

   \mathrm{TCL} = \mathrm{TCB}
     - \frac{1}{c^2}\!\int \!\Big(\tfrac{1}{2}v_M^2 + w_\mathrm{ext}(x_M)\Big)d\mathrm{TCB}
     - \frac{1}{c^4}\!\int (\cdots)\,d\mathrm{TCB}
     - \frac{v_M\!\cdot r_M}{c^2} - \cdots

where :math:`r_M = x_\mathrm{bcrs}-x_M`, :math:`v_M = dx_M/d\mathrm{TCB}`, and
:math:`w_\mathrm{ext}(x_M)=\sum_{B\neq M} GM_B/|x_M-x_B|` sums the Sun,
Mercury, Venus, Earth, Mars, and the outer planets
(``kTclExternalBodies``).  Passing no position places the clock at the lunar
center (:math:`v_M\!\cdot r_M = 0`).  ``TclToTcb`` inverts by Newton
iteration.

.. code-block:: cpp

   // time_conversions.cc :: TcbToTcl  (Turyshev Eq. 23/25)
   auto [i2, i4] = IntegrateTcbToTclTerms(t_tcb);
   return t_tcb - i2/(C*C) - i4/(C*C*C*C) + TcbToTclPositionTerm(t_tcb, x_bcrs);

**TCL <-> LT** is the lunar analogue of TCG->TT, a defining rate on the
selenoid (thesis Eqs. 2.48-2.49) with
:math:`L_L = 3.13905\times10^{-11}`:

.. math::

   \mathrm{LT} = \mathrm{TCL} - L_L\,(\mathrm{TCL}-T_{L0}),
   \qquad
   \frac{d\,\mathrm{LT}}{d\,\mathrm{TCL}} = 1 - L_L .

.. code-block:: cpp

   // time_conversions.cc :: TclToLt  (T_L0 assumed = T_0)
   Real lt_tcl = -L_L * (mjd_tcl - MJD_COORDINATE_TT_TCG_TCB) * SECS_DAY;
   return t_tcl + lt_tcl;

.. note::

   :math:`L_L` is not yet internationally standardized; the code uses the
   Turyshev selenoid value :math:`3.13905\times10^{-11}` (potential
   :math:`\Phi_L = 2.82123744381\times10^{6}\ \mathrm{m^2/s^2}`), whereas
   Kopeikin adopts :math:`3.14027\times10^{-11}` (thesis Ch. 2.3.5).

Lunar Time as a Function of TT (direct path)
-------------------------------------------------------------------

Rather than chaining TT->TDB->TCB->TCL->LT, LuPNT implements the closed-form
:math:`\mathrm{TL}-\mathrm{TT}(\mathrm{TDB})` of Turyshev 2026 Eq. 57 (the
TDB-argument form of thesis Eq. 2.50) in ``TdbToLtMinusTt``.  The integral is
a trapezoidal sweep at :math:`dt = 0.01` day.  For repeated queries in a
window, ``InitLtMinusTtFit`` fits piecewise Chebyshev segments (degree 12,
1-day segments) over the window in a single forward sweep and caches them in
``s_lt_minus_tt_fit``; ``TdbToLtMinusTt`` then uses the fast path
(``ChebyshevFitModel::Eval``, ``numerics/cheby_fit.h``).  ``TdbToLt`` composes
:math:`\mathrm{TL} = \mathrm{TT}(\mathrm{TDB}) + (\mathrm{TL}-\mathrm{TT})`;
``LtToTdb``/``LtToTt`` invert by 10-step Newton iteration.

.. code-block:: cpp

   // time_conversions.cc :: TdbToLt / TtToLt
   Real TdbToLt(Real t_tdb) { return TDBToTt(t_tdb) + TdbToLtMinusTt(t_tdb); }
   Real TtToLt(Real t_tt)   { return TdbToLt(TtToTdb(t_tt)); }

A separate lunar-surface proper-time term is provided by
``GetProperTimeCorrectionTcl(t_tcg, x_mci)``, the first-order clock
correction :math:`-v_{LE}\cdot(r_x - r_{LE})/c^2` for a clock at MCI
position :math:`x_\mathrm{mci}` relative to the Moon's geocentric motion.

Conversion of Physical Quantities and Constants
-------------------------------------------------------------------

Because the scaled coordinate times (TT, TDB, LT) run at rates
:math:`(1-L_G)`, :math:`(1-L_B)`, :math:`(1-L_L)` relative to their unscaled
partners (TCG, TCB, TCL), spatial coordinates and gravitational parameters
must be rescaled to keep :math:`c` and the equations of motion invariant
(thesis Eqs. 2.52-2.57):

.. math::

   \mathcal{X}_\mathrm{TT} = (1-L_G)\,\mathcal{X}_\mathrm{TCG},\quad
   \mu_\mathrm{TT} = (1-L_G)\,\mu_\mathrm{TCG},

.. math::

   \mathcal{X}_\mathrm{TDB} = (1-L_B)\,\mathcal{X}_\mathrm{TCB},\quad
   \mathcal{X}_\mathrm{LT}  = (1-L_L)\,\mathcal{X}_\mathrm{TCL}.

LuPNT encodes these as the ``CoordinateScale`` factors
:math:`1-L_x` (``constants.h :: CoordinateScaleFactor``), with
``CoordinateScaleRatio(from, to)`` returning
:math:`(1-L_\mathrm{to})/(1-L_\mathrm{from})` for the rescale, guarded by
``AreCoordinateScalesConvertible`` (only barycentric<->barycentric,
geocentric<->geocentric, lunicentric<->lunicentric pairs share a constant
factor).

.. code-block:: cpp

   // constants.h :: CoordinateScaleFactor
   case CoordinateScale::TDB: return 1.0 - L_B;
   case CoordinateScale::TT:  return 1.0 - L_G;
   case CoordinateScale::TL:  return 1.0 - L_L;

Model Boundaries
-------------------------------------------------------------------

* All API times are seconds from the per-scale J2000 origin; calendar I/O
  goes through ``GregorianToTime`` / ``TimeToGregorianString``.
* Only GPS is implemented among GNSS scales; GLONASS/Galileo/BeiDou offsets
  are documented but not exposed.
* The scalar TT<->TDB and TDB->TT paths are periodic approximations
  (:math:`\lesssim 2` ms); sub-microsecond and position-dependent work must
  use the ``x_bcrs`` overloads via ``ConvertCoordinateTime``.
* TCL/TDB position-aware transforms integrate from :math:`T_0` per call
  unless a Chebyshev fit window is installed (LT only, via
  ``InitLtMinusFit``).
* :math:`L_L` (and hence LT) is provisional pending international
  standardization.

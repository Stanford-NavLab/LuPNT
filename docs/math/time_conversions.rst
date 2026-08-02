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
   \texttt{TCG},\ \texttt{TCB},\ \texttt{GPS},\ \texttt{TCL},\ \texttt{LT}.

A Julian Date is a *representation* of an instant rather than a time scale, so
it is not a member of the enum; use ``JdToTime``/``TimeToJd`` with the relevant
scale.

Conversion Graph
-------------------------------------------------------------------

The time-scale relationships live in ``Epoch`` (``conversions/epoch.cc``) as a
registered set of directed edges.  Each edge returns the **offset**
:math:`(\text{to}-\text{from})` in seconds, given the epoch's reading in the
``from`` scale.  A conversion finds the shortest registered route
(``numerics/graphs.h :: FindShortestPath``) and sums the offsets along it.

.. math::

   \Delta_{A\to B}(t) \;=\; \sum_{k} \delta_{s_k \to s_{k+1}}\!\left(t_k\right),
   \qquad s_0 = A,\; s_n = B .

Because only small offsets are summed, two large absolute epochs are never
differenced and the result carries full double precision.

Edge costs span many orders of magnitude, which is why the route matters:

.. list-table::
   :header-rows: 1

   * - Class
     - Edges
     - Cost
   * - constant
     - TAI↔TT, TAI↔GPS
     - arithmetic
   * - table
     - TAI↔UTC (leap seconds), UTC↔UT1 (EOP)
     - lookup
   * - linear
     - TT↔TCG, TDB↔TCB, TCL↔LT
     - rescaling by :math:`L\sim10^{-8}`
   * - model
     - TT↔TDB
     - Chebyshev fit / analytic series
   * - integral
     - TDB↔TCL
     - trapezoidal sweep from :math:`T_0`

TCL↔LT and TDB↔TCB are single linear edges and are reached without touching the
TDB↔TCL integral.  The search minimises hops, which coincides with minimum cost
for this edge set because the one expensive edge is also the only bridge to the
lunicentric scales.

Precision contract
~~~~~~~~~~~~~~~~~~~

An absolute epoch is a ``float64`` count of seconds from J2000, so at
present-day dates :math:`|t|\sim10^{9}` s one unit in the last place is

.. math::

   |t|\,2^{-52} \;\approx\; 2.45\times10^{-7}\ \mathrm{s} \;\approx\; 245\ \mathrm{ns}
   \;\approx\; 73\ \mathrm{m}\times c .

This bounds any API that *returns* an absolute epoch:

.. list-table::
   :header-rows: 1

   * - Entry point
     - Returns
     - Accuracy (TDB↔TT, 2020–2035)
   * - ``Epoch::To`` / ``TimeScaleOffset``
     - offset
     - :math:`\sim3\times10^{-4}` ns
   * - ``TtMinusTdb``, ``TdbMinusTcl``, ``TdbMinusLt``
     - offset
     - :math:`\sim4\times10^{-4}` ns
   * - ``ConvertTime``
     - absolute epoch
     - :math:`\sim38` ns rms, 115 ns peak

``ConvertTime(t, from, to)`` delegates to ``Epoch``; its floor is a property of
the return type, not of the model.  Work that multiplies a time by :math:`c`, or
that differences two epochs, must use ``Epoch`` or the offset accessors.

.. code-block:: python

   import pylupnt as pnt

   e     = pnt.Epoch.from_seconds(t_tai, pnt.Time.TAI)
   e_tdb = e.to(pnt.Time.TDB)                          # exact
   dt    = pnt.time_scale_offset(e, pnt.Time.TCL)      # offset, sub-ps

   t_tt  = pnt.convert_time(t_tai, pnt.Time.TAI, pnt.Time.TT)   # ~245 ns floor

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

TT <-> TDB model chain
~~~~~~~~~~~~~~~~~~~~~~~

``TtMinusTdb(t)`` is the single implementation of the TT–TDB relation; both
``TtToTdb`` and ``TDBToTt`` delegate to it, so the pair round-trips exactly
whichever model is active.  It selects, in order:

1. **Chebyshev fit of the DE440t TT–TDB ephemeris.**  ``InitTtMinusTdbFit``
   samples the ``de440t.bsp`` time-ephemeris segment as a clean offset and fits
   piecewise Chebyshev segments (16-day segments, degree 12).  Reproduces the
   kernel to :math:`\sim4\times10^{-4}` ns.

2. **DE440 Eq. (3) relativistic integral**, opt-in via
   ``SetTtTdbModel(TtTdbModel::DE440_INTEGRAL)``.  Trapezoidal sweep of the
   :math:`c^{-2}` and :math:`c^{-4}` integrands from :math:`T_0`; agrees with the
   DE440t ephemeris to :math:`\sim7` ns rms with a :math:`+0.14` ns/yr secular
   term (see :ref:`model-boundaries`).

3. **Auto-fit**: if no fit covers the requested epoch, one decade-wide window is
   built on demand (``SetTtTdbAutoFit``, default on).  This is what makes the
   *default* accuracy :math:`\sim4\times10^{-4}` ns.

4. **Analytic fallback**, used only when the DE440t segment is unavailable — the
   two-term IAU/IERS TN 36 series

   .. math::

      \mathrm{TDB}-\mathrm{TT} = 0.001658\sin g + 0.000014\sin 2g,\qquad
      g = 357.53^\circ + 0.9856003^\circ\, d_\mathrm{J2000},

   good to :math:`\lesssim2` ms.

.. note::

   The **position-aware** overloads ``TtToTdb(t, x_bcrs)`` /
   ``TDBToTt(t, x_bcrs)`` evaluate the full GCRS↔BCRS 4-D transform
   (Turyshev Eq. 21/22): a trapezoidal integral of the :math:`c^{-2}` and
   :math:`c^{-4}` integrands built from the Earth's barycentric velocity and the
   external Solar-System potential
   :math:`w_\mathrm{ext}=\sum_{B\neq E} GM_B/r_{EB}`, plus a position term
   :math:`-v_E\!\cdot\!r_E/c^2`.  ``TtToTdb(t, x_bcrs)`` inverts by 10-step
   Newton iteration.

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
:math:`L_L = 3.139054\times10^{-11}`:

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
   Turyshev selenoid value :math:`3.139054\times10^{-11}` (potential
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

.. _model-boundaries:

Model Boundaries
-------------------------------------------------------------------

* All API times are seconds from the per-scale J2000 origin; calendar I/O
  goes through ``GregorianToTime`` / ``TimeToGregorianString``.
* Only GPS is implemented among GNSS scales; GLONASS/Galileo/BeiDou offsets
  are documented but not exposed.
* ``ConvertTime``/``ConvertCoordinateTime`` return absolute epochs and are
  therefore bounded at :math:`\sim245` ns regardless of model quality.  Use
  ``Epoch`` or the offset accessors below that level.
* The DE440 Eq. (3) integral carries a :math:`+0.14` ns/yr secular difference
  against the DE440t ephemeris.  ``w_{0E}`` sums the ten ephemeris bodies plus
  ring models of the main asteroid belt and the Kuiper belt
  (``GM_ASTEROID_BELT``, ``GM_KUIPER_BELT``); DE440 integrates those populations
  as 343 + 30 discrete bodies, and the published masses account for
  :math:`\sim80\%` of the difference.  Use the Chebyshev fit when absolute
  agreement with JPL matters.
* ``TdbMinusTcl`` and ``TdbToLtMinusTt`` integrate from :math:`T_0` per call.
  Both auto-fit a decade-wide Chebyshev window on demand
  (``SetTdbTclAutoFit``, ``InitLtMinusTtFit``); with auto-fitting disabled a
  single lunar conversion costs seconds.
* :math:`L_L` (and hence LT) is provisional pending international
  standardization.

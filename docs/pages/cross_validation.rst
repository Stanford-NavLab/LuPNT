Cross-validation against Orekit and GMAT
=========================================

LuPNT's time-scale conversions, sidereal angles, reference-frame
conversions (Earth and Moon), Sun/Moon ephemerides, and orbit propagation
are cross-validated against two independent, widely used astrodynamics
tools:

* `Orekit <https://www.orekit.org/>`_ v13.1 (JVM-based library)
* `GMAT <https://sourceforge.net/projects/gmat/>`_ R2022a (NASA's General
  Mission Analysis Tool)

This page summarizes **how** the comparison works, the **measured
differences** per category, and the **known modeling differences** the
comparison uncovered.

Methodology
-----------

Neither tool is a build, run, or test dependency of LuPNT. Instead:

1. A developer-only generator script drives the external tool once and
   stores its outputs in a checked-in JSON fixture:

   * ``cpp/test/orekit/gen_orekit_reference.py`` (run via the optional
     ``dev`` pixi environment) writes
     ``cpp/test/orekit/data/orekit_reference.json``.
   * ``cpp/test/gmat/gen_gmat_reference.py`` (requires a local GMAT desktop
     install) generates and runs a GMAT script and writes
     ``cpp/test/gmat/data/gmat_reference.json`` (the exact GMAT script is
     checked in alongside as ``gmat_reference.script``).

2. Regular Catch2 tests under ``cpp/test/orekit/`` and ``cpp/test/gmat/``
   load the fixtures and compare LuPNT's own computations against them.
   They run as part of ``pixi run test-cpp`` on every machine, with no
   external-tool dependency.

Both fixtures use LuPNT's physical constants
(``GM_EARTH = 398600.435507e9``, ``GM_MOON = 4902.800118e9``,
``R_EARTH = 6378137``, ``J2_EARTH = 1.08262668e-3``; see
``cpp/lupnt/core/constants.h``) on both sides of the comparison, so the
observed differences are *algorithmic/data* differences, not differences in
adopted constants. Where the comparisons overlap (frame-conversion input
states, orbit cases, epochs), the two fixtures use identical inputs, so the
Orekit and GMAT columns below are directly comparable.

The test tolerances are *calibrated*: each tolerance is set a small factor
(typically 2-20x) above the measured difference, with a comment in the test
file explaining the difference's origin. A tolerance that merely "passes"
is not the goal -- the point is that each residual difference is
*understood*.

What is compared
----------------

.. list-table::
   :header-rows: 1
   :widths: 28 36 36

   * - Category
     - vs Orekit
     - vs GMAT
   * - Time scales
     - TAI-UTC, TT-TAI, TDB-TAI, GPS-TAI, UT1-UTC at 6 epochs bracketing
       the 2012/2016 leap seconds
     - TAI-UTC, TT-TAI, TDB-TAI at the same epochs (GMAT exposes no
       UT1/GPS)
   * - Sidereal angles
     - GMST (IAU-82 / IERS 1996) and Earth Rotation Angle (IAU 2000)
     - --
   * - Earth frames
     - GCRF <-> EME2000, GCRF <-> ITRF for LEO/MEO/GEO/trans-lunar states
       at 3 epochs
     - same states/epochs via GMAT's ICRF / MJ2000Eq / BodyFixed systems
   * - Moon frames
     - GCRF <-> MOON_CI (DE440 translation); MOON_ME vs the analytical IAU
       pole model
     - GCRF <-> MOON_CI (DE421); MOON_PA vs GMAT's SPICE-kernel Luna
       BodyFixed (a principal-axes frame)
   * - Ephemerides
     - Earth->Moon and Earth->Sun position/velocity (DE440) at 3 epochs
     - (implicit in the Moon-frame translations)
   * - Kepler propagation
     - 6 orbits -- Earth MEO-HEO/LEO/GEO/Molniya + lunar LLO/ELFO --
       propagated 0.25/0.5/0.75/1.5 periods (analytical
       ``KeplerianPropagator``)
     - same 6 orbits (numerical RK89 point-mass propagation)
   * - J2 acceleration
     - formula-level check at 5 positions
       (``J2OnlyPerturbation.computeAccelerationInJ2Frame``)
     - --
   * - J2 propagation
     - 3 orbits, J2 about the **inertial** Z axis (same model as LuPNT) --
       isolates integrator differences
     - 2 orbits, J2 in the **body-fixed** frame -- exposes the modeling
       difference

Measured differences
--------------------

The numbers below were measured when the fixtures were generated
(June 2026); regenerate the fixtures and re-derive them if the comparison
data changes. "Tolerance" is the value asserted by the tests.

Time scales and sidereal angles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 22 22 22 34

   * - Quantity
     - vs Orekit (max)
     - vs GMAT (max)
     - Tolerance / origin of difference
   * - TAI-UTC, TT-TAI, GPS-TAI
     - ~2e-8 s
     - ~6e-7 s
     - 1e-6 s. Exactly defined offsets; the residuals are floating-point
       representation artifacts (LuPNT stores epochs as seconds-since-J2000
       doubles; GMAT reports ~30000-day Modified Julian dates whose
       subtraction amplifies the ULP).
   * - TDB-TAI
     - ~4e-8 s
     - ~5e-7 s
     - 1e-6 s (Orekit) / 1e-5 s (GMAT). The ~1.7 ms periodic relativistic
       term is implemented with different truncated series.
   * - UT1-UTC
     - <= 0.024 s (away from leap seconds)
     - --
     - 0.1 s. EOP-table vintage difference (different IERS bulletin
       snapshots). See the leap-second note below.
   * - GMST, ERA
     - <= 1.8e-6 rad
     - --
     - 1e-5 rad. Identical formulas on both sides (IAU-82 GMST, IAU-2000
       ERA); the residual is the UT1 difference times the Earth rotation
       rate.

Frames
^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 22 22 22 34

   * - Conversion
     - vs Orekit (max)
     - vs GMAT (max)
     - Tolerance / origin of difference
   * - GCRF <-> EME2000
     - 1.8e-4 m (trans-lunar), i.e. ~5e-13 * \|r\|
     - 54 m (trans-lunar), i.e. ~1.6e-7 * \|r\|
     - 1e-3 m (Orekit); max(5 m, 3e-7 * \|r\|) (GMAT). LuPNT and Orekit
       implement the same IERS frame-bias matrix; GMAT's ICRF <->
       MJ2000Eq rotation differs by ~1e-7 rad (slightly epoch-dependent).
   * - GCRF <-> ITRF
     - 4.4 km @ trans-lunar 2025 epoch (~1.3e-5 * \|r\|); 14-98 m @ LEO
     - 1.6 km @ trans-lunar (~4.5e-6 * \|r\|); 34-36 m @ LEO
     - max(100 m, 2e-5 * \|r\|), velocity 1 m/s. Earth-orientation
       (UT1/polar-motion) table vintage: a ~tens-of-ms UT1 difference
       rotates positions by ~1e-6..1e-5 rad about the polar axis, scaling
       with radius. Differences grow toward 2025 epochs where one or both
       tables are predictions.
   * - GCRF -> MOON_CI
     - 0.33 m
     - 80 m
     - 1 m (Orekit) / 200 m (GMAT). Pure DE-ephemeris translation: LuPNT
       and Orekit both evaluate DE440 (residual = distribution/reader
       differences, ~5e-11 relative); GMAT shipped DE421, so its column
       measures the DE421-vs-DE440 lunar ephemeris offset.
   * - Moon body-fixed
     - MOON_ME vs IAU model: 57-200 m (~3e-5 * \|r\|)
     - MOON_PA vs Luna BodyFixed: 56-86 m (translation-dominated)
     - 3e-4 * \|r\| (Orekit) / 300 m (GMAT). Orekit's body-oriented Moon
       frame is the analytical IAU pole model, which approximates the DE
       *mean-Earth* axes to ~3e-5 rad. GMAT loads a SPICE Luna
       principal-axes kernel, validating LuPNT's MOON_PA: after removing
       the DE421-vs-DE440 translation, the PA orientations agree at the
       sub-arcsecond level. Cross-check: swapping the frames (PA vs IAU,
       or ME vs GMAT) increases the differences 10-40x (0.6-3.4 km),
       confirming each identification.
   * - Ephemerides (Earth->Moon, Earth->Sun)
     - Moon 0.33 m / Sun 7.2 m (both ~5e-11 relative)
     - --
     - max(1 m, 1e-10 * \|r\|), velocity 1e-5 m/s. Same DE440 data,
       independently distributed and independently interpolated.

Orbit propagation
^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 22 22 22 34

   * - Quantity
     - vs Orekit (max)
     - vs GMAT (max)
     - Tolerance / origin of difference
   * - Kepler: initial COE -> Cartesian
     - 1.5e-8 m
     - 1.3e-8 m
     - 1e-3 m. Pure conversion math; machine precision.
   * - Kepler: propagated Cartesian
     - 5.7e-8 m / 5e-12 m/s
     - 1.8e-2 m / 8e-6 m/s
     - 1e-2 m (Orekit) / 0.1 m (GMAT). Orekit's propagator is analytical
       (like LuPNT's), so agreement is machine precision; GMAT numerically
       integrates the point-mass problem (RK89, accuracy 1e-13), so its
       column shows integrator drift.
   * - Kepler: propagated angles (M/E/nu)
     - <= 3e-15 rad
     - <= 5e-8 rad
     - 1e-9 rad (Orekit) / 1e-6 rad (GMAT). Same analytical-vs-numerical
       distinction.
   * - J2 acceleration
     - <= 1e-15 m/s^2 (vs 1e-12 tolerance)
     - --
     - Identical closed-form J2 formula.
   * - J2 propagation (inertial-axis J2)
     - 2.2e-5 m / 1.9e-9 m/s over 1.5 orbits
     - --
     - 1e-3 m / 1e-6 m/s. Same force model on both sides (J2 about the
       inertial Z axis); the residual is LuPNT's fixed-step RK4 (1 s) vs
       Orekit's adaptive DP853.
   * - J2 propagation (body-fixed J2)
     - --
     - 22-182 m / 0.005-0.099 m/s over 0.25-1.5 orbits
     - 300 m / 0.15 m/s. **Modeling difference**: LuPNT's
       ``JToCartTwoBodyDynamics`` takes the inertial Z axis as the
       oblateness symmetry axis, GMAT evaluates the zonal harmonic in the
       Earth body-fixed frame. The Orekit row above (same dynamics, same
       axis convention, ~1e-5 m agreement) isolates the axis choice as the
       cause.

Known modeling differences and limitations
------------------------------------------

The cross-validation surfaced three noteworthy items, all documented in the
corresponding test files:

**UT1 interpolation across leap seconds (LuPNT limitation).**
LuPNT's EOP table stores daily UT1-UTC values and linearly interpolates
*across* the 1 s leap-second jump. For epochs within ~1 day of an
insertion this produces spurious UT1-UTC: measured +0.938 s one hour
before the 2016-12-31 leap and +0.094 s six hours after the 2012-06-30
leap. Orekit interpolates the continuous UT1-TAI and is unaffected. The
fixtures flag such epochs (``near_leap_second``) and the tests skip the
UT1-dependent checks (UT1-UTC, GMST, ERA) there; TAI/TT/TDB/GPS are still
verified at those epochs.

**J2 symmetry axis (intentional simplification).**
``JToCartTwoBodyDynamics`` applies the J2 acceleration about the inertial
Z axis rather than the body-fixed pole. The GMAT comparison quantifies the
effect of this choice for Earth orbits: tens to ~200 m over up to 1.5
orbits for the tested LEO/MEO-HEO cases. Use the full gravity-field
dynamics when this matters.

**GMAT's ICRF/MJ2000Eq frame bias.**
GMAT's ICRF <-> MJ2000Eq rotation differs from the IERS frame-bias matrix
used by LuPNT and Orekit (which agree to sub-mm) by ~1e-7 rad. This is a
GMAT-side convention difference, visible only at trans-lunar ranges
(~50 m).

Regenerating the fixtures
-------------------------

See `cpp/test/orekit/README.md
<https://github.com/Stanford-NavLab/LuPNT/blob/main/cpp/test/orekit/README.md>`_
and `cpp/test/gmat/README.md
<https://github.com/Stanford-NavLab/LuPNT/blob/main/cpp/test/gmat/README.md>`_
for the per-tool regeneration workflows (optional ``dev`` pixi environment
for Orekit; a local GMAT R2022a desktop install for GMAT) and the full JSON
data dictionaries.

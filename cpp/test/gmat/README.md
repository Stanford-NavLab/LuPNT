# GMAT cross-validation

This directory holds the GMAT cross-validation tests and their
**pre-generated, checked-in fixture** (`data/gmat_reference.json`), used to
cross-check LuPNT's time-scale conversions, frame conversions (Earth and
Moon), two-body Keplerian propagation (Earth and lunar orbits), and
J2-perturbed propagation against
[GMAT](https://sourceforge.net/projects/gmat/) (R2022a).

**Normal users and CI never need GMAT.** The tests in this directory:

- `test_time_conversions_gmat.cc` -- TAI/TT/TDB time scales
- `test_frame_conversions_gmat.cc` -- GCRF/EME2000/ITRF and Moon frames
- `test_orbit_dynamics_gmat.cc` -- Kepler (Earth + lunar) and J2 propagation

only read `data/gmat_reference.json` (via `LoadTestJson()` in
`cpp/test/utils.cc`) and run as part of the regular `pixi run test-cpp`
suite, with no GMAT dependency at build or run time.

Unlike Orekit (an optional `dev` pixi/JVM dependency), **GMAT is a
separately installed desktop application**, not something pixi can install.
It is only needed by developers who want to *regenerate* the fixture.

A summary of the observed LuPNT-vs-GMAT differences per category is in the
documentation page [`docs/pages/cross_validation.rst`](../../../docs/pages/cross_validation.rst).

## What's in `data/gmat_reference.json`

| Section          | Contents | Consumed by |
|-------------------|----------|-------------|
| `_meta`           | Provenance: generator script, GMAT version, planetary-ephemeris source (DE421), the physical constants used (`GM_EARTH`, `GM_MOON`, `R_EARTH`, `J2_EARTH`), and notes documenting the GMAT-specific caveats below. | (informational) |
| `time_scales`     | For six epochs (the same as the Orekit fixture, bracketing the 2012/2016 leap seconds), the offsets `TAI-UTC`, `TT-TAI`, `TDB-TAI` as read from a GMAT `Spacecraft`'s `*ModJulian` epoch parameters. | `test_time_conversions_gmat.cc` |
| `frames`          | For three epochs and four position/velocity states (LEO, MEO, GEO, trans-lunar; same inputs as the Orekit fixture), the same state in GMAT's `EarthICRFCS` (ICRF axes = LuPNT's GCRF), `EarthMJ2000Eq` (= EME2000), and `EarthFixedCS` (BodyFixed = ITRF). | `test_frame_conversions_gmat.cc` |
| `moon_frames`     | For two epochs and two Moon-centered states (same lunar-centric offsets as the Orekit fixture), the state in Earth-ICRF, Luna-centered ICRF (`r_moon_ci`), and GMAT's Luna `BodyFixed` system (`r_moon_fixed`, a DE-era lunar *principal axes* frame loaded from a SPICE kernel). | `test_frame_conversions_gmat.cc` |
| `kepler`          | Six cases -- four Earth orbits (MEO-HEO, LEO, GEO, Molniya; point-mass-only `TwoBodyFM`) and two lunar orbits (LLO, ELFO; `LunaTwoBodyFM` with the overridden `Luna.Mu`) -- with initial Cartesian state, anomalies, period, and the propagated state at 0.25/0.5/0.75/1.5 orbital periods. | `test_orbit_dynamics_gmat.cc` (`dynamics.kepler_gmat_reference`) |
| `j2_propagation`  | The first two Earth cases propagated with `J2FM` (point mass + J2 only, via the custom `lupnt_j2_earth.cof` potential file). | `test_orbit_dynamics_gmat.cc` (`dynamics.j2_propagation_gmat_reference`) |

Within `kepler[i]`/`j2_propagation[i]`:

- `gm` (and `j2`/`r_earth` for `j2_propagation`) are the exact constants
  from `cpp/lupnt/core/constants.h`; for the Earth they are also baked into
  `lupnt_j2_earth.cof`.
- `coe0` is `[a, e, i, Omega, w, M0]` (SI units / radians, **mean**
  anomaly); `body` is `"EARTH"` or `"MOON"`.
- `cart0`/`propagated[].cart` are GMAT's Cartesian states (`[x..vz]`,
  m / m/s), in `EarthMJ2000Eq` axes for the Earth cases and Luna-centered
  ICRF axes for the lunar cases.

`data/gmat_reference.script` is a checked-in copy of the exact GMAT script
the generator ran (with machine-specific paths replaced by `<...>` tokens),
for documentation and debugging.

## GMAT-specific caveats

These are also recorded in `_meta.notes` of `data/gmat_reference.json`:

1. **No UT1 or GPS time scale.** GMAT R2022a does not expose a
   `UT1ModJulian` or `GPSModJulian` spacecraft parameter, so `time_scales`
   only has `tai_minus_utc`, `tt_minus_tai`, and `tdb_minus_tai`.

2. **`coe0`/anomalies/`period` are computed in Python, not by GMAT.**
   GMAT's Keplerian display state only accepts **True Anomaly** as the 6th
   element, so `gen_gmat_reference.py` solves Kepler's equation itself and
   uses the resulting True Anomaly only to configure each spacecraft's
   initial state. `cart0` and every `propagated` entry are genuine GMAT
   outputs.

3. **DE421, not DE440.** GMAT R2022a ships DE405/DE421/DE424; the fixture
   was generated with `SolarSystem.EphemerisSource = 'DE421'`, the closest
   to LuPNT's DE440. The Moon-centered comparisons therefore absorb a
   ~55-85 m DE421-vs-DE440 lunar-ephemeris offset (measured at the
   2024/2025 fixture epochs) on top of any algorithmic difference.

4. **GMAT's Luna `BodyFixed` is a principal-axes frame.** GMAT loads a
   SPICE Luna frame kernel for the lunar orientation, which is a DE-era
   *principal axes* (PA) frame. `test_frame_conversions_gmat.cc` therefore
   compares it against LuPNT's `MOON_PA` (DE440 PA): the observed ~55-90 m
   difference is dominated by the DE421-vs-DE440 translation, i.e. the PA
   orientations agree at the sub-arcsecond level. (Comparing against
   `MOON_ME` instead fails by ~10-40x more, ~0.6-3.4 km.) This is
   complementary to the Orekit comparison, which validates `MOON_ME`
   against the analytical IAU model.

5. **Custom potential file (`lupnt_j2_earth.cof`).** GMAT's default Earth
   gravity model uses NASA/EGM constants that differ slightly from LuPNT's.
   `lupnt_j2_earth.cof` is a minimal potential file setting `GM`, `R_eq`,
   and `C20` (= `-J2` normalized) to LuPNT's exact values:
   - `TwoBodyFM` uses `Degree=0, Order=0` (pure point mass).
   - `J2FM` uses `Degree=2, Order=0` (point mass + J2 zonal term only).

   `Earth.Mu`/`Earth.EquatorialRadius` (and `Luna.Mu`) are also overridden
   globally so GMAT's internal Cartesian <-> Keplerian conversions use the
   same constants.

6. **The "combined Propagate" bug.** Combining two spacecraft/propagators
   in one synced command, e.g.
   `Propagate Prop1(SatA) Prop2(SatB) {SatA.ElapsedSecs = X};`, was found
   to desync the second spacecraft's force model in GMAT R2022a (its state
   ended up ~60x further from the two-body solution than the expected J2
   perturbation -- clearly a GMAT bug, not a real J2 effect). The generator
   always emits separate, sequential
   `Propagate <Prop>(<Sat>) {<Sat>.ElapsedSecs = <dt>};` commands, one per
   spacecraft/propagator pair. If you extend the generator, preserve this
   pattern.

7. **GMAT's "two-body" propagation is numerical.** Unlike Orekit's
   analytical `KeplerianPropagator`, GMAT numerically integrates even the
   point-mass problem (RungeKutta89, `Accuracy = 1e-13`), so the `kepler`
   comparisons carry integrator drift (~2e-2 m, ~5e-8 rad in the fast
   angles over 1.5 orbits) and the test tolerances are correspondingly
   looser than the Orekit ones.

8. **J2 modeling axis difference (`j2_propagation`).** LuPNT's
   `JToCartTwoBodyDynamics` evaluates the J2 acceleration about the
   *inertial* (GCRF) Z axis, whereas GMAT's `GravityField` evaluates the
   degree-2 zonal term in the Earth *body-fixed* frame. Over 0.25-1.5
   orbits this expected modeling difference accumulates to ~20-180 m /
   ~0.005-0.1 m/s (the same comparison against Orekit with an
   inertial-frame J2 agrees to ~1e-5 m, isolating the axis choice as the
   cause).

## Regenerating the fixture (developers only)

`data/gmat_reference.json` is produced by `gen_gmat_reference.py`, which
*generates* the GMAT script from its case tables, runs it, and parses the
report output. You only need to run this if you're adding new
cross-validation cases, or deliberately updating the reference values.

1. Install [GMAT](https://sourceforge.net/projects/gmat/) R2022a (or a
   compatible version) -- a normal desktop application install, not a
   pixi/conda package. On macOS the default install location is
   `/Applications/GMAT_R2022a/`.

2. Point the generator at `GmatConsole`, either via
   `GMAT_CONSOLE=<full path to the executable>` or
   `GMAT_HOME`/`GMAT_ROOT_DIR=<GMAT install directory>`. Common default
   install locations are searched automatically.

3. Run the generator from the repo root (only the Python standard library
   is used):
   ```bash
   GMAT_CONSOLE=/Applications/GMAT_R2022a/bin/GmatConsole-R2022a \
       python3 cpp/test/gmat/gen_gmat_reference.py
   ```
   A successful run prints
   `*** GMAT Integration test (Console version) successful! ***` and
   overwrites `data/gmat_reference.json` and `data/gmat_reference.script`.

4. Review the diff, then run `pixi run test-cpp` to confirm the
   `*_gmat_reference` tests still pass (or update the affected
   tests/tolerances if you intentionally changed something).

## Tolerance rationale

The comparisons in the `*_gmat_reference` tests are not exact-equality
checks; each test file documents its tolerances next to the assertions, and
`docs/pages/cross_validation.rst` tabulates the observed differences. In
brief:

- **Time scales**: ~2-6e-7 s artifacts from GMAT's ModJulian-day
  representation (subtracting ~30000-day MJDs amplifies the double ULP);
  1 us tolerance for TAI/TT, 10 us for TDB (different truncated series for
  the periodic relativistic term).
- **GCRF <-> EME2000**: GMAT's ICRF->MJ2000Eq rotation differs from
  LuPNT/Orekit's frame bias by ~1e-7 rad (slightly epoch-dependent), so the
  position tolerance scales as `3e-7 * |r|` (5 m floor).
- **GCRF <-> ITRF**: EOP-vintage differences, same `max(100 m, 2e-5 * |r|)`
  scaling as the Orekit comparison.
- **Moon frames**: dominated by the DE421-vs-DE440 lunar ephemeris offset
  (~55-85 m); 200 m (MOON_CI) and 300 m (MOON_PA vs Luna BodyFixed)
  tolerances.
- **Kepler propagation**: GMAT RK89 integrator drift; 0.1 m / 1e-4 m/s for
  the propagated Cartesian state, 1e-6 rad for the fast angles.
- **J2 propagation**: the body-fixed-vs-inertial J2 axis modeling
  difference; flat 300 m / 0.15 m/s.

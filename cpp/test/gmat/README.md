# GMAT cross-validation

This directory holds the GMAT cross-validation tests and their
**pre-generated, checked-in fixture** (`data/gmat_reference.json`), used to
cross-check LuPNT's time-scale conversions, frame conversions (Earth and
Moon), two-body Keplerian propagation (Earth and lunar orbits), and
J2-perturbed propagation against
[GMAT](https://sourceforge.net/projects/gmat/) (R2026a).

A companion fixture (`data/gmat_forces_reference.json`, produced by
`gen_gmat_forces_reference.py`) extends the cross-validation to three higher
force models -- spherical-harmonic gravity, third-body point mass, and solar
radiation pressure -- and is consumed by `test_force_models_gmat.cc`. See
[What's in `data/gmat_forces_reference.json`](#whats-in-datagmat_forces_referencejson)
below.

**Normal users and CI never need GMAT.** The tests in this directory:

- `test_time_conversions_gmat.cc` -- TAI/TT/TDB time scales
- `test_frame_conversions_gmat.cc` -- GCRF/EME2000/ITRF and Moon frames
- `test_orbit_dynamics_gmat.cc` -- Kepler (Earth + lunar) and J2 propagation
- `test_force_models_gmat.cc` -- gravity 8x8, third-body, and SRP propagation
  (reads `data/gmat_forces_reference.json`)

only read the checked-in JSON fixtures (via `LoadTestJson()` in
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
| `_meta`           | Provenance: generator script, GMAT version (R2026a), planetary-ephemeris source (DE440, supplied as LuPNT's own `de440.bsp` via SPICE), the physical constants used (`GM_EARTH`, `GM_MOON`, `R_EARTH`, `J2_EARTH`), and notes documenting the GMAT-specific caveats below. | (informational) |
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

1. **No UT1 or GPS time scale.** GMAT R2026a does not expose a
   `UT1ModJulian` or `GPSModJulian` spacecraft parameter, so `time_scales`
   only has `tai_minus_utc`, `tt_minus_tai`, and `tdb_minus_tai`.

2. **`coe0`/anomalies/`period` are computed in Python, not by GMAT.**
   GMAT's Keplerian display state only accepts **True Anomaly** as the 6th
   element, so `gen_gmat_reference.py` solves Kepler's equation itself and
   uses the resulting True Anomaly only to configure each spacecraft's
   initial state. `cart0` and every `propagated` entry are genuine GMAT
   outputs.

3. **DE440 via SPICE (both sides).** GMAT ships DE405/DE421/DE424 and still
   has no DE440 even in R2026a, so it is pointed at LuPNT's own `de440.bsp`
   with `SolarSystem.EphemerisSource = 'SPICE'` and
   `SolarSystem.SPKFilename = <de440.bsp>`. Both sides therefore evaluate the
   same DE440 planetary ephemeris; the former DE421-vs-DE440 Moon-frame offset
   is gone, and the residual in the Moon-centered comparisons is now at the
   ephemeris distribution/interpolation level (GMAT's SPICE reader vs LuPNT's),
   not a DE-version offset.

4. **GMAT's Luna `BodyFixed` is a principal-axes frame.** GMAT loads a
   SPICE Luna frame kernel for the lunar orientation, which is a DE-era
   *principal axes* (PA) frame. `test_frame_conversions_gmat.cc` therefore
   compares it against LuPNT's `MOON_PA` (DE440 PA): with both sides on DE440
   the Moon-position translation agrees at the ephemeris distribution level and
   the PA orientations agree at the sub-arcsecond level, so the residual is well
   within tolerance. (Comparing against `MOON_ME` instead fails by ~10-40x
   more, ~0.6-3.4 km.) This is complementary to the Orekit comparison, which
   validates `MOON_ME` against the analytical IAU model.

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
   to desync the second spacecraft's force model in GMAT (its state
   ended up ~60x further from the two-body solution than the expected J2
   perturbation -- clearly a GMAT bug, not a real J2 effect; the R2026a
   generators keep the same workaround). The generator
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

## What's in `data/gmat_forces_reference.json`

This companion fixture cross-validates the three higher force models that the
`gmat_reference.json` scenarios above do not exercise. It is **propagation
level**: GMAT reports spacecraft states, not accelerations, so each case
propagates a spacecraft with a *single* force model active and stores the
Cartesian state after a fixed 2-hour arc. `test_force_models_gmat.cc` re-runs
the same arc through LuPNT's `NBodyDynamics` (RK8, 10 s step) and compares. This
isolates the assembled propagation (force model + integrator + frame
realization), and complements the Orekit checks in
`cpp/test/orekit/test_force_models_orekit.cc`, which validate the same LuPNT
acceleration *algorithms* at the formula level.

All cases share one epoch, **2024-03-15T12:00:00 UTC**, and are reported at
elapsed times `dt = 1800, 3600, 7200 s`. Initial and propagated states are
GMAT's `EarthMJ2000Eq` Cartesian states (`[x..vz]`, m / m/s); LuPNT consumes
them as GCRF (the frame-bias difference is sub-metre at these radii and is
absorbed by the tolerances).

| Section                  | Contents | Consumed by |
|--------------------------|----------|-------------|
| `_meta`                  | Provenance: generator, `gmat_version` (`R2026a`), `ephemeris` (`DE440 (LuPNT's, via SPICE)`), and the constants forced onto GMAT (`GM_EARTH`, `R_EARTH`, `GM_SUN`, `GM_MOON`, `solar_flux`). | (informational) |
| `gravity_propagation`    | 8x8 spherical-harmonic Earth gravity. `cof_file` (`EGM96.cof`, the same file LuPNT loads), `n_max`, `m_max`, and `cases[]`. | `dynamics.gravity_propagation_gmat_reference` |
| `third_body_propagation` | Sun + Moon third-body point-mass perturbation on an Earth orbit. `perturbers` (`["SUN","MOON"]`), `GM_sun`, `GM_moon`, and `cases[]`. | `dynamics.third_body_propagation_gmat_reference` |
| `srp_propagation`        | Cannonball solar radiation pressure on a high, continuously sunlit orbit. `Cr`, `area`, `mass`, `solar_flux`, and `cases[]`. | `dynamics.srp_propagation_gmat_reference` |

Every `cases[i]` has:

- `cart0` -- the initial Cartesian state (`[x..vz]`, m / m/s).
- `propagated[]` -- one entry per elapsed time, each `{dt, cart}` with `dt` in
  seconds and `cart` the GMAT state at that time.

The force models are configured to match LuPNT exactly: `GravFM` uses the shared
`EGM96.cof` at degree/order 8 with LuPNT's `GM`/`R_eq`; `TbFM` uses
`PointMasses = {Earth, Sun, Luna}` on DE440; `SrpFM` uses a `Spherical`
(cannonball) SRP with `Flux = 1361 W/m^2`, `Cr = 1.5`, `area = 4 m^2`,
`mass = 500 kg`. The SRP orbit is deliberately high and sunward so the
satellite is never in eclipse -- LuPNT's apparent-disk penumbra and GMAT's
shadow model differ, and a sunlit arc keeps that difference out of the
cannonball comparison. Because LuPNT's SRP takes the Sun position from SPICE
internally, the Sun is *not* added as a gravitating body in the SRP case (which
would add third-body gravity GMAT's `SrpFM` does not have).

Measured worst-case LuPNT-vs-GMAT R2026a differences over the 2-hour arc:

| Force model | Observed (position) | Tolerance |
|-------------|---------------------|-----------|
| Gravity 8x8 | 0.19 m              | 3 m / 5e-3 m/s |
| Third body (Sun + Moon) | 0.004 m | 0.2 m / 1e-3 m/s |
| SRP (cannonball, sunlit) | 0.003 m | 0.2 m / 1e-3 m/s |

The residuals are integrator and frame-realization (ITRF vs GMAT `EarthFixed`)
differences only: LuPNT's constants, the shared `EGM96.cof`, and the shared
DE440 remove every adopted-constant and ephemeris-version difference.

`data/gmat_forces_reference.script` is the checked-in copy of the exact GMAT
script the generator ran (machine paths tokenised). Like `gen_gmat_reference.py`
it emits one `Propagate <Prop>(<Sat>) {...}` per spacecraft-propagator pair
(never a combined `Propagate` -- see caveat 6 above).

## Regenerating the fixtures (developers only)

`data/gmat_reference.json` is produced by `gen_gmat_reference.py`, and
`data/gmat_forces_reference.json` by `gen_gmat_forces_reference.py`. Each
*generates* a GMAT script from its case tables, runs it, and parses the report
output. You only need to run these if you're adding new cross-validation cases,
or deliberately updating the reference values.

1. Install [GMAT](https://sourceforge.net/projects/gmat/) R2026a (or a
   compatible version) -- a normal desktop application install, not a
   pixi/conda package. On macOS the default install location is
   `/Applications/GMAT_R2026a/`. GMAT ships DE405/DE421/DE424 only, so both
   generators point it at LuPNT's own `de440.bsp` via
   `SolarSystem.EphemerisSource = 'SPICE'` /
   `SolarSystem.SPKFilename = <de440.bsp>`, giving both sides the same DE440.

2. Point the generator at `GmatConsole`, either via
   `GMAT_CONSOLE=<full path to the executable>` or
   `GMAT_HOME`/`GMAT_ROOT_DIR=<GMAT install directory>`. Common default
   install locations are searched automatically.

3. Run the generators from the repo root (only the Python standard library
   is used):
   ```bash
   GMAT_CONSOLE=/Applications/GMAT_R2026a/bin/GmatConsole-R2026a \
       python3 cpp/test/gmat/gen_gmat_reference.py
   GMAT_CONSOLE=/Applications/GMAT_R2026a/bin/GmatConsole-R2026a \
       python3 cpp/test/gmat/gen_gmat_forces_reference.py
   ```
   The first (`gen_gmat_reference.py`) overwrites `data/gmat_reference.json`
   and `data/gmat_reference.script`; the second
   (`gen_gmat_forces_reference.py`) overwrites
   `data/gmat_forces_reference.json` and `data/gmat_forces_reference.script`.
   The forces generator locates `de440.bsp` under `LuPNT_data/` by default;
   override it with `LUPNT_DE440=<path>`.

4. Review the diff, then run `pixi run test-cpp` to confirm the
   `*_gmat_reference` tests (including the new force-model cases) still pass
   (or update the affected tests/tolerances if you intentionally changed
   something).

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
- **Moon frames**: with both sides on DE440 (GMAT via `de440.bsp`/SPICE) the
  residual is now the ephemeris distribution/interpolation level plus the
  sub-arcsecond PA orientation agreement; the 200 m (MOON_CI) and 300 m
  (MOON_PA vs Luna BodyFixed) tolerances are unchanged and pass with extra
  margin.
- **Kepler propagation**: GMAT RK89 integrator drift; 0.1 m / 1e-4 m/s for
  the propagated Cartesian state, 1e-6 rad for the fast angles.
- **J2 propagation**: the body-fixed-vs-inertial J2 axis modeling
  difference; flat 300 m / 0.15 m/s.
- **Gravity 8x8 propagation** (`gmat_forces_reference.json`): shared
  `EGM96.cof` and DE440 leave only integrator + Earth-orientation (ITRF vs
  GMAT `EarthFixed`) differences; observed 0.19 m over 2 h, 3 m / 5e-3 m/s
  tolerance.
- **Third-body propagation**: Sun + Moon point masses on shared DE440;
  observed 0.004 m over 2 h, 0.2 m / 1e-3 m/s tolerance.
- **SRP propagation**: cannonball on a continuously sunlit orbit (shadow
  models never engage), matched flux/Cr/area/mass; observed 0.003 m over 2 h,
  0.2 m / 1e-3 m/s tolerance.

# Orekit cross-validation

This directory holds the Orekit cross-validation tests and their
**pre-generated, checked-in fixture** (`data/orekit_reference.json`), used to
cross-check LuPNT's time-scale conversions, sidereal/Earth-rotation angles,
frame conversions (Earth and Moon), Sun/Moon ephemerides, two-body Keplerian
propagation, J2 acceleration/propagation, and the spherical-harmonic gravity,
third-body point-mass and solar-radiation-pressure force models against
[Orekit](https://www.orekit.org/) (v13.1).

**Normal users and CI never need Orekit/Java.** The tests in this directory:

- `test_time_conversions_orekit.cc` -- time scales, GMST, ERA
- `test_frame_conversions_orekit.cc` -- GCRF/EME2000/ITRF and Moon frames
- `test_ephemeris_orekit.cc` -- Sun/Moon DE440 ephemerides
- `test_orbit_dynamics_orekit.cc` -- Kepler propagation, J2 acceleration and
  propagation
- `test_force_models_orekit.cc` -- spherical-harmonic gravity, third-body
  point-mass and SRP accelerations (formula level)

only read `data/orekit_reference.json` (via `LoadTestJson()` in
`cpp/test/utils.cc`) and run as part of the regular `pixi run test-cpp`
suite, with no JVM/Orekit dependency at build or run time.

A summary of the observed LuPNT-vs-Orekit differences per category is in the
documentation page [`docs/pages/cross_validation.rst`](../../../docs/pages/cross_validation.rst).

## What's in `data/orekit_reference.json`

| Section          | Contents | Consumed by |
|-------------------|----------|-------------|
| `_meta`           | Provenance: generator script, Orekit version, `orekit-data` source, and the physical constants (`GM_EARTH`, `GM_MOON`, `R_EARTH`, `J2_EARTH`) used so the comparison isolates *algorithmic* differences from differences in adopted constants. | (informational) |
| `time_scales`     | For six epochs (bracketing the 2012-06-30 and 2016-12-31 leap seconds), the offsets `TAI-UTC`, `TT-TAI`, `TDB-TAI`, `GPS-TAI`, and `UT1-UTC` from Orekit's `TimeScalesFactory`. Epochs within ~1 day of a leap insertion carry `near_leap_second: true` (see "Known limitation" below). | `test_time_conversions_orekit.cc` |
| `sidereal`        | GMST (IERS 1996 conventions == the classic IAU-82 formula LuPNT implements) and the IAU-2000 Earth Rotation Angle at the same six epochs. | `test_time_conversions_orekit.cc` |
| `frames`          | For three epochs and four position/velocity states (LEO, MEO, GEO, trans-lunar), the same state in `GCRF`, `EME2000`, and `ITRF` (IERS 2010, no `simpleEOP`). | `test_frame_conversions_orekit.cc` |
| `moon_frames`     | For two epochs and two Moon-centered states (LLO-like and ELFO-like), the same state in GCRF, Moon-centered inertial (`r_moon_ci`, a pure translation by the DE440 Earth->Moon vector), and Orekit's body-oriented Moon frame (`r_moon_fixed_iau`, the analytical IAU pole model). | `test_frame_conversions_orekit.cc` |
| `ephemerides`     | Position/velocity of the Moon and the Sun relative to the Earth (GCRF axes, DE440) at three epochs. | `test_ephemeris_orekit.cc` |
| `kepler`          | Six `ClassicalOE` cases -- four Earth orbits (MEO-HEO, LEO, GEO, Molniya) and two lunar orbits (LLO, ELFO, with `GM_MOON`) -- with the initial Cartesian state, eccentric/true anomalies, orbital period, and the propagated state at 0.25/0.5/0.75/1.5 orbital periods via Orekit's `KeplerianOrbit`/`KeplerianPropagator`. | `test_orbit_dynamics_orekit.cc` (`dynamics.kepler_orekit_reference`) |
| `j2_acceleration` | For five sample positions, the J2 perturbing acceleration from Orekit's `J2OnlyPerturbation.computeAccelerationInJ2Frame`. | `test_orbit_dynamics_orekit.cc` (`dynamics.j2_acceleration_orekit_reference`) |
| `j2_propagation`  | Three cases (MEO-HEO, LEO, Molniya) numerically propagated with `NumericalPropagator` + `NewtonianAttraction` + `J2OnlyPerturbation` with **GCRF** as the J2 frame -- the same inertial-Z-axis modeling choice as LuPNT's `JToCartTwoBodyDynamics`, so this comparison isolates pure numerical-integration differences (unlike the GMAT one, which exposes the body-fixed-vs-inertial J2 modeling difference). | `test_orbit_dynamics_orekit.cc` (`dynamics.j2_propagation_orekit_reference`) |
| `gravity_acceleration` | Spherical-harmonic gravity acceleration, **formula level**. Keyed by body (`EARTH`, 8x8 EGM96; `MOON`, 12x12 GRGM1200B); each holds `cof_file`/`n_max`/`m_max`/`GM`/`R` and `cases[]` of `{r_bf, a_bf}` (body-fixed position and Orekit's `HolmesFeatherstoneAttractionModel` full-field acceleration, with the two-body term added back). Both sides use the **same** LuPNT `.cof`: the generator parses it and writes an Orekit-readable ICGEM copy, so the check isolates the harmonic recursion + normalization convention, not the coefficients. | `test_force_models_orekit.cc` (`dynamics.gravity_field_orekit_reference`) |
| `third_body_acceleration` | Third-body point-mass acceleration, **formula level**. `EARTH_CENTERED[]` of `{center, r, perturbers[]}`; each perturber is `{body, GM, s, a}` -- the perturber's DE440 position `s` (in the central-body inertial frame) and Orekit's `ThirdBodyAttraction` acceleration `a`. Storing `s` decouples the ephemeris (checked by `ephemerides`) from the perturbation formula. Sun and Moon on two Earth orbits (4 cases). | `test_force_models_orekit.cc` (`dynamics.third_body_orekit_reference`) |
| `srp_acceleration` | Cannonball solar-radiation-pressure acceleration, **formula level**, fully sunlit (Orekit `SolarRadiationPressure`, lighting ratio asserted `== 1`). `cases[]` of `{r, r_sun, Cr, area, mass, P0, AU, a}`; `P0` is Orekit's reference pressure at 1 AU (stored so LuPNT uses the same value), so the check isolates the cannonball formula and its Sun-direction / inverse-square / AU^2 convention, not the flux. | `test_force_models_orekit.cc` (`dynamics.srp_orekit_reference`) |

The last three sections are **formula-level** force-model checks: LuPNT's
`AccelarationGravityField`, `AccelerationPointMass` and
`AccelerationSolarRadiation` are fed the exact inputs Orekit's force models used
and the acceleration vector is compared directly, so there is no integrator or
frame noise -- only the force-model algorithm is exercised. This is the same
style as `j2_acceleration`, and it complements the GMAT force-model checks in
`cpp/test/gmat/`, which are propagation-level (GMAT reports states, not
accelerations). The observed worst-case agreements are at machine precision:

| Section | Orekit model | Observed | Tolerance |
|---------|--------------|----------|-----------|
| `gravity_acceleration` (Earth 8x8 EGM96, Moon 12x12 GRGM1200B, 42 cases) | `HolmesFeatherstoneAttractionModel` | 6.2e-15 m/s^2 | 1e-9 m/s^2 |
| `third_body_acceleration` (Sun + Moon, 4 cases) | `ThirdBodyAttraction` | 1.9e-18 m/s^2 | 1e-12 m/s^2 |
| `srp_acceleration` (cannonball, sunlit, 3 cases) | `SolarRadiationPressure` | 6.7e-24 m/s^2 | 1e-15 m/s^2 |

## Regenerating the fixture (developers only)

`data/orekit_reference.json` is produced by `gen_orekit_reference.py`. You
only need to run this if you're adding new cross-validation cases, or
deliberately updating the reference values (e.g. after an intentional
algorithm change).

1. Make sure the optional `dev` pixi environment is installed. It adds the
   JVM-based `orekit` Python package (via JCC) on top of the default
   environment:
   ```bash
   pixi install -e dev
   ```
2. Run the generator from the repo root:
   ```bash
   pixi run -e dev python cpp/test/orekit/gen_orekit_reference.py
   ```
   On first run, this downloads the small
   [`orekit-data`](https://gitlab.orekit.org/orekit/orekit-data) bundle
   (leap seconds, EOP, DE440 ephemerides, ...) into
   `data/orekit-data-main/` (under the repo's top-level `data/` directory).
   That directory is matched by the `**/*orekit-data*/` pattern in
   `.gitignore`, so it is never committed.
3. The script overwrites `cpp/test/orekit/data/orekit_reference.json`.
   Review the diff, then run `pixi run test-cpp` to confirm the
   `*_orekit_reference` tests still pass (or update the affected
   tests/tolerances if you intentionally changed something).

## Known limitation: UT1 near leap seconds

LuPNT's EOP table stores daily UT1-UTC values and linearly interpolates
*across* the 1 s leap-second jump, producing spurious UT1-UTC (up to
~0.94 s one hour before an insertion, ~0.09 s six hours after) for epochs
within ~1 day of a leap second. Orekit interpolates the continuous UT1-TAI
instead and is unaffected. The fixture flags such epochs with
`near_leap_second: true`, and the tests skip the UT1-dependent checks
(UT1-UTC, GMST, ERA) there while still checking the constant-offset scales
(TAI/TT/TDB/GPS). See `docs/pages/cross_validation.rst` for details.

## Tolerance rationale

The comparisons in the `*_orekit_reference` tests are not exact-equality
checks. The main sources of unavoidable disagreement are documented in the
test files themselves and summarized (with the observed magnitudes) in
`docs/pages/cross_validation.rst`. In brief:

- **Floating-point representation of large epochs** (time-scale tests):
  ~1e-7 s artifacts from representing epochs as seconds-since-J2000 doubles;
  1 us tolerances.
- **EOP table "vintage"** (UT1, GMST/ERA, and ITRF tests): LuPNT's and
  Orekit's EOP tables are different (but both valid) IERS snapshots,
  disagreeing by ~tens of ms in UT1; the ITRF position tolerance scales
  with the orbit radius (`max(100 m, 2e-5 * |r|)`).
- **DE440 distribution differences** (ephemerides, MOON_CI): ~5e-11
  relative (~0.1-0.35 m for the Moon, ~1-7 m for the Sun).
- **IAU pole model vs DE440 principal axes** (Moon body-fixed test):
  Orekit's body-oriented Moon frame is the analytical IAU model, which
  approximates the DE mean-Earth axes to ~3e-5 rad; the tolerance scales as
  3e-4 * |r_moon|. (LuPNT's MOON_PA is validated separately against GMAT's
  DE-kernel Luna frame -- see `cpp/test/gmat/README.md`.)
- **Numerical-integration differences** (j2_propagation): observed ~1e-5 m
  over 1.5 orbits between LuPNT's fixed-step RK4 (1 s) and Orekit's
  adaptive DP853.
- **Force-model formulas** (gravity_acceleration, third_body_acceleration,
  srp_acceleration): matched inputs on both sides (the shared `.cof`/ICGEM
  coefficients for gravity, the stored perturber position for third body, the
  stored `P0` for SRP), so the residual is pure floating-point evaluation of
  the same algorithm: observed 6.2e-15 / 1.9e-18 / 6.7e-24 m/s^2, with
  1e-9 / 1e-12 / 1e-15 m/s^2 tolerances -- far above the observed values yet
  still many orders of magnitude below any real force-model error.

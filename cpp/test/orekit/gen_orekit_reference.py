#!/usr/bin/env python3
"""Generate reference vectors from Orekit for cross-validating LuPNT.

This script is a *developer-only* tool. It is **not** required to build or
test LuPNT normally -- the values it produces are checked into
`data/orekit_reference.json` and consumed by the `*_orekit_reference` Catch2
tests (in this directory) without any dependency on Orekit/Java.

Run it (from the repo root) with the optional `dev` pixi environment, which
provides the JVM-based `orekit` Python wrapper:

    pixi run -e dev python cpp/test/orekit/gen_orekit_reference.py

The first run downloads the small "orekit-data" bundle (leap seconds, EOP,
DE440 ephemerides, ...) into `data/orekit-data-main/` (already covered by
`.gitignore`'s `**/*orekit-data*/` pattern, so it is never committed).

See `cpp/test/orekit/README.md` for the rationale and what each section of
the generated JSON covers.
"""

from __future__ import annotations

import json
import math
import os
import zipfile
from pathlib import Path
from urllib.request import urlretrieve

OUTPUT_PATH = Path(__file__).resolve().parent / "data" / "orekit_reference.json"

# Physical constants -- copied from cpp/lupnt/core/constants.h so the
# comparison below isolates *algorithmic* differences from differences in
# adopted constants.
GM_EARTH = 398600.435507e9  # [m^3/s^2]
GM_MOON = 4902.800118e9  # [m^3/s^2]
R_EARTH = 6378.137e3  # [m]
J2_EARTH = 1.08262668e-3  # [-]

# Epochs used for the scalar (time-scale / sidereal) sections. The list
# deliberately brackets the 2012-06-30 and 2016-12-31 leap seconds so the
# leap-second table lookup is exercised on both sides of an insertion.
#
# Epochs within ~1 day of a leap-second insertion are flagged
# `near_leap_second`: LuPNT's EOP table stores daily UT1-UTC values and
# interpolates *across* the 1 s leap jump (yielding up to ~0.94 s of
# spurious UT1-UTC one hour before an insertion), while Orekit interpolates
# the continuous UT1-TAI. The C++ tests therefore skip the UT1-dependent
# checks (UT1-UTC, GMST, ERA) at flagged epochs; the constant-offset scales
# (TAI/TT/TDB/GPS) are still checked there.
TIME_EPOCHS = [
    # (y, mo, d, h, mi, s, near_leap_second)   TAI-UTC
    (2000, 1, 1, 12, 0, 0.0, False),  # 32 s
    (2005, 8, 10, 7, 15, 0.0, False),  # 32 s
    (2012, 7, 1, 6, 0, 0.0, True),  # 35 s (6 h after the 2012-06-30 leap)
    (2016, 12, 31, 23, 0, 0.0, True),  # 36 s (1 h before the 2016-12-31 leap)
    (2017, 1, 2, 0, 0, 0.0, False),  # 37 s (a full day after)
    (2024, 3, 15, 12, 0, 0.0, False),  # 37 s
]


def epoch_str(y: int, mo: int, d: int, h: int, mi: int, s: float) -> str:
    return f"{y:04d}-{mo:02d}-{d:02d}T{h:02d}:{mi:02d}:{s:09.6f}"


def vec(v3) -> list[float]:
    return [v3.getX(), v3.getY(), v3.getZ()]


def pv_dict(pv, prefix: str) -> dict:
    return {
        f"r_{prefix}": vec(pv.getPosition()),
        f"v_{prefix}": vec(pv.getVelocity()),
    }


def ensure_orekit_data() -> Path:
    """Download/extract the orekit-data bundle next to LUPNT_DATA_PATH if needed."""
    lupnt_data_path = Path(os.environ["LUPNT_DATA_PATH"])
    base_path = lupnt_data_path.parent
    zip_path = base_path / "orekit-data.zip"
    folder_path = base_path / "orekit-data-main"

    if not folder_path.exists():
        if not zip_path.exists():
            url = "https://gitlab.orekit.org/orekit/orekit-data/-/archive/main/orekit-data-main.zip"
            print(f"Downloading {url} -> {zip_path}")
            urlretrieve(url, zip_path)
        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(base_path)

    return folder_path


def main() -> None:
    import orekit
    from orekit import JArray_double

    orekit.initVM()

    from orekit.pyhelpers import setup_orekit_curdir

    setup_orekit_curdir(filename=str(ensure_orekit_data()))

    from org.hipparchus.geometry.euclidean.threed import Vector3D
    from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
    from org.orekit.bodies import CelestialBodyFactory
    from org.orekit.forces.gravity import J2OnlyPerturbation, NewtonianAttraction
    from org.orekit.frames import FramesFactory
    from org.orekit.orbits import CartesianOrbit, KeplerianOrbit, OrbitType, PositionAngleType
    from org.orekit.propagation import SpacecraftState
    from org.orekit.propagation.analytical import KeplerianPropagator
    from org.orekit.propagation.numerical import NumericalPropagator
    from org.orekit.time import AbsoluteDate, TimeScalesFactory
    from org.orekit.utils import IERSConventions, PVCoordinates, PVCoordinatesProvider

    UTC = TimeScalesFactory.getUTC()
    TAI = TimeScalesFactory.getTAI()
    TT = TimeScalesFactory.getTT()
    TDB = TimeScalesFactory.getTDB()
    GPS = TimeScalesFactory.getGPS()
    UT1 = TimeScalesFactory.getUT1(IERSConventions.IERS_2010, True)

    GCRF = FramesFactory.getGCRF()
    EME2000 = FramesFactory.getEME2000()
    ITRF = FramesFactory.getITRF(IERSConventions.IERS_2010, False)

    moon = CelestialBodyFactory.getMoon()
    sun = CelestialBodyFactory.getSun()
    moon_pv_provider = PVCoordinatesProvider.cast_(moon)
    sun_pv_provider = PVCoordinatesProvider.cast_(sun)

    out: dict = {
        "_meta": {
            "description": (
                "Reference vectors generated by Orekit for cross-checking LuPNT's "
                "time-scale conversions, sidereal angles, frame conversions, "
                "Sun/Moon ephemerides, two-body Keplerian orbit propagation, and "
                "J2 acceleration/propagation. See cpp/test/orekit/README.md."
            ),
            "generator": "cpp/test/orekit/gen_orekit_reference.py",
            "orekit_version": "13.1",
            "orekit_data": "https://gitlab.orekit.org/orekit/orekit-data (main)",
            "constants": {
                "GM_EARTH": GM_EARTH,
                "GM_MOON": GM_MOON,
                "R_EARTH": R_EARTH,
                "J2_EARTH": J2_EARTH,
            },
        }
    }

    # ------------------------------------------------------------------
    # Time scales
    # ------------------------------------------------------------------
    time_scales = []
    for y, mo, d, h, mi, s, near_leap in TIME_EPOCHS:
        t0 = AbsoluteDate(y, mo, d, h, mi, s, UTC)
        time_scales.append(
            {
                "epoch_utc": epoch_str(y, mo, d, h, mi, s),
                "near_leap_second": near_leap,
                "tai_minus_utc": -UTC.offsetFromTAI(t0),
                "tt_minus_tai": TT.offsetFromTAI(t0),
                "tdb_minus_tai": TDB.offsetFromTAI(t0),
                "gps_minus_tai": GPS.offsetFromTAI(t0),
                "ut1_minus_utc": UT1.offsetFromTAI(t0) - UTC.offsetFromTAI(t0),
            }
        )
    out["time_scales"] = time_scales

    # ------------------------------------------------------------------
    # Sidereal / Earth-rotation angles: GMST (IERS 1996 conventions ==
    # the classic IAU-82 formula LuPNT implements) and the IAU-2000
    # Earth Rotation Angle, both functions of UT1.
    # ------------------------------------------------------------------
    try:
        gmst_func = IERSConventions.IERS_1996.getGMSTFunction(UT1)
        era_func = IERSConventions.IERS_2010.getEarthOrientationAngleFunction(UT1)
        sidereal = []
        for y, mo, d, h, mi, s, near_leap in TIME_EPOCHS:
            t0 = AbsoluteDate(y, mo, d, h, mi, s, UTC)
            sidereal.append(
                {
                    "epoch_utc": epoch_str(y, mo, d, h, mi, s),
                    "near_leap_second": near_leap,
                    "gmst": gmst_func.value(t0) % (2 * math.pi),
                    "era": era_func.value(t0) % (2 * math.pi),
                }
            )
        out["sidereal"] = sidereal
    except Exception as exc:  # pragma: no cover - API availability guard
        print(f"WARNING: skipping sidereal section ({exc})")

    # ------------------------------------------------------------------
    # Frame conversions: GCRF <-> EME2000 (fixed frame bias) and
    # GCRF <-> ITRF (full EOP-dependent rotation), across orbital regimes
    # from LEO to trans-lunar distances.
    # ------------------------------------------------------------------
    frames = []
    sample_states = [
        # LEO
        (Vector3D(7000e3, 1000e3, 2000e3), Vector3D(-1.0e3, 7.0e3, 1.5e3)),
        # MEO (GPS-like radius)
        (Vector3D(-4000e3, 30000e3, 6000e3), Vector3D(-2.0e3, -0.4e3, 1.0e3)),
        # GEO
        (Vector3D(40000e3, -13000e3, 50e3), Vector3D(0.95e3, 2.92e3, 0.01e3)),
        # Trans-lunar distance
        (Vector3D(250000e3, 200000e3, 100000e3), Vector3D(-0.5e3, 0.8e3, 0.3e3)),
    ]
    for y, mo, d, h, mi, s in [
        (2010, 6, 1, 0, 0, 0.0),
        (2024, 3, 15, 12, 0, 0.0),
        (2025, 9, 20, 6, 30, 0.0),
    ]:
        t0 = AbsoluteDate(y, mo, d, h, mi, s, UTC)
        for pos, vel in sample_states:
            pv = PVCoordinates(pos, vel)
            pv_eme = GCRF.getTransformTo(EME2000, t0).transformPVCoordinates(pv)
            pv_itrf = GCRF.getTransformTo(ITRF, t0).transformPVCoordinates(pv)
            frames.append(
                {
                    "epoch_utc": epoch_str(y, mo, d, h, mi, s),
                    "r_gcrf": vec(pos),
                    "v_gcrf": vec(vel),
                    **pv_dict(pv_eme, "eme2000"),
                    **pv_dict(pv_itrf, "itrf"),
                }
            )
    out["frames"] = frames

    # ------------------------------------------------------------------
    # Moon-centered frames: GCRF -> MOON_CI (Moon-centered inertial,
    # GCRF-aligned axes: a pure translation by the DE440 Earth->Moon
    # vector) and GCRF -> Moon body-fixed (Orekit's IAU-pole
    # body-oriented frame; LuPNT's MOON_PA/MOON_ME come from the DE440
    # principal-axes kernel, so this comparison is approximate -- see
    # the README for the tolerance rationale).
    # ------------------------------------------------------------------
    moon_body_frame = moon.getBodyOrientedFrame()
    moon_frames = []
    for y, mo, d, h, mi, s in [
        (2024, 3, 15, 12, 0, 0.0),
        (2025, 9, 20, 6, 30, 0.0),
    ]:
        t0 = AbsoluteDate(y, mo, d, h, mi, s, UTC)
        moon_pv = moon_pv_provider.getPVCoordinates(t0, GCRF)
        sample_pvs = [
            # Low lunar orbiter: ~2000 km from the Moon's center
            PVCoordinates(
                moon_pv.getPosition().add(Vector3D(1500e3, 1000e3, 800e3)),
                moon_pv.getVelocity().add(Vector3D(0.5e3, -1.2e3, 0.3e3)),
            ),
            # ELFO-like: a few thousand km from the Moon
            PVCoordinates(
                moon_pv.getPosition().add(Vector3D(-4000e3, 2500e3, 5000e3)),
                moon_pv.getVelocity().add(Vector3D(-0.3e3, 0.6e3, -0.2e3)),
            ),
        ]
        for pv in sample_pvs:
            pv_moon_ci = PVCoordinates(
                pv.getPosition().subtract(moon_pv.getPosition()),
                pv.getVelocity().subtract(moon_pv.getVelocity()),
            )
            pv_moon_fixed = GCRF.getTransformTo(moon_body_frame, t0).transformPVCoordinates(pv)
            moon_frames.append(
                {
                    "epoch_utc": epoch_str(y, mo, d, h, mi, s),
                    "r_gcrf": vec(pv.getPosition()),
                    "v_gcrf": vec(pv.getVelocity()),
                    **pv_dict(pv_moon_ci, "moon_ci"),
                    **pv_dict(pv_moon_fixed, "moon_fixed_iau"),
                }
            )
    out["moon_frames"] = moon_frames

    # ------------------------------------------------------------------
    # Sun/Moon ephemerides (DE440): position/velocity of the Moon and the
    # Sun relative to the Earth, in GCRF axes.
    # ------------------------------------------------------------------
    ephemerides = []
    for y, mo, d, h, mi, s in [
        (2010, 6, 1, 0, 0, 0.0),
        (2024, 3, 15, 12, 0, 0.0),
        (2025, 9, 20, 6, 30, 0.0),
    ]:
        t0 = AbsoluteDate(y, mo, d, h, mi, s, UTC)
        for name, body in [("MOON", moon_pv_provider), ("SUN", sun_pv_provider)]:
            pv = body.getPVCoordinates(t0, GCRF)
            ephemerides.append(
                {
                    "epoch_utc": epoch_str(y, mo, d, h, mi, s),
                    "target": name,
                    "center": "EARTH",
                    "frame": "GCRF",
                    **pv_dict(pv, "gcrf"),
                }
            )
    out["ephemerides"] = ephemerides

    # ------------------------------------------------------------------
    # Two-body Keplerian propagation + COE <-> Cartesian + anomalies,
    # for Earth orbits (MEO-HEO/LEO/GEO/Molniya) and lunar orbits
    # (low lunar orbit and ELFO, with GM_MOON). The `gm` field is read by
    # the C++ test, so Earth and Moon cases share the same schema.
    # ------------------------------------------------------------------
    kepler = []
    t0 = AbsoluteDate(2024, 3, 15, 12, 0, 0.0, UTC)
    coe_cases = [
        # gm, body, a [m], e, i [deg], raan [deg], argp [deg], M0 [deg]
        (GM_EARTH, "EARTH", 24396000.0, 0.10, 30.0, 50.0, 60.0, 10.0),
        (GM_EARTH, "EARTH", 7000000.0, 0.001, 98.0, 120.0, 0.0, 0.0),
        (GM_EARTH, "EARTH", 42164000.0, 0.0005, 0.05, 10.0, 20.0, 30.0),  # GEO
        (GM_EARTH, "EARTH", 26562000.0, 0.74, 63.4, 200.0, 270.0, 45.0),  # Molniya
        (GM_MOON, "MOON", 2038000.0, 0.01, 92.0, 40.0, 270.0, 20.0),  # LLO
        (GM_MOON, "MOON", 6541400.0, 0.60, 56.2, 0.0, 90.0, 0.0),  # ELFO
    ]
    for gm, body, a, e, i_deg, raan_deg, argp_deg, m0_deg in coe_cases:
        orbit0 = KeplerianOrbit(
            a,
            e,
            math.radians(i_deg),
            math.radians(argp_deg),
            math.radians(raan_deg),
            math.radians(m0_deg),
            PositionAngleType.MEAN,
            GCRF,
            t0,
            gm,
        )
        pv0 = orbit0.getPVCoordinates()
        period = orbit0.getKeplerianPeriod()
        prop = KeplerianPropagator(orbit0)

        case = {
            "gm": gm,
            "body": body,
            "epoch_utc": "2024-03-15T12:00:00.000000",
            "coe0": [
                a,
                e,
                math.radians(i_deg),
                math.radians(raan_deg),
                math.radians(argp_deg),
                math.radians(m0_deg),
            ],
            "cart0": vec(pv0.getPosition()) + vec(pv0.getVelocity()),
            "ecc_anomaly0": orbit0.getEccentricAnomaly(),
            "true_anomaly0": orbit0.getTrueAnomaly(),
            "period": period,
            "propagated": [],
        }

        for frac in [0.25, 0.5, 0.75, 1.5]:
            dt = period * frac
            state = prop.propagate(t0.shiftedBy(dt))
            pv = state.getPVCoordinates(GCRF)
            orbit = KeplerianOrbit(pv, GCRF, state.getDate(), gm)
            case["propagated"].append(
                {
                    "dt": dt,
                    "cart": vec(pv.getPosition()) + vec(pv.getVelocity()),
                    "coe": [
                        orbit.getA(),
                        orbit.getE(),
                        orbit.getI(),
                        orbit.getRightAscensionOfAscendingNode(),
                        orbit.getPerigeeArgument(),
                        orbit.getMeanAnomaly(),
                    ],
                    "ecc_anomaly": orbit.getEccentricAnomaly(),
                    "true_anomaly": orbit.getTrueAnomaly(),
                }
            )

        kepler.append(case)
    out["kepler"] = kepler

    # ------------------------------------------------------------------
    # J2 acceleration (formula-level check, frame-independent)
    # ------------------------------------------------------------------
    j2 = []
    positions = [
        (7000e3, 0.0, 0.0),
        (7000e3, 1000e3, 2000e3),
        (-4000e3, 5000e3, -3000e3),
        (24396000.0 * 0.9, 0.0, 6000e3),
        (0.0, 0.0, 8000e3),
    ]
    for x, y, z in positions:
        p = Vector3D(x, y, z)
        acc = J2OnlyPerturbation.computeAccelerationInJ2Frame(p, GM_EARTH, R_EARTH, J2_EARTH)
        j2.append({"r": [x, y, z], "a_j2": vec(acc)})
    out["j2_acceleration"] = j2

    # ------------------------------------------------------------------
    # J2-perturbed numerical propagation. The J2OnlyPerturbation force
    # model is evaluated with **GCRF** (inertial) as the J2 frame, i.e.
    # the oblateness symmetry axis is the inertial Z axis -- the *same*
    # modeling choice as LuPNT's JToCartTwoBodyDynamics, so unlike the
    # GMAT j2_propagation comparison (body-fixed J2), this one isolates
    # pure numerical-integration differences and should agree tightly.
    # ------------------------------------------------------------------
    def make_j2_propagator(orbit0):
        tol = NumericalPropagator.tolerances(1e-6, orbit0, OrbitType.CARTESIAN)
        integrator = DormandPrince853Integrator(
            1e-6, 300.0, JArray_double.cast_(tol[0]), JArray_double.cast_(tol[1])
        )
        propagator = NumericalPropagator(integrator)
        propagator.setOrbitType(OrbitType.CARTESIAN)
        propagator.addForceModel(NewtonianAttraction(GM_EARTH))
        propagator.addForceModel(J2OnlyPerturbation(GM_EARTH, R_EARTH, J2_EARTH, GCRF))
        propagator.setInitialState(SpacecraftState(orbit0))
        return propagator

    j2_propagation = []
    j2_cases = [
        # a [m], e, i [deg], raan [deg], argp [deg], M0 [deg]
        (24396000.0, 0.10, 30.0, 50.0, 60.0, 10.0),
        (7000000.0, 0.001, 98.0, 120.0, 0.0, 0.0),
        (26562000.0, 0.74, 63.4, 200.0, 270.0, 45.0),  # Molniya
    ]
    for a, e, i_deg, raan_deg, argp_deg, m0_deg in j2_cases:
        kep0 = KeplerianOrbit(
            a,
            e,
            math.radians(i_deg),
            math.radians(argp_deg),
            math.radians(raan_deg),
            math.radians(m0_deg),
            PositionAngleType.MEAN,
            GCRF,
            t0,
            GM_EARTH,
        )
        orbit0 = CartesianOrbit(kep0.getPVCoordinates(), GCRF, t0, GM_EARTH)
        period = kep0.getKeplerianPeriod()

        case = {
            "gm": GM_EARTH,
            "j2": J2_EARTH,
            "r_earth": R_EARTH,
            "epoch_utc": "2024-03-15T12:00:00.000000",
            "cart0": vec(orbit0.getPVCoordinates().getPosition())
            + vec(orbit0.getPVCoordinates().getVelocity()),
            "period": period,
            "propagated": [],
        }
        for frac in [0.25, 0.5, 0.75, 1.5]:
            dt = period * frac
            propagator = make_j2_propagator(orbit0)
            state = propagator.propagate(t0.shiftedBy(dt))
            pv = state.getPVCoordinates(GCRF)
            case["propagated"].append(
                {"dt": dt, "cart": vec(pv.getPosition()) + vec(pv.getVelocity())}
            )
        j2_propagation.append(case)
    out["j2_propagation"] = j2_propagation

    OUTPUT_PATH.write_text(json.dumps(out, indent=2) + "\n")
    print(f"Wrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Generate reference vectors from GMAT for cross-validating LuPNT.

This script is a *developer-only* tool. It is **not** required to build or
test LuPNT normally -- the values it produces are checked into
`data/gmat_reference.json` and consumed by the `*_gmat_reference` Catch2
tests (in this directory) without any dependency on GMAT.

Unlike Orekit (an optional pixi/JVM dependency), GMAT is a separately
installed desktop application -- see
https://sourceforge.net/projects/gmat/ -- and is **not** distributed via
pixi/conda. To run this script you need a local GMAT R2026a (or compatible)
installation and must point it at the `GmatConsole` executable, either by
setting the `GMAT_CONSOLE` environment variable to the full path of the
executable, or `GMAT_HOME`/`GMAT_ROOT_DIR` to the GMAT install directory
(the script will look for `bin/GmatConsole*` under it). Common default
install locations are also searched automatically.

Usage (from the repo root; only the Python standard library is used):

    GMAT_CONSOLE=/path/to/GMAT/R2026a/bin/GmatConsole-R2026a \\
        python cpp/test/gmat/gen_gmat_reference.py

The script *builds* a GMAT script from the case tables below (substituting
absolute paths for the custom potential file and the report outputs), runs
it through `GmatConsole -r <script>`, parses the resulting `ReportFile`
outputs, and overwrites `cpp/test/gmat/data/gmat_reference.json`. A copy of
the generated script (with machine-specific paths replaced by placeholder
tokens) is saved to `cpp/test/gmat/data/gmat_reference.script` for
documentation.

IMPORTANT (GMAT Propagate workaround): each generated `Propagate` command
propagates exactly ONE spacecraft with ONE propagator. Combining two
spacecraft/propagators in a single synced
`Propagate Prop1(SatA) Prop2(SatB) {SatA.ElapsedSecs = X};` command was
found (originally on R2022a) to desync the second spacecraft's force model
(J2Sat ended up ~60x further from the two-body solution than the expected J2
perturbation); the one-spacecraft-per-command form is kept regardless.

See `cpp/test/gmat/README.md` for the JSON schema, the GMAT-specific
limitations (no UT1/GPS time scales, custom potential file, DE440 via SPICE
since GMAT ships only up to DE424, etc.), and the tolerance rationale used by
the `*_gmat_reference` tests.
"""

from __future__ import annotations

import glob
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
COF_PATH = HERE / "lupnt_j2_earth.cof"
OUTPUT_PATH = HERE / "data" / "gmat_reference.json"
SCRIPT_DOC_PATH = HERE / "data" / "gmat_reference.script"

# LuPNT's own DE440 kernel. GMAT ships only up to DE424, so it is pointed at
# LuPNT's DE440 via SPICE -- the two then use the *same* planetary ephemeris,
# removing the former DE421-vs-DE440 offset that dominated the Moon-frame
# residuals. Override with LUPNT_DE440 if the data lives elsewhere.
DE440_PATH = os.environ.get(
    "LUPNT_DE440",
    str(Path(__file__).resolve().parents[3] / "data" / "LuPNT_data" / "ephemeris" / "de440.bsp"),
)

# Physical constants -- copied from cpp/lupnt/core/constants.h (GM/R/J2 for
# the Earth are also baked into lupnt_j2_earth.cof) so the comparison
# isolates *algorithmic* differences from differences in adopted constants.
GM_EARTH = 398600.435507e9  # [m^3/s^2]
GM_MOON = 4902.800118e9  # [m^3/s^2]
R_EARTH = 6378.137e3  # [m]
J2_EARTH = 1.08262668e-3  # [-]

GMAT_SUCCESS_MARKER = "GMAT Integration test (Console version) successful"

MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

# ----------------------------------------------------------------------------
# Case tables (kept in sync with cpp/test/orekit/gen_orekit_reference.py so
# the two fixtures share input values where the comparison overlaps).
# ----------------------------------------------------------------------------

# (y, mo, d, h, mi, s) -- same epochs as the Orekit fixture, bracketing the
# 2012-06-30 and 2016-12-31 leap seconds. GMAT exposes no UT1, so there is
# no near-leap-second concern here (TAI/TT/TDB are leap-jump-free).
TIME_EPOCHS = [
    (2000, 1, 1, 12, 0, 0.0),
    (2005, 8, 10, 7, 15, 0.0),
    (2012, 7, 1, 6, 0, 0.0),
    (2016, 12, 31, 23, 0, 0.0),
    (2017, 1, 2, 0, 0, 0.0),
    (2024, 3, 15, 12, 0, 0.0),
]

FRAME_EPOCHS = [
    (2010, 6, 1, 0, 0, 0.0),
    (2024, 3, 15, 12, 0, 0.0),
    (2025, 9, 20, 6, 30, 0.0),
]

# (r [km], v [km/s]) in Earth ICRF axes -- LEO, MEO, GEO, trans-lunar.
FRAME_STATES = [
    ((7000.0, 1000.0, 2000.0), (-1.0, 7.0, 1.5)),
    ((-4000.0, 30000.0, 6000.0), (-2.0, -0.4, 1.0)),
    ((40000.0, -13000.0, 50.0), (0.95, 2.92, 0.01)),
    ((250000.0, 200000.0, 100000.0), (-0.5, 0.8, 0.3)),
]

MOON_EPOCHS = [
    (2024, 3, 15, 12, 0, 0.0),
    (2025, 9, 20, 6, 30, 0.0),
]

# (r [km], v [km/s]) in Luna-centered ICRF axes -- same lunar-centric
# offsets as the Orekit moon_frames cases.
MOON_STATES = [
    ((1500.0, 1000.0, 800.0), (0.5, -1.2, 0.3)),
    ((-4000.0, 2500.0, 5000.0), (-0.3, 0.6, -0.2)),
]

KEPLER_EPOCH = (2024, 3, 15, 12, 0, 0.0)

# Earth orbits: a [m], e, i, raan, aop, M0 [deg]; `j2` adds a twin
# spacecraft propagated with the PointMass+J2 force model.
EARTH_KEPLER_CASES = [
    dict(
        sat="KepSat1",
        j2sat="J2Sat1",
        a=24396000.0,
        e=0.10,
        i=30.0,
        raan=50.0,
        aop=60.0,
        m0=10.0,
        j2=True,
    ),
    dict(
        sat="KepSat2",
        j2sat="J2Sat2",
        a=7000000.0,
        e=0.001,
        i=98.0,
        raan=120.0,
        aop=0.0,
        m0=0.0,
        j2=True,
    ),
    dict(
        sat="KepSat3",
        j2sat=None,
        a=42164000.0,
        e=0.0005,
        i=0.05,
        raan=10.0,
        aop=20.0,
        m0=30.0,
        j2=False,
    ),  # GEO
    dict(
        sat="KepSat4",
        j2sat=None,
        a=26562000.0,
        e=0.74,
        i=63.4,
        raan=200.0,
        aop=270.0,
        m0=45.0,
        j2=False,
    ),  # Molniya
]

# Lunar orbits (Luna point-mass only): LLO and ELFO, matching the Orekit
# lunar kepler cases.
MOON_KEPLER_CASES = [
    dict(sat="LunaKepSat1", a=2038000.0, e=0.01, i=92.0, raan=40.0, aop=270.0, m0=20.0),
    dict(sat="LunaKepSat2", a=6541400.0, e=0.60, i=56.2, raan=0.0, aop=90.0, m0=0.0),
]

DT_FRACTIONS = [0.25, 0.5, 0.75, 1.5]


# ----------------------------------------------------------------------------
# Locating the GmatConsole executable
# ----------------------------------------------------------------------------
def find_gmat_console() -> Path:
    env = os.environ.get("GMAT_CONSOLE")
    if env:
        p = Path(env).expanduser()
        if p.is_file():
            return p
        raise FileNotFoundError(f"GMAT_CONSOLE={env!r} does not point to a file")

    candidates: list[str] = []
    home = os.environ.get("GMAT_HOME") or os.environ.get("GMAT_ROOT_DIR")
    if home:
        candidates += sorted(glob.glob(str(Path(home).expanduser() / "bin" / "GmatConsole*")))

    # Common default install locations (macOS, Linux, user-local).
    candidates += sorted(glob.glob("/Applications/GMAT_*/bin/GmatConsole-*"))
    candidates += sorted(glob.glob("/usr/local/GMAT*/bin/GmatConsole*"))
    candidates += sorted(glob.glob("/opt/GMAT*/bin/GmatConsole*"))
    candidates += sorted(glob.glob(str(Path.home() / "GMAT*" / "bin" / "GmatConsole*")))

    for c in candidates:
        p = Path(c)
        if p.is_file() and os.access(p, os.X_OK):
            return p

    raise FileNotFoundError(
        "Could not find a GmatConsole executable. Set the GMAT_CONSOLE "
        "environment variable to its full path (e.g. "
        "/path/to/GMAT/R2026a/bin/GmatConsole-R2026a), or GMAT_HOME to "
        "the GMAT install directory."
    )


# ----------------------------------------------------------------------------
# Kepler-equation helpers (pure Python, no GMAT/LuPNT dependency)
# ----------------------------------------------------------------------------
def mean_to_ecc_anomaly(M: float, e: float, tol: float = 1e-15, max_iter: int = 100) -> float:
    M = ((M + math.pi) % (2 * math.pi)) - math.pi
    E = M if e < 0.8 else math.pi
    for _ in range(max_iter):
        f = E - e * math.sin(E) - M
        fp = 1 - e * math.cos(E)
        dE = -f / fp
        E += dE
        if abs(dE) < tol:
            break
    return E


def ecc_to_true_anomaly(E: float, e: float) -> float:
    return 2 * math.atan2(math.sqrt(1 + e) * math.sin(E / 2), math.sqrt(1 - e) * math.cos(E / 2))


def true_to_ecc_anomaly(nu: float, e: float) -> float:
    return 2 * math.atan2(math.sqrt(1 - e) * math.sin(nu / 2), math.sqrt(1 + e) * math.cos(nu / 2))


# ----------------------------------------------------------------------------
# GMAT script generation
# ----------------------------------------------------------------------------
def gmat_epoch(epoch: tuple) -> str:
    y, mo, d, h, mi, s = epoch
    return f"{d:02d} {MONTHS[mo - 1]} {y:04d} {h:02d}:{mi:02d}:{s:06.3f}"


def epoch_str(epoch: tuple) -> str:
    y, mo, d, h, mi, s = epoch
    return f"{y:04d}-{mo:02d}-{d:02d}T{h:02d}:{mi:02d}:{s:09.6f}"


def cart_spacecraft(name: str, epoch: tuple, cs: str, r_km, v_km) -> str:
    return f"""Create Spacecraft {name};
GMAT {name}.DateFormat = UTCGregorian;
GMAT {name}.Epoch = '{gmat_epoch(epoch)}';
GMAT {name}.CoordinateSystem = {cs};
GMAT {name}.DisplayStateType = Cartesian;
GMAT {name}.X = {r_km[0]!r};
GMAT {name}.Y = {r_km[1]!r};
GMAT {name}.Z = {r_km[2]!r};
GMAT {name}.VX = {v_km[0]!r};
GMAT {name}.VY = {v_km[1]!r};
GMAT {name}.VZ = {v_km[2]!r};
"""


def kep_spacecraft(name: str, epoch: tuple, cs: str, case: dict) -> str:
    # GMAT's Keplerian display state only accepts True Anomaly as the 6th
    # element, so convert the case's Mean Anomaly.
    E0 = mean_to_ecc_anomaly(math.radians(case["m0"]), case["e"])
    ta_deg = math.degrees(ecc_to_true_anomaly(E0, case["e"])) % 360.0
    return f"""Create Spacecraft {name};
GMAT {name}.DateFormat = UTCGregorian;
GMAT {name}.Epoch = '{gmat_epoch(epoch)}';
GMAT {name}.CoordinateSystem = {cs};
GMAT {name}.DisplayStateType = Keplerian;
GMAT {name}.SMA = {case['a'] / 1e3!r};
GMAT {name}.ECC = {case['e']!r};
GMAT {name}.INC = {case['i']!r};
GMAT {name}.RAAN = {case['raan']!r};
GMAT {name}.AOP = {case['aop']!r};
GMAT {name}.TA = {ta_deg!r};
"""


def report_file(name: str, path: str) -> str:
    return f"""Create ReportFile {name};
GMAT {name}.Filename = '{path}';
GMAT {name}.Precision = 16;
GMAT {name}.WriteHeaders = true;
GMAT {name}.ColumnWidth = 24;
GMAT {name}.SolverIterations = None;
"""


def cs_cart_params(sat: str, cs: str) -> str:
    return " ".join(f"{sat}.{cs}.{c}" for c in ("X", "Y", "Z", "VX", "VY", "VZ"))


def earth_kep_params(sat: str) -> str:
    return " ".join(
        f"{sat}.{p}"
        for p in ("SMA", "ECC", "INC", "RAAN", "AOP", "MA", "TA", "X", "Y", "Z", "VX", "VY", "VZ")
    )


def luna_kep_params(sat: str) -> str:
    # SMA/ECC/MA/TA are origin-dependent parameters (qualified with the
    # central body), INC/RAAN/AOP and the Cartesian state are
    # coordinate-system-dependent (qualified with the coordinate system).
    p = [
        f"{sat}.Luna.SMA",
        f"{sat}.Luna.ECC",
        f"{sat}.LunaICRFCS.INC",
        f"{sat}.LunaICRFCS.RAAN",
        f"{sat}.LunaICRFCS.AOP",
        f"{sat}.Luna.MA",
        f"{sat}.Luna.TA",
    ]
    p += [f"{sat}.LunaICRFCS.{c}" for c in ("X", "Y", "Z", "VX", "VY", "VZ")]
    return " ".join(p)


def case_period(case: dict, gm: float) -> float:
    return 2 * math.pi * math.sqrt(case["a"] ** 3 / gm)


def segment_dts(period: float) -> list[float]:
    """Per-segment ElapsedSecs so the cumulative time hits DT_FRACTIONS."""
    cum = [period * f for f in DT_FRACTIONS]
    return [cum[0]] + [cum[i] - cum[i - 1] for i in range(1, len(cum))]


def build_script(cof_path: str, rf_paths: dict[str, str]) -> str:
    s: list[str] = []
    s.append(
        f"""%----------------------------------------------------------------------------
% LuPNT cross-validation reference script (GMAT R2026a)
%
% GENERATED by cpp/test/gmat/gen_gmat_reference.py -- do not edit by hand;
% edit the case tables in the generator instead.
%
% Sections (one ReportFile each, parsed into gmat_reference.json):
%   A. Time scale offsets (TAI/TT/TDB - UTC)                 -> rfTime
%   B. Earth frame conversions (ICRF/MJ2000Eq/BodyFixed)     -> rfFrames
%   C. Moon-centered frames (LunaICRF/LunaFixed)             -> rfMoon
%   D. Earth two-body Keplerian (+J2) propagation            -> rfKep1..4
%   E. Lunar two-body Keplerian propagation                  -> rfLuna1..2
%
% IMPORTANT: each `Propagate` command below propagates exactly ONE
% spacecraft with ONE propagator. Combining two spacecraft/propagators in a
% single synced `Propagate Prop1(SatA) Prop2(SatB) {{SatA.ElapsedSecs = X}};`
% command was found (originally on R2022a) to desync the second spacecraft's
% force model. Always use separate sequential `Propagate` commands.
%----------------------------------------------------------------------------

% Use LuPNT's own DE440 kernel via SPICE, so both sides share one planetary
% ephemeris (GMAT ships only up to DE424). This removes the former DE421-vs-DE440
% offset that dominated the Moon-frame residuals.
GMAT SolarSystem.EphemerisSource = 'SPICE';
GMAT SolarSystem.SPKFilename = '{DE440_PATH}';

% Use the exact GM/R values adopted by LuPNT (cpp/lupnt/core/constants.h) so
% the comparison isolates *algorithmic* differences from differences in
% adopted constants.
GMAT Earth.Mu = {GM_EARTH / 1e9!r};
GMAT Earth.EquatorialRadius = {R_EARTH / 1e3!r};
GMAT Luna.Mu = {GM_MOON / 1e9!r};

%----------------------------------------------------------------------------
% Coordinate systems
%----------------------------------------------------------------------------
Create CoordinateSystem EarthICRFCS;
GMAT EarthICRFCS.Origin = Earth;
GMAT EarthICRFCS.Axes = ICRF;

Create CoordinateSystem EarthFixedCS;
GMAT EarthFixedCS.Origin = Earth;
GMAT EarthFixedCS.Axes = BodyFixed;

Create CoordinateSystem LunaICRFCS;
GMAT LunaICRFCS.Origin = Luna;
GMAT LunaICRFCS.Axes = ICRF;

Create CoordinateSystem LunaFixedCS;
GMAT LunaFixedCS.Origin = Luna;
GMAT LunaFixedCS.Axes = BodyFixed;
"""
    )

    # --- Section A: time-scale spacecraft ---------------------------------
    s.append("%---------------------------- Section A: time scales\n")
    for k, ep in enumerate(TIME_EPOCHS, 1):
        s.append(
            f"""Create Spacecraft TimeSat{k};
GMAT TimeSat{k}.DateFormat = UTCGregorian;
GMAT TimeSat{k}.Epoch = '{gmat_epoch(ep)}';
GMAT TimeSat{k}.CoordinateSystem = EarthMJ2000Eq;
"""
        )
    s.append(report_file("rfTime", rf_paths["rfTime"]))

    # --- Section B: Earth frame spacecraft --------------------------------
    s.append("%---------------------------- Section B: Earth frames\n")
    frame_sats = []
    for ei, ep in enumerate(FRAME_EPOCHS, 1):
        for si, (r, v) in enumerate(FRAME_STATES, 1):
            name = f"Sat_e{ei}_s{si}"
            frame_sats.append(name)
            s.append(cart_spacecraft(name, ep, "EarthICRFCS", r, v))
    s.append(report_file("rfFrames", rf_paths["rfFrames"]))

    # --- Section C: Moon frame spacecraft ---------------------------------
    s.append("%---------------------------- Section C: Moon frames\n")
    moon_sats = []
    for ei, ep in enumerate(MOON_EPOCHS, 1):
        for si, (r, v) in enumerate(MOON_STATES, 1):
            name = f"MoonSat_e{ei}_s{si}"
            moon_sats.append(name)
            s.append(cart_spacecraft(name, ep, "LunaICRFCS", r, v))
    s.append(report_file("rfMoon", rf_paths["rfMoon"]))

    # --- Sections D/E: propagation spacecraft, force models, propagators ---
    s.append("%---------------------------- Sections D/E: propagation\n")
    for case in EARTH_KEPLER_CASES:
        s.append(kep_spacecraft(case["sat"], KEPLER_EPOCH, "EarthMJ2000Eq", case))
        if case["j2"]:
            s.append(kep_spacecraft(case["j2sat"], KEPLER_EPOCH, "EarthMJ2000Eq", case))
    for case in MOON_KEPLER_CASES:
        s.append(kep_spacecraft(case["sat"], KEPLER_EPOCH, "LunaICRFCS", case))

    s.append(
        f"""% PointMass-only force model (Degree=0, Order=0 in the custom potential
% file): matches LuPNT's analytical two-body KeplerianDynamics.
Create ForceModel TwoBodyFM;
GMAT TwoBodyFM.CentralBody = Earth;
GMAT TwoBodyFM.PrimaryBodies = {{Earth}};
GMAT TwoBodyFM.Drag = None;
GMAT TwoBodyFM.SRP = Off;
GMAT TwoBodyFM.GravityField.Earth.Degree = 0;
GMAT TwoBodyFM.GravityField.Earth.Order = 0;
GMAT TwoBodyFM.GravityField.Earth.PotentialFile = '{cof_path}';

% PointMass + J2-only force model (Degree=2, Order=0): matches LuPNT's
% JToCartTwoBodyDynamics(GM_EARTH, J2_EARTH, R_EARTH) modulo the J2 axis
% (GMAT: body-fixed; LuPNT: inertial Z -- see README).
Create ForceModel J2FM;
GMAT J2FM.CentralBody = Earth;
GMAT J2FM.PrimaryBodies = {{Earth}};
GMAT J2FM.Drag = None;
GMAT J2FM.SRP = Off;
GMAT J2FM.GravityField.Earth.Degree = 2;
GMAT J2FM.GravityField.Earth.Order = 0;
GMAT J2FM.GravityField.Earth.PotentialFile = '{cof_path}';

% Luna point-mass-only force model (uses the overridden Luna.Mu).
Create ForceModel LunaTwoBodyFM;
GMAT LunaTwoBodyFM.CentralBody = Luna;
GMAT LunaTwoBodyFM.PointMasses = {{Luna}};
GMAT LunaTwoBodyFM.Drag = None;
GMAT LunaTwoBodyFM.SRP = Off;

Create Propagator TwoBodyProp;
GMAT TwoBodyProp.FM = TwoBodyFM;
GMAT TwoBodyProp.Type = RungeKutta89;
GMAT TwoBodyProp.InitialStepSize = 60;
GMAT TwoBodyProp.Accuracy = 1e-13;
GMAT TwoBodyProp.MinStep = 0;
GMAT TwoBodyProp.MaxStep = 600;

Create Propagator J2Prop;
GMAT J2Prop.FM = J2FM;
GMAT J2Prop.Type = RungeKutta89;
GMAT J2Prop.InitialStepSize = 60;
GMAT J2Prop.Accuracy = 1e-13;
GMAT J2Prop.MinStep = 0;
GMAT J2Prop.MaxStep = 600;

Create Propagator LunaTwoBodyProp;
GMAT LunaTwoBodyProp.FM = LunaTwoBodyFM;
GMAT LunaTwoBodyProp.Type = RungeKutta89;
GMAT LunaTwoBodyProp.InitialStepSize = 60;
GMAT LunaTwoBodyProp.Accuracy = 1e-13;
GMAT LunaTwoBodyProp.MinStep = 0;
GMAT LunaTwoBodyProp.MaxStep = 600;
"""
    )
    for k in range(1, len(EARTH_KEPLER_CASES) + 1):
        s.append(report_file(f"rfKep{k}", rf_paths[f"rfKep{k}"]))
    for k in range(1, len(MOON_KEPLER_CASES) + 1):
        s.append(report_file(f"rfLuna{k}", rf_paths[f"rfLuna{k}"]))

    # --- Mission sequence ---------------------------------------------------
    s.append(
        """%----------------------------------------------------------------------------
% Mission sequence
%----------------------------------------------------------------------------
BeginMissionSequence;

% --- Section A: time scales
"""
    )
    for k in range(1, len(TIME_EPOCHS) + 1):
        params = " ".join(
            f"TimeSat{k}.{p}"
            for p in ("UTCModJulian", "TAIModJulian", "TTModJulian", "TDBModJulian")
        )
        s.append(f"Report rfTime {params};\n")

    s.append("\n% --- Section B: Earth frames\n")
    for name in frame_sats:
        params = " ".join(
            [
                cs_cart_params(name, "EarthICRFCS"),
                cs_cart_params(name, "EarthMJ2000Eq"),
                cs_cart_params(name, "EarthFixedCS"),
            ]
        )
        s.append(f"Report rfFrames {params};\n")

    s.append("\n% --- Section C: Moon frames\n")
    for name in moon_sats:
        params = " ".join(
            [
                cs_cart_params(name, "EarthICRFCS"),
                cs_cart_params(name, "LunaICRFCS"),
                cs_cart_params(name, "LunaFixedCS"),
            ]
        )
        s.append(f"Report rfMoon {params};\n")

    s.append("\n% --- Section D: Earth Keplerian (+J2) propagation\n")
    for k, case in enumerate(EARTH_KEPLER_CASES, 1):
        sat, j2sat = case["sat"], case["j2sat"]
        params = earth_kep_params(sat)
        if case["j2"]:
            params += " " + earth_kep_params(j2sat)
        period = case_period(case, GM_EARTH)
        s.append(f"\n% case {k}: a = {case['a'] / 1e3} km, e = {case['e']}\n")
        s.append(f"Report rfKep{k} {params};\n")
        for dt in segment_dts(period):
            s.append(f"Propagate TwoBodyProp({sat}) {{{sat}.ElapsedSecs = {dt!r}}};\n")
            if case["j2"]:
                s.append(f"Propagate J2Prop({j2sat}) {{{j2sat}.ElapsedSecs = {dt!r}}};\n")
            s.append(f"Report rfKep{k} {params};\n")

    s.append("\n% --- Section E: lunar Keplerian propagation\n")
    for k, case in enumerate(MOON_KEPLER_CASES, 1):
        sat = case["sat"]
        params = luna_kep_params(sat)
        period = case_period(case, GM_MOON)
        s.append(f"\n% case {k}: a = {case['a'] / 1e3} km, e = {case['e']} (Luna)\n")
        s.append(f"Report rfLuna{k} {params};\n")
        for dt in segment_dts(period):
            s.append(f"Propagate LunaTwoBodyProp({sat}) {{{sat}.ElapsedSecs = {dt!r}}};\n")
            s.append(f"Report rfLuna{k} {params};\n")

    return "".join(s)


# ----------------------------------------------------------------------------
# ReportFile parsing
#
# Each `Report rf <params...>;` call writes a header line (the parameter
# names, space-padded to ColumnWidth) followed by a single data line. With N
# `Report` calls referencing the same ReportFile, the file therefore consists
# of N repeated (header, data) line pairs.
# ----------------------------------------------------------------------------
def parse_report_blocks(path: Path) -> list[dict[str, float]]:
    lines = [line for line in path.read_text().splitlines() if line.strip()]
    if len(lines) % 2 != 0:
        raise ValueError(f"{path}: expected an even number of non-blank lines, got {len(lines)}")
    blocks = []
    for i in range(0, len(lines), 2):
        header = lines[i].split()
        values = [float(x) for x in lines[i + 1].split()]
        if len(header) != len(values):
            raise ValueError(f"{path}: header/value column count mismatch at line {i + 1}")
        blocks.append(dict(zip(header, values)))
    return blocks


def read_cart(blk: dict[str, float], sat: str, cs: str) -> list[float]:
    """Cartesian state [m, m/s] reported in coordinate system `cs`."""
    return [blk[f"{sat}.{cs}.{c}"] * 1e3 for c in ("X", "Y", "Z", "VX", "VY", "VZ")]


def extract_kep_state(blk: dict[str, float], sat: str, origin: str = "", cs: str = ""):
    """Extract (cart [m, m/s], coe [m, rad], E [rad], nu [rad]).

    `origin`/`cs` are the parameter-dependency qualifiers: empty for Earth
    spacecraft (plain `Sat.SMA`, `Sat.X`), `"Luna."`/`"LunaICRFCS."` for
    lunar spacecraft (`Sat.Luna.SMA`, `Sat.LunaICRFCS.X`, ...).
    """
    cart = [blk[f"{sat}.{cs}{c}"] * 1e3 for c in ("X", "Y", "Z", "VX", "VY", "VZ")]
    a = blk[f"{sat}.{origin}SMA"] * 1e3
    e = blk[f"{sat}.{origin}ECC"]
    coe = [
        a,
        e,
        math.radians(blk[f"{sat}.{cs}INC"]),
        math.radians(blk[f"{sat}.{cs}RAAN"]),
        math.radians(blk[f"{sat}.{cs}AOP"]),
        math.radians(blk[f"{sat}.{origin}MA"]),
    ]
    nu = math.radians(blk[f"{sat}.{origin}TA"])
    E = true_to_ecc_anomaly(nu, e)
    return cart, coe, E, nu


def kepler_case_json(case: dict, gm: float, body: str, blocks: list, origin: str, cs: str) -> dict:
    M0 = math.radians(case["m0"])
    E0 = mean_to_ecc_anomaly(M0, case["e"])
    nu0 = ecc_to_true_anomaly(E0, case["e"])
    period = case_period(case, gm)
    coe0 = [
        case["a"],
        case["e"],
        math.radians(case["i"]),
        math.radians(case["raan"]),
        math.radians(case["aop"]),
        M0,
    ]

    cart0, _, _, _ = extract_kep_state(blocks[0], case["sat"], origin, cs)
    out = {
        "gm": gm,
        "body": body,
        "epoch_utc": epoch_str(KEPLER_EPOCH),
        "coe0": coe0,
        "cart0": cart0,
        "ecc_anomaly0": E0,
        "true_anomaly0": nu0,
        "period": period,
        "propagated": [],
    }
    for frac, blk in zip(DT_FRACTIONS, blocks[1:]):
        cart, coe, E, nu = extract_kep_state(blk, case["sat"], origin, cs)
        out["propagated"].append(
            {"dt": period * frac, "cart": cart, "coe": coe, "ecc_anomaly": E, "true_anomaly": nu}
        )
    return out


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main() -> None:
    gmat_console = find_gmat_console()
    print(f"Using GmatConsole: {gmat_console}")

    rf_names = (
        ["rfTime", "rfFrames", "rfMoon"]
        + [f"rfKep{k}" for k in range(1, len(EARTH_KEPLER_CASES) + 1)]
        + [f"rfLuna{k}" for k in range(1, len(MOON_KEPLER_CASES) + 1)]
    )

    with tempfile.TemporaryDirectory(prefix="lupnt_gmat_ref_") as tmpdir:
        tmp = Path(tmpdir)
        rf_paths = {name: str(tmp / f"{name}.txt") for name in rf_names}
        script = build_script(str(COF_PATH.resolve()), rf_paths)

        # Documentation copy with machine-specific paths replaced by tokens.
        doc_script = script.replace(str(COF_PATH.resolve()), "<LUPNT_J2_EARTH_COF>")
        doc_script = doc_script.replace(DE440_PATH, "<DE440.bsp>")
        for name, p in rf_paths.items():
            doc_script = doc_script.replace(p, f"<{name.upper()}_PATH>")
        SCRIPT_DOC_PATH.write_text(doc_script)

        script_path = tmp / "lupnt_gmat_reference.script"
        script_path.write_text(script)

        # GMAT resets its working directory to <gmat_install>/bin, so all
        # paths above (and the -r argument) must be absolute.
        print(f"Running GMAT on {script_path} ...")
        result = subprocess.run(
            [str(gmat_console), "-r", str(script_path)],
            cwd=str(gmat_console.parent),
            capture_output=True,
            text=True,
        )
        print(result.stdout)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        if result.returncode != 0 or GMAT_SUCCESS_MARKER not in result.stdout:
            raise RuntimeError(f"GMAT run failed (exit code {result.returncode}); see output above")

        blocks = {name: parse_report_blocks(Path(rf_paths[name])) for name in rf_names}

    out: dict = {
        "_meta": {
            "description": (
                "Reference vectors generated by GMAT for cross-checking LuPNT's "
                "time-scale conversions, frame conversions (Earth and Moon), and "
                "two-body Keplerian + J2 orbit propagation. See "
                "cpp/test/gmat/README.md."
            ),
            "generator": "cpp/test/gmat/gen_gmat_reference.py",
            "gmat_version": "R2026a",
            "gmat_ephemeris": "DE440 (LuPNT's, via SPICE)",
            "constants": {
                "GM_EARTH": GM_EARTH,
                "GM_MOON": GM_MOON,
                "R_EARTH": R_EARTH,
                "J2_EARTH": J2_EARTH,
            },
            "notes": [
                "GMAT R2026a does not expose a UT1 or GPS time scale as a "
                "spacecraft parameter, so `time_scales` only has TAI/TT/TDB "
                "offsets from UTC (orekit_reference.json additionally has "
                "gps_minus_tai and ut1_minus_utc).",
                "`coe0`, `ecc_anomaly0`, `true_anomaly0`, and `period` in "
                "`kepler`/`j2_propagation` are computed directly in this "
                "script via a standard Kepler-equation solver, independent "
                "of GMAT (GMAT's Keplerian display state only accepts True "
                "Anomaly as the 6th element, so the script solves for the "
                "equivalent True Anomaly to set up each spacecraft). "
                "`cart0` and all `propagated` entries are genuine GMAT "
                "outputs.",
                "GMAT R2026a ships only up to DE424 (no DE440), so this script "
                "points GMAT at LuPNT's own de440.bsp via "
                "SolarSystem.EphemerisSource = 'SPICE'; the `moon_frames` "
                "section is therefore generated on DE440, the same ephemeris "
                "LuPNT and Orekit use, and the former DE421-vs-DE440 lunar "
                "offset is gone -- what remains is an ephemeris "
                "distribution/interpolation residual. `r_moon_fixed` is GMAT's "
                "Luna BodyFixed coordinate system.",
                "The J2 propagation in `j2_propagation` differs from "
                "LuPNT's JToCartTwoBodyDynamics by tens to a couple hundred "
                "meters after 0.25-1.5 orbits. Both evaluate the degree-2 "
                "zonal term in the Earth body-fixed frame (LuPNT's "
                "JToCartTwoBodyDynamics defaults to the body-fixed frame, "
                "matching GMAT's GravityField), so the residual is "
                "Earth-orientation (ITRF-vs-GMAT-EarthFixed realisation) and "
                "numerical-integrator differences, not an axis-convention "
                "difference. See cpp/test/gmat/README.md for details.",
            ],
        }
    }

    # --- Section A: time scales -------------------------------------------
    DAY = 86400.0
    time_scales = []
    for k, (ep, blk) in enumerate(zip(TIME_EPOCHS, blocks["rfTime"]), 1):
        utc = blk[f"TimeSat{k}.UTCModJulian"]
        tai = blk[f"TimeSat{k}.TAIModJulian"]
        tt = blk[f"TimeSat{k}.TTModJulian"]
        tdb = blk[f"TimeSat{k}.TDBModJulian"]
        time_scales.append(
            {
                "epoch_utc": epoch_str(ep),
                "tai_minus_utc": (tai - utc) * DAY,
                "tt_minus_tai": (tt - tai) * DAY,
                "tdb_minus_tai": (tdb - tai) * DAY,
            }
        )
    out["time_scales"] = time_scales

    # --- Section B: Earth frames --------------------------------------------
    frames = []
    blk_iter = iter(blocks["rfFrames"])
    for ei, ep in enumerate(FRAME_EPOCHS, 1):
        for si in range(1, len(FRAME_STATES) + 1):
            name = f"Sat_e{ei}_s{si}"
            blk = next(blk_iter)
            cart_icrf = read_cart(blk, name, "EarthICRFCS")
            cart_eme = read_cart(blk, name, "EarthMJ2000Eq")
            cart_itrf = read_cart(blk, name, "EarthFixedCS")
            frames.append(
                {
                    "epoch_utc": epoch_str(ep),
                    # GMAT's "ICRF" axes correspond to LuPNT/Orekit's GCRF.
                    "r_gcrf": cart_icrf[:3],
                    "v_gcrf": cart_icrf[3:],
                    # GMAT's "EarthMJ2000Eq" axes correspond to LuPNT's EME2000.
                    "r_eme2000": cart_eme[:3],
                    "v_eme2000": cart_eme[3:],
                    # GMAT's "BodyFixed" axes correspond to LuPNT's ITRF.
                    "r_itrf": cart_itrf[:3],
                    "v_itrf": cart_itrf[3:],
                }
            )
    out["frames"] = frames

    # --- Section C: Moon frames ---------------------------------------------
    moon_frames = []
    blk_iter = iter(blocks["rfMoon"])
    for ei, ep in enumerate(MOON_EPOCHS, 1):
        for si in range(1, len(MOON_STATES) + 1):
            name = f"MoonSat_e{ei}_s{si}"
            blk = next(blk_iter)
            cart_gcrf = read_cart(blk, name, "EarthICRFCS")
            cart_ci = read_cart(blk, name, "LunaICRFCS")
            cart_fixed = read_cart(blk, name, "LunaFixedCS")
            moon_frames.append(
                {
                    "epoch_utc": epoch_str(ep),
                    "r_gcrf": cart_gcrf[:3],
                    "v_gcrf": cart_gcrf[3:],
                    "r_moon_ci": cart_ci[:3],
                    "v_moon_ci": cart_ci[3:],
                    "r_moon_fixed": cart_fixed[:3],
                    "v_moon_fixed": cart_fixed[3:],
                }
            )
    out["moon_frames"] = moon_frames

    # --- Sections D/E: Keplerian (+J2) propagation ---------------------------
    kepler = []
    j2_propagation = []
    for k, case in enumerate(EARTH_KEPLER_CASES, 1):
        kep_blocks = blocks[f"rfKep{k}"]
        kepler.append(kepler_case_json(case, GM_EARTH, "EARTH", kep_blocks, "", ""))
        if case["j2"]:
            period = case_period(case, GM_EARTH)
            cart0_j2, _, _, _ = extract_kep_state(kep_blocks[0], case["j2sat"], "", "")
            case_j2 = {
                "gm": GM_EARTH,
                "j2": J2_EARTH,
                "r_earth": R_EARTH,
                "epoch_utc": epoch_str(KEPLER_EPOCH),
                "cart0": cart0_j2,
                "period": period,
                "propagated": [],
            }
            for frac, blk in zip(DT_FRACTIONS, kep_blocks[1:]):
                cart, _, _, _ = extract_kep_state(blk, case["j2sat"], "", "")
                case_j2["propagated"].append({"dt": period * frac, "cart": cart})
            j2_propagation.append(case_j2)

    for k, case in enumerate(MOON_KEPLER_CASES, 1):
        kepler.append(
            kepler_case_json(case, GM_MOON, "MOON", blocks[f"rfLuna{k}"], "Luna.", "LunaICRFCS.")
        )

    out["kepler"] = kepler
    out["j2_propagation"] = j2_propagation

    OUTPUT_PATH.write_text(json.dumps(out, indent=2) + "\n")
    print(f"Wrote {OUTPUT_PATH}")
    print(f"Wrote {SCRIPT_DOC_PATH} (documentation copy)")


if __name__ == "__main__":
    main()

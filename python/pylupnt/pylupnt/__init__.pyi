from __future__ import annotations
from matplotlib.gridspec import GridSpec
from matplotlib import pyplot as plt
import numpy as np
import os as os
import pandas as pd
import pickle as pickle
from pylupnt._pylupnt import AzElRange2Cart
from pylupnt._pylupnt import Body
from pylupnt._pylupnt import Cart2AzElRange
from pylupnt._pylupnt import Cart2EastNorthUp
from pylupnt._pylupnt import Cart2LatLonAlt
from pylupnt._pylupnt import CartesianOrbitState
from pylupnt._pylupnt import CartesianTwoBodyDynamics
from pylupnt._pylupnt import ClassicalOE
from pylupnt._pylupnt import Decibel2Decimal
from pylupnt._pylupnt import Decimal2Decibel
from pylupnt._pylupnt import DegMinSec2Degrees
from pylupnt._pylupnt import Degrees2DegMinSec
from pylupnt._pylupnt import EastNorthUp2Cart
from pylupnt._pylupnt import EquinoctialOE
from pylupnt._pylupnt import Frame
from pylupnt._pylupnt import KeplerianDynamics
from pylupnt._pylupnt import LatLonAlt2Cart
from pylupnt._pylupnt import NBodyDynamics
from pylupnt._pylupnt import NaifId
from pylupnt._pylupnt import NumericalOrbitDynamics
from pylupnt._pylupnt import OrbitState
from pylupnt._pylupnt import OrbitStateRepres
from pylupnt._pylupnt import QuasiNonsingOE
from pylupnt._pylupnt import QuasiNonsingROE
from pylupnt._pylupnt import SingularROE
from pylupnt._pylupnt import Wrap2Pi
from pylupnt._pylupnt import Wrap2TwoPi
from pylupnt._pylupnt import cartesian_to_classical
from pylupnt._pylupnt import classical_to_cartesian
from pylupnt._pylupnt import classical_to_delaunay
from pylupnt._pylupnt import classical_to_equinoctial
from pylupnt._pylupnt import classical_to_quasi_nonsingular
from pylupnt._pylupnt import compute_occultation
from pylupnt._pylupnt import convert_frame
from pylupnt._pylupnt import convert_orbit_state
from pylupnt._pylupnt import delaunay_to_classical
from pylupnt._pylupnt import eccentric_to_mean_anomaly
from pylupnt._pylupnt import eccentric_to_true_anomaly
from pylupnt._pylupnt import equinoctial_to_classical
from pylupnt._pylupnt import mean_to_eccentric_anomaly
from pylupnt._pylupnt import mean_to_true_anomaly
from pylupnt._pylupnt import quasi_nonsingular_to_classical
from pylupnt._pylupnt import relative_quasi_nonsingular_to_classical
from pylupnt._pylupnt import true_to_eccentric_anomaly
from pylupnt._pylupnt import true_to_mean_anomaly
from pylupnt.math_utils import arr_to_mat_idx
from pylupnt.math_utils import cross_norm
from pylupnt.math_utils import i_to_arr_idxs
from pylupnt.math_utils import mat_to_arr_idx
from pylupnt.math_utils import wrap2Pi
from pylupnt.utils import dump_pickle
from pylupnt.utils import find_file
from pylupnt.utils import format_element
from pylupnt.utils import get_basepath
from pylupnt.utils import load_data
from pylupnt.utils import load_pickle
from pylupnt.utils import plot_RTN
from pylupnt.utils import print_aligned
from pylupnt.utils import set_axes_equal
from pylupnt.utils import timed
from pylupnt.utils import timer_func
import time as time
from . import _pylupnt
from . import math_utils
from . import plot
from . import render
from . import scenarios
from . import utils

__all__ = [
    "A1_TAI_OFFSET",
    "ARCSEC_DEG",
    "ARCSEC_RAD",
    "AzElRange2Cart",
    "Body",
    "CARTESIAN",
    "CLASSICAL_OE",
    "Cart2AzElRange",
    "Cart2EastNorthUp",
    "Cart2LatLonAlt",
    "CartesianOrbitState",
    "CartesianTwoBodyDynamics",
    "ClassicalOE",
    "DAYS_CENTURY",
    "DAYS_SEC",
    "DAYS_WEEK",
    "DAYS_YEAR",
    "DEG",
    "DEG_ARCSEC",
    "DEIMOS",
    "DELAUNAY_OE",
    "Decibel2Decimal",
    "Decimal2Decibel",
    "DegMinSec2Degrees",
    "Degrees2DegMinSec",
    "E",
    "EARTH",
    "EARTH_MOON_BARYCENTER",
    "ECEF",
    "ECI",
    "EMB",
    "EME",
    "EMR",
    "EPS",
    "EQUINOTICAL_OE",
    "EastNorthUp2Cart",
    "EquinoctialOE",
    "FOOT_M",
    "Frame",
    "GCRF",
    "GM_EARTH",
    "GM_MERCURY",
    "GM_SUN",
    "GM_VENUS",
    "GPS",
    "GSE",
    "GridSpec",
    "HOURS_DAY",
    "ICRF",
    "INCH_M",
    "ITRF",
    "JD_J2000",
    "JD_JAN_5_1941",
    "JD_MJD_OFFSET",
    "JD_NOV_17_1858",
    "JD_T0",
    "JD_TDB",
    "JD_TT",
    "JULIAN_DATE_OF_010541",
    "JUPITER",
    "JUPITER_BARYCENTER",
    "KM_M",
    "KeplerianDynamics",
    "LBM_TO_KG",
    "LUPNT_DATA_PATH",
    "L_B",
    "L_G",
    "LatLonAlt2Cart",
    "MARS",
    "MARS_BARYCENTER",
    "MARS_FIXED",
    "MERCURY",
    "MERCURY_BARYCENTER",
    "MILE_M",
    "MINS_DAY",
    "MINS_HOUR",
    "MJD_J2000",
    "MOD",
    "MOON",
    "MOON_CI",
    "MOON_ME",
    "MOON_OP",
    "MOON_PA",
    "M_KM",
    "NBodyDynamics",
    "NEPTUNE_BARYCENTER",
    "NUM_SECS",
    "NaifId",
    "NumericalOrbitDynamics",
    "OrbitState",
    "OrbitStateRepres",
    "PHOBOS",
    "PI",
    "PI_OVER_TWO",
    "PLUTO_BARYCENTER",
    "QUASINONSINGULAR_ROE",
    "QUASI_NONSINGULAR_OE",
    "QuasiNonsingOE",
    "QuasiNonsingROE",
    "RAD",
    "RAD_ARCSEC",
    "R_EARTH",
    "R_MOON",
    "SATURN_BARYCENTER",
    "SECS_DAY",
    "SECS_HOUR",
    "SECS_MINUTE",
    "SER",
    "SINGULAR_ROE",
    "SLUG_TO_KG",
    "SOLAR_SYSTEM_BARYCENTER",
    "SSB",
    "SUN",
    "SingularROE",
    "TAI",
    "TCB",
    "TCG",
    "TDB",
    "TIME_OF_J2000",
    "TOD",
    "TT",
    "TT_TAI_OFFSET",
    "TWO_PI",
    "URANUS_BARYCENTER",
    "UT1",
    "UTC",
    "VENUS",
    "VENUS_BARYCENTER",
    "VENUS_FIXED",
    "Wrap2Pi",
    "Wrap2TwoPi",
    "arr_to_mat_idx",
    "cartesian_to_classical",
    "classical_to_cartesian",
    "classical_to_delaunay",
    "classical_to_equinoctial",
    "classical_to_quasi_nonsingular",
    "compute_occultation",
    "convert_frame",
    "convert_orbit_state",
    "cross_norm",
    "delaunay_to_classical",
    "dump_pickle",
    "eccentric_to_mean_anomaly",
    "eccentric_to_true_anomaly",
    "equinoctial_to_classical",
    "find_file",
    "format_element",
    "get_basepath",
    "i_to_arr_idxs",
    "load_data",
    "load_pickle",
    "mat_to_arr_idx",
    "math_utils",
    "mean_to_eccentric_anomaly",
    "mean_to_true_anomaly",
    "np",
    "os",
    "pd",
    "pickle",
    "plot",
    "plot_RTN",
    "plt",
    "print_aligned",
    "quasi_nonsingular_to_classical",
    "relative_quasi_nonsingular_to_classical",
    "render",
    "scenarios",
    "set_axes_equal",
    "time",
    "timed",
    "timer_func",
    "true_to_eccentric_anomaly",
    "true_to_mean_anomaly",
    "utils",
    "wrap2Pi",
]
A1_TAI_OFFSET: float = 0.0343817
ARCSEC_DEG: float = 3600.0
ARCSEC_RAD: float = 206264.80624709636
CARTESIAN: _pylupnt.OrbitStateRepres  # value = <OrbitStateRepres.CARTESIAN: 0>
CLASSICAL_OE: _pylupnt.OrbitStateRepres  # value = <OrbitStateRepres.CLASSICAL_OE: 1>
DAYS_CENTURY: float = 36525.0
DAYS_SEC: float = 1.1574074074074073e-05
DAYS_WEEK: float = 7.0
DAYS_YEAR: float = 365.25
DEG: float = 57.29577951308232
DEG_ARCSEC: float = 0.0002777777777777778
DEIMOS: _pylupnt.NaifId  # value = <NaifId.DEIMOS: 402>
DELAUNAY_OE: _pylupnt.OrbitStateRepres  # value = <OrbitStateRepres.DELAUNAY_OE: 6>
E: float = 2.718281828459045
EARTH: _pylupnt.NaifId  # value = <NaifId.EARTH: 399>
EARTH_MOON_BARYCENTER: _pylupnt.NaifId  # value = <NaifId.EMB: 3>
ECEF: _pylupnt.Frame  # value = <Frame.ITRF: 0>
ECI: _pylupnt.Frame  # value = <Frame.ECI: 2>
EMB: _pylupnt.NaifId  # value = <NaifId.EMB: 3>
EME: _pylupnt.Frame  # value = <Frame.ECI: 2>
EMR: _pylupnt.Frame  # value = <Frame.EMR: 8>
EPS: float = 1e-16
EQUINOTICAL_OE: (
    _pylupnt.OrbitStateRepres
)  # value = <OrbitStateRepres.EQUINOTICAL_OE: 5>
FOOT_M: float = 0.3048
GCRF: _pylupnt.Frame  # value = <Frame.GCRF: 1>
GM_EARTH: float = 398600.435507
GM_MERCURY: float = 22031.868551
GM_SUN: float = 132712440041.27942
GM_VENUS: float = 324858.592
GPS: str = "GPS"
GSE: _pylupnt.Frame  # value = <Frame.GSE: 5>
HOURS_DAY: float = 24.0
ICRF: _pylupnt.Frame  # value = <Frame.ICRF: 3>
INCH_M: float = 0.0254
ITRF: _pylupnt.Frame  # value = <Frame.ITRF: 0>
JD_J2000: float = 2451545.0
JD_JAN_5_1941: float = 2430000.0
JD_MJD_OFFSET: float = 2400000.5
JD_NOV_17_1858: float = 2400000.5
JD_T0: float = 2443144.5003725
JD_TDB: str = "JDTDB"
JD_TT: str = "JDTDT"
JULIAN_DATE_OF_010541: int = 2430000
JUPITER: _pylupnt.NaifId  # value = <NaifId.JUPITER: 599>
JUPITER_BARYCENTER: _pylupnt.NaifId  # value = <NaifId.JUPITER_BARYCENTER: 5>
KM_M: float = 0.001
LBM_TO_KG: float = 0.45359237
LUPNT_DATA_PATH: str = (
    "/Users/guillemcv/Library/CloudStorage/GoogleDrive-guillemc@stanford.edu/My Drive/LunaNet Drive/LuPNT_data"
)
L_B: float = 1.550505e-08
L_G: float = 6.969290134e-10
MARS: _pylupnt.NaifId  # value = <NaifId.MARS: 499>
MARS_BARYCENTER: _pylupnt.NaifId  # value = <NaifId.MARS_BARYCENTER: 4>
MARS_FIXED: _pylupnt.Frame  # value = <Frame.MARS_FIXED: 13>
MERCURY: _pylupnt.NaifId  # value = <NaifId.MERCURY: 199>
MERCURY_BARYCENTER: _pylupnt.NaifId  # value = <NaifId.MERCURY_BARYCENTER: 1>
MILE_M: float = 1609.344
MINS_DAY: float = 1440.0
MINS_HOUR: float = 60.0
MJD_J2000: float = 51544.5
MOD: _pylupnt.Frame  # value = <Frame.MOD: 6>
MOON: _pylupnt.NaifId  # value = <NaifId.MOON: 301>
MOON_CI: _pylupnt.Frame  # value = <Frame.MOON_CI: 9>
MOON_ME: _pylupnt.Frame  # value = <Frame.MOON_ME: 11>
MOON_OP: _pylupnt.Frame  # value = <Frame.MOON_OP: 12>
MOON_PA: _pylupnt.Frame  # value = <Frame.MOON_PA: 10>
M_KM: float = 1000.0
NEPTUNE_BARYCENTER: _pylupnt.NaifId  # value = <NaifId.NEPTUNE_BARYCENTER: 8>
NUM_SECS: float = 86400.0
PHOBOS: _pylupnt.NaifId  # value = <NaifId.PHOBOS: 401>
PI: float = 3.141592653589793
PI_OVER_TWO: float = 1.5707963267948966
PLUTO_BARYCENTER: _pylupnt.NaifId  # value = <NaifId.PLUTO_BARYCENTER: 9>
QUASINONSINGULAR_ROE: (
    _pylupnt.OrbitStateRepres
)  # value = <OrbitStateRepres.QUASINONSINGULAR_ROE: 9>
QUASI_NONSINGULAR_OE: (
    _pylupnt.OrbitStateRepres
)  # value = <OrbitStateRepres.QUASI_NONSINGULAR_OE: 2>
RAD: float = 0.017453292519943295
RAD_ARCSEC: float = 4.84813681109536e-06
R_EARTH: float = 6378.137
R_MOON: float = 1737.4
SATURN_BARYCENTER: _pylupnt.NaifId  # value = <NaifId.SATURN_BARYCENTER: 6>
SECS_DAY: float = 86400.0
SECS_HOUR: float = 3600.0
SECS_MINUTE: float = 60.0
SER: _pylupnt.Frame  # value = <Frame.SER: 4>
SINGULAR_ROE: _pylupnt.OrbitStateRepres  # value = <OrbitStateRepres.SINGULAR_ROE: 3>
SLUG_TO_KG: float = 14.59390294
SOLAR_SYSTEM_BARYCENTER: _pylupnt.NaifId  # value = <NaifId.SOLAR_SYSTEM_BARYCENTER: 0>
SSB: _pylupnt.NaifId  # value = <NaifId.SOLAR_SYSTEM_BARYCENTER: 0>
SUN: _pylupnt.NaifId  # value = <NaifId.SUN: 10>
TAI: str = "TAI"
TCB: str = "TCB"
TCG: str = "TCG"
TDB: str = "TDB"
TIME_OF_J2000: float = 883655990.85
TOD: _pylupnt.Frame  # value = <Frame.TOD: 7>
TT: str = "TT"
TT_TAI_OFFSET: float = 32.184
TWO_PI: float = 6.283185307179586
URANUS_BARYCENTER: _pylupnt.NaifId  # value = <NaifId.URANUS_BARYCENTER: 7>
UT1: str = "UT1"
UTC: str = "UTC"
VENUS: _pylupnt.NaifId  # value = <NaifId.VENUS: 299>
VENUS_BARYCENTER: _pylupnt.NaifId  # value = <NaifId.VENUS_BARYCENTER: 2>
VENUS_FIXED: _pylupnt.Frame  # value = <Frame.VENUS_FIXED: 14>

from __future__ import annotations
from matplotlib.gridspec import GridSpec
from matplotlib import pyplot as plt
import numpy as np
import os as os
import pandas as pd
import pickle as pickle
import pkg_resources as pkg_resources
from pylupnt.math_utils import arr_to_mat_idx
from pylupnt.math_utils import cross_norm
from pylupnt.math_utils import i_to_arr_idxs
from pylupnt.math_utils import mat_to_arr_idx
from pylupnt.math_utils import wrapToPi
from pylupnt.pylupnt_pybind import Body
from pylupnt.pylupnt_pybind import CartesianOrbitState
from pylupnt.pylupnt_pybind import CartesianTwoBodyDynamics
from pylupnt.pylupnt_pybind import ClassicalOE
from pylupnt.pylupnt_pybind import CoordConverter
from pylupnt.pylupnt_pybind import CoordSystem
from pylupnt.pylupnt_pybind import EquinoctialOE
from pylupnt.pylupnt_pybind import KeplerianDynamics
from pylupnt.pylupnt_pybind import NBodyDynamics
from pylupnt.pylupnt_pybind import NaifId
from pylupnt.pylupnt_pybind import NumericalOrbitDynamics
from pylupnt.pylupnt_pybind import OrbitState
from pylupnt.pylupnt_pybind import OrbitStateRepres
from pylupnt.pylupnt_pybind import QuasiNonsingularOE
from pylupnt.pylupnt_pybind import QuasiNonsingularROE
from pylupnt.pylupnt_pybind import SingularROE
from pylupnt.pylupnt_pybind import SpiceInterface
from pylupnt.pylupnt_pybind import azimuth_elevation_range_to_cartesian
from pylupnt.pylupnt_pybind import cartesian_to_azimuth_elevation_range
from pylupnt.pylupnt_pybind import cartesian_to_classical
from pylupnt.pylupnt_pybind import cartesian_to_east_north_up
from pylupnt.pylupnt_pybind import cartesian_to_geographical
from pylupnt.pylupnt_pybind import cartesian_to_spherical
from pylupnt.pylupnt_pybind import classical_to_cartesian
from pylupnt.pylupnt_pybind import classical_to_delaunay
from pylupnt.pylupnt_pybind import classical_to_equinoctial
from pylupnt.pylupnt_pybind import classical_to_quasi_nonsingular
from pylupnt.pylupnt_pybind import compute_occultation
from pylupnt.pylupnt_pybind import convert_orbit_state
from pylupnt.pylupnt_pybind import dB2decimal
from pylupnt.pylupnt_pybind import decimal2dB
from pylupnt.pylupnt_pybind import degrees2dms
from pylupnt.pylupnt_pybind import delaunay_to_classical
from pylupnt.pylupnt_pybind import dms2degrees
from pylupnt.pylupnt_pybind import east_north_up_to_cartesian
from pylupnt.pylupnt_pybind import eccentric_to_mean_anomaly
from pylupnt.pylupnt_pybind import eccentric_to_true_anomaly
from pylupnt.pylupnt_pybind import equinoctial_to_classical
from pylupnt.pylupnt_pybind import geographical_to_cartesian
from pylupnt.pylupnt_pybind import mean_to_eccentric_anomaly
from pylupnt.pylupnt_pybind import mean_to_true_anomaly
from pylupnt.pylupnt_pybind import quasi_nonsingular_to_classical
from pylupnt.pylupnt_pybind import relative_quasi_nonsingular_to_classical
from pylupnt.pylupnt_pybind import spherical_to_cartesian
from pylupnt.pylupnt_pybind import true_to_eccentric_anomaly
from pylupnt.pylupnt_pybind import true_to_mean_anomaly
from pylupnt.pylupnt_pybind import wrapTo2Pi
from pylupnt.utils import dump_pickle
from pylupnt.utils import format_element
from pylupnt.utils import get_basepath
from pylupnt.utils import load_data
from pylupnt.utils import load_pickle
from pylupnt.utils import plot_RTN
from pylupnt.utils import print_aligned
from pylupnt.utils import set_axes_equal
from pylupnt.utils import timed
from pylupnt.utils import timer_func
from time import time
from . import math_utils
from . import plotting
from . import pylupnt_pybind
from . import utils
__all__ = ['A1MJD_OF_J2000', 'A1_TAI_OFFSET', 'AU', 'Body', 'C', 'C22_MOON', 'CARTESIAN', 'CLASSICAL_OE', 'CartesianOrbitState', 'CartesianTwoBodyDynamics', 'ClassicalOE', 'CoordConverter', 'CoordSystem', 'DAYS_PER_JULIAN_CENTURY', 'DAYS_PER_SEC', 'DAYS_PER_YEAR', 'DEG_PER_RAD', 'DEIMOS', 'DELAUNAY_OE', 'E', 'EARTH', 'EARTH_BARYCENTER', 'EARTH_MOON_BARYCENTER', 'ECEF', 'ECI', 'EME', 'EMR', 'EQUINOTICAL_OE', 'EquinoctialOE', 'GCRF', 'GSE', 'GridSpec', 'ICRF', 'ITRF', 'J2_EARTH', 'J2_MOON', 'JD_JAN_5_1941', 'JD_MJD_OFFSET', 'JD_NOV_17_1858', 'JD_OF_J2000', 'JUPITER', 'JUPITER_BARYCENTER', 'KeplerianDynamics', 'LUPNT_DATA_PATH', 'MARS', 'MARS_BARYCENTER', 'ME', 'MERCURY', 'MERCURY_BARYCENTER', 'MI', 'MJD_OF_J2000', 'MOD', 'MOON', 'MU_EARTH', 'MU_MOON', 'NBodyDynamics', 'NEPTUNE_BARYCENTER', 'NaifId', 'NumericalOrbitDynamics', 'OMEGA_E_M', 'OP', 'OrbitState', 'OrbitStateRepres', 'PA', 'PHOBOS', 'PI', 'PI_OVER_TWO', 'PLUTO_BARYCENTER', 'P_SUN', 'QUASINONSINGULAR_ROE', 'QUASI_NONSINGULAR_OE', 'QuasiNonsingularOE', 'QuasiNonsingularROE', 'RAD_PER_DEG', 'R_EARTH', 'R_MOON', 'SATURN_BARYCENTER', 'SECS_PER_DAY', 'SECS_PER_HOUR', 'SECS_PER_MINUTE', 'SER', 'SINGULAR_ROE', 'SOLAR_SYSTEM_BARYCENTER', 'SUN', 'S_AU', 'SingularROE', 'SpiceInterface', 'TIME_OF_J2000', 'TOD', 'TT_TAI_OFFSET', 'TWO_PI', 'URANUS_BARYCENTER', 'VENUS', 'VENUS_BARYCENTER', 'arr_to_mat_idx', 'azimuth_elevation_range_to_cartesian', 'cartesian_to_azimuth_elevation_range', 'cartesian_to_classical', 'cartesian_to_east_north_up', 'cartesian_to_geographical', 'cartesian_to_spherical', 'classical_to_cartesian', 'classical_to_delaunay', 'classical_to_equinoctial', 'classical_to_quasi_nonsingular', 'compute_occultation', 'convert_orbit_state', 'cross_norm', 'dB2decimal', 'd_E_EMB', 'd_E_M', 'd_M_EMB', 'decimal2dB', 'degrees2dms', 'delaunay_to_classical', 'dms2degrees', 'dump_pickle', 'east_north_up_to_cartesian', 'eccentric_to_mean_anomaly', 'eccentric_to_true_anomaly', 'equinoctial_to_classical', 'format_element', 'geographical_to_cartesian', 'get_basepath', 'i_to_arr_idxs', 'load_data', 'load_pickle', 'mat_to_arr_idx', 'math_utils', 'mean_to_eccentric_anomaly', 'mean_to_true_anomaly', 'np', 'os', 'pd', 'pickle', 'pkg_resources', 'plot_RTN', 'plotting', 'plt', 'print_aligned', 'pylupnt_pybind', 'quasi_nonsingular_to_classical', 'relative_quasi_nonsingular_to_classical', 'set_axes_equal', 'spherical_to_cartesian', 'time', 'timed', 'timer_func', 'true_to_eccentric_anomaly', 'true_to_mean_anomaly', 'utils', 'wrapTo2Pi', 'wrapToPi']
A1MJD_OF_J2000: float = 21545.0
A1_TAI_OFFSET: float = 0.0343817
AU: float = 149597970.0
C: float = 299792.458
C22_MOON: float = 3.470983013194e-05
CARTESIAN: pylupnt_pybind.OrbitStateRepres  # value = <OrbitStateRepres.CARTESIAN: 0>
CLASSICAL_OE: pylupnt_pybind.OrbitStateRepres  # value = <OrbitStateRepres.CLASSICAL_OE: 1>
DAYS_PER_JULIAN_CENTURY: float = 36525.0
DAYS_PER_SEC: float = 1.1574074074074073e-05
DAYS_PER_YEAR: float = 365.25
DEG_PER_RAD: float = 57.29577951308232
DEIMOS: pylupnt_pybind.NaifId  # value = <NaifId.DEIMOS: 402>
DELAUNAY_OE: pylupnt_pybind.OrbitStateRepres  # value = <OrbitStateRepres.DELAUNAY_OE: 6>
E: float = 2.718281828459045
EARTH: pylupnt_pybind.NaifId  # value = <NaifId.EARTH: 399>
EARTH_BARYCENTER: pylupnt_pybind.NaifId  # value = <NaifId.EARTH_BARYCENTER: 3>
EARTH_MOON_BARYCENTER: pylupnt_pybind.NaifId  # value = <NaifId.EARTH_BARYCENTER: 3>
ECEF: pylupnt_pybind.CoordSystem  # value = <CoordSystem.ITRF: 0>
ECI: pylupnt_pybind.CoordSystem  # value = <CoordSystem.GCRF: 1>
EME: pylupnt_pybind.CoordSystem  # value = <CoordSystem.EME: 5>
EMR: pylupnt_pybind.CoordSystem  # value = <CoordSystem.EMR: 8>
EQUINOTICAL_OE: pylupnt_pybind.OrbitStateRepres  # value = <OrbitStateRepres.EQUINOTICAL_OE: 5>
GCRF: pylupnt_pybind.CoordSystem  # value = <CoordSystem.GCRF: 1>
GSE: pylupnt_pybind.CoordSystem  # value = <CoordSystem.GSE: 4>
ICRF: pylupnt_pybind.CoordSystem  # value = <CoordSystem.ICRF: 2>
ITRF: pylupnt_pybind.CoordSystem  # value = <CoordSystem.ITRF: 0>
J2_EARTH: float = 0.00108262668
J2_MOON: float = 9.09427845027e-05
JD_JAN_5_1941: float = 2430000.0
JD_MJD_OFFSET: float = 2400000.5
JD_NOV_17_1858: float = 2400000.5
JD_OF_J2000: float = 2451545.0
JUPITER: pylupnt_pybind.NaifId  # value = <NaifId.JUPITER: 599>
JUPITER_BARYCENTER: pylupnt_pybind.NaifId  # value = <NaifId.JUPITER_BARYCENTER: 5>
LUPNT_DATA_PATH: str = '/Users/guillemcv/Development/NavLab/LuPNT/data'
MARS: pylupnt_pybind.NaifId  # value = <NaifId.MARS: 499>
MARS_BARYCENTER: pylupnt_pybind.NaifId  # value = <NaifId.MARS_BARYCENTER: 4>
ME: pylupnt_pybind.CoordSystem  # value = <CoordSystem.ME: 11>
MERCURY: pylupnt_pybind.NaifId  # value = <NaifId.MERCURY: 199>
MERCURY_BARYCENTER: pylupnt_pybind.NaifId  # value = <NaifId.MERCURY_BARYCENTER: 1>
MI: pylupnt_pybind.CoordSystem  # value = <CoordSystem.MI: 9>
MJD_OF_J2000: float = 21545.0
MOD: pylupnt_pybind.CoordSystem  # value = <CoordSystem.MOD: 6>
MOON: pylupnt_pybind.NaifId  # value = <NaifId.MOON: 301>
MU_EARTH: float = 398600.4418
MU_MOON: float = 4902.800066
NEPTUNE_BARYCENTER: pylupnt_pybind.NaifId  # value = <NaifId.NEPTUNE_BARYCENTER: 8>
OMEGA_E_M: float = 2.6617e-06
OP: pylupnt_pybind.CoordSystem  # value = <CoordSystem.OP: 12>
PA: pylupnt_pybind.CoordSystem  # value = <CoordSystem.PA: 10>
PHOBOS: pylupnt_pybind.NaifId  # value = <NaifId.PHOBOS: 401>
PI: float = 3.141592653589793
PI_OVER_TWO: float = 1.5707963267948966
PLUTO_BARYCENTER: pylupnt_pybind.NaifId  # value = <NaifId.PLUTO_BARYCENTER: 9>
P_SUN: float = 4.529800412790905e-09
QUASINONSINGULAR_ROE: pylupnt_pybind.OrbitStateRepres  # value = <OrbitStateRepres.QUASINONSINGULAR_ROE: 9>
QUASI_NONSINGULAR_OE: pylupnt_pybind.OrbitStateRepres  # value = <OrbitStateRepres.QUASI_NONSINGULAR_OE: 2>
RAD_PER_DEG: float = 0.017453292519943295
R_EARTH: float = 6378.137
R_MOON: float = 1737.4
SATURN_BARYCENTER: pylupnt_pybind.NaifId  # value = <NaifId.SATURN_BARYCENTER: 6>
SECS_PER_DAY: float = 86400.0
SECS_PER_HOUR: float = 3600.0
SECS_PER_MINUTE: float = 60.0
SER: pylupnt_pybind.CoordSystem  # value = <CoordSystem.SER: 3>
SINGULAR_ROE: pylupnt_pybind.OrbitStateRepres  # value = <OrbitStateRepres.SINGULAR_ROE: 3>
SOLAR_SYSTEM_BARYCENTER: pylupnt_pybind.NaifId  # value = <NaifId.SOLAR_SYSTEM_BARYCENTER: 0>
SUN: pylupnt_pybind.NaifId  # value = <NaifId.SUN: 10>
S_AU: float = 0.0013579999999999998
TIME_OF_J2000: float = 883655990.85
TOD: pylupnt_pybind.CoordSystem  # value = <CoordSystem.TOD: 7>
TT_TAI_OFFSET: float = 32.184
TWO_PI: float = 6.283185307179586
URANUS_BARYCENTER: pylupnt_pybind.NaifId  # value = <NaifId.URANUS_BARYCENTER: 7>
VENUS: pylupnt_pybind.NaifId  # value = <NaifId.VENUS: 299>
VENUS_BARYCENTER: pylupnt_pybind.NaifId  # value = <NaifId.VENUS_BARYCENTER: 2>
__version__: str = '@PROJECT_VERSION@'
__warningregistry__: dict = {'version': 121}
d_E_EMB: float = 4671.0
d_E_M: float = 384400.0
d_M_EMB: float = 379729.0

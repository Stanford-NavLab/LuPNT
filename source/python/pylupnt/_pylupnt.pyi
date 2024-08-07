from __future__ import annotations
import numpy
import pylupnt
import typing
__all__ = ['A1_TAI_OFFSET', 'ARCSEC_DEG', 'ARCSEC_RAD', 'AzElRange2Cart', 'Body', 'CARTESIAN', 'CLASSICAL_OE', 'Cart2AzElRange', 'Cart2EastNorthUp', 'Cart2LatLonAlt', 'CartesianOrbitState', 'CartesianTwoBodyDynamics', 'ClassicalOE', 'DAYS_CENTURY', 'DAYS_SEC', 'DAYS_WEEK', 'DAYS_YEAR', 'DEG', 'DEG_ARCSEC', 'DEIMOS', 'DELAUNAY_OE', 'Decibel2Decimal', 'Decimal2Decibel', 'DegMinSec2Degrees', 'Degrees2DegMinSec', 'E', 'EARTH', 'EARTH_MOON_BARYCENTER', 'ECEF', 'ECI', 'EMB', 'EME', 'EMR', 'EPS', 'EQUINOTICAL_OE', 'EastNorthUp2Cart', 'EquinoctialOE', 'FOOT_M', 'Frame', 'GCRF', 'GM_EARTH', 'GM_MERCURY', 'GM_MOON', 'GM_SUN', 'GM_VENUS', 'GPS', 'GSE', 'HOURS_DAY', 'ICRF', 'INCH_M', 'ITRF', 'JD_J2000', 'JD_JAN_5_1941', 'JD_MJD_OFFSET', 'JD_NOV_17_1858', 'JD_T0', 'JD_TDB', 'JD_TT', 'JULIAN_DATE_OF_010541', 'JUPITER', 'JUPITER_BARYCENTER', 'KM_M', 'KeplerianDynamics', 'LBM_TO_KG', 'L_B', 'L_G', 'LatLonAlt2Cart', 'MARS', 'MARS_BARYCENTER', 'MARS_FIXED', 'MERCURY', 'MERCURY_BARYCENTER', 'MILE_M', 'MINS_DAY', 'MINS_HOUR', 'MJD_J2000', 'MOD', 'MOON', 'MOON_CI', 'MOON_ME', 'MOON_OP', 'MOON_PA', 'M_KM', 'NBodyDynamics', 'NEPTUNE_BARYCENTER', 'NUM_SECS', 'NaifId', 'NumericalOrbitDynamics', 'OrbitState', 'OrbitStateRepres', 'PHOBOS', 'PI', 'PI_OVER_TWO', 'PLUTO_BARYCENTER', 'QUASINONSINGULAR_ROE', 'QUASI_NONSINGULAR_OE', 'QuasiNonsingOE', 'QuasiNonsingROE', 'RAD', 'RAD_ARCSEC', 'R_EARTH', 'R_MOON', 'SATURN_BARYCENTER', 'SECS_DAY', 'SECS_HOUR', 'SECS_MINUTE', 'SER', 'SINGULAR_ROE', 'SLUG_TO_KG', 'SOLAR_SYSTEM_BARYCENTER', 'SSB', 'SUN', 'SingularROE', 'SpiceInterface', 'TAI', 'TCB', 'TCG', 'TDB', 'TIME_OF_J2000', 'TOD', 'TT', 'TT_TAI_OFFSET', 'TWO_PI', 'URANUS_BARYCENTER', 'UT1', 'UTC', 'VENUS', 'VENUS_BARYCENTER', 'VENUS_FIXED', 'Wrap2Pi', 'Wrap2TwoPi', 'cartesian_to_classical', 'classical_to_cartesian', 'classical_to_delaunay', 'classical_to_equinoctial', 'classical_to_quasi_nonsingular', 'compute_occultation', 'convert_frame', 'convert_orbit_state', 'delaunay_to_classical', 'eccentric_to_mean_anomaly', 'eccentric_to_true_anomaly', 'equinoctial_to_classical', 'mean_to_eccentric_anomaly', 'mean_to_true_anomaly', 'quasi_nonsingular_to_classical', 'relative_quasi_nonsingular_to_classical', 'true_to_eccentric_anomaly', 'true_to_mean_anomaly']
class Body:
    @staticmethod
    def Earth(n_max: int = 0, m_max: int = 0, gravity_file: str = 'grgm900c.cof') -> Body:
        ...
    @staticmethod
    def Mars(n_max: int = 0, m_max: int = 0, gravity_file: str = 'GMM1.cof') -> Body:
        ...
    @staticmethod
    def Moon(n_max: int = 0, m_max: int = 0, gravity_file: str = 'EGM96.cof') -> Body:
        ...
    @staticmethod
    def Sun() -> Body:
        ...
    @staticmethod
    def Venus(n_max: int = 0, m_max: int = 0, gravity_file: str = 'MGN75HSAAP.cof') -> Body:
        ...
    def __init__(self) -> None:
        ...
    def __repr__(self) -> str:
        ...
class CartesianOrbitState(OrbitState):
    r: numpy.ndarray[numpy.float64[3, 1]]
    v: numpy.ndarray[numpy.float64[3, 1]]
    def __init__(self, rv: numpy.ndarray[numpy.float64[6, 1]], frame: Frame = pylupnt.Frame.GCRF) -> None:
        ...
    def __repr__(self) -> str:
        ...
class CartesianTwoBodyDynamics(NumericalOrbitDynamics):
    def __init__(self, GM: float, integrator: str = 'RK4') -> None:
        ...
    def __repr__(self) -> str:
        ...
class ClassicalOE(OrbitState):
    M: float
    Omega: float
    a: float
    e: float
    i: float
    w: float
    @staticmethod
    def __init__(*args, **kwargs) -> None:
        ...
    def __repr__(self) -> str:
        ...
class EquinoctialOE(OrbitState):
    a: float
    h: float
    k: float
    lon: float
    p: float
    q: float
    def __init__(self, rv: numpy.ndarray[numpy.float64[6, 1]], frame: Frame = pylupnt.Frame.GCRF) -> None:
        ...
    def __repr__(self) -> str:
        ...
class Frame:
    """
    Members:
    
      ITRF
    
      ECEF
    
      GCRF
    
      ECI
    
      ICRF
    
      SER
    
      GSE
    
      EME
    
      MOD
    
      TOD
    
      EMR
    
      MOON_CI
    
      MOON_PA
    
      MOON_ME
    
      MOON_OP
    
      MARS_FIXED
    
      VENUS_FIXED
    """
    ECEF: typing.ClassVar[Frame]  # value = <Frame.ITRF: 0>
    ECI: typing.ClassVar[Frame]  # value = <Frame.ECI: 2>
    EME: typing.ClassVar[Frame]  # value = <Frame.ECI: 2>
    EMR: typing.ClassVar[Frame]  # value = <Frame.EMR: 8>
    GCRF: typing.ClassVar[Frame]  # value = <Frame.GCRF: 1>
    GSE: typing.ClassVar[Frame]  # value = <Frame.GSE: 5>
    ICRF: typing.ClassVar[Frame]  # value = <Frame.ICRF: 3>
    ITRF: typing.ClassVar[Frame]  # value = <Frame.ITRF: 0>
    MARS_FIXED: typing.ClassVar[Frame]  # value = <Frame.MARS_FIXED: 13>
    MOD: typing.ClassVar[Frame]  # value = <Frame.MOD: 6>
    MOON_CI: typing.ClassVar[Frame]  # value = <Frame.MOON_CI: 9>
    MOON_ME: typing.ClassVar[Frame]  # value = <Frame.MOON_ME: 11>
    MOON_OP: typing.ClassVar[Frame]  # value = <Frame.MOON_OP: 12>
    MOON_PA: typing.ClassVar[Frame]  # value = <Frame.MOON_PA: 10>
    SER: typing.ClassVar[Frame]  # value = <Frame.SER: 4>
    TOD: typing.ClassVar[Frame]  # value = <Frame.TOD: 7>
    VENUS_FIXED: typing.ClassVar[Frame]  # value = <Frame.VENUS_FIXED: 14>
    __members__: typing.ClassVar[dict[str, Frame]]  # value = {'ITRF': <Frame.ITRF: 0>, 'ECEF': <Frame.ITRF: 0>, 'GCRF': <Frame.GCRF: 1>, 'ECI': <Frame.ECI: 2>, 'ICRF': <Frame.ICRF: 3>, 'SER': <Frame.SER: 4>, 'GSE': <Frame.GSE: 5>, 'EME': <Frame.ECI: 2>, 'MOD': <Frame.MOD: 6>, 'TOD': <Frame.TOD: 7>, 'EMR': <Frame.EMR: 8>, 'MOON_CI': <Frame.MOON_CI: 9>, 'MOON_PA': <Frame.MOON_PA: 10>, 'MOON_ME': <Frame.MOON_ME: 11>, 'MOON_OP': <Frame.MOON_OP: 12>, 'MARS_FIXED': <Frame.MARS_FIXED: 13>, 'VENUS_FIXED': <Frame.VENUS_FIXED: 14>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class KeplerianDynamics:
    def __init__(self, arg0: float) -> None:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def propagate(self, state: ClassicalOE, dt: float) -> None:
        ...
    @typing.overload
    def propagate(self, state: QuasiNonsingOE, dt: float) -> None:
        ...
    @typing.overload
    def propagate(self, state: EquinoctialOE, dt: float) -> None:
        ...
    @typing.overload
    def propagate_with_stm(self, state: ClassicalOE, dt: float) -> numpy.ndarray[numpy.float64[6, 6]]:
        ...
    @typing.overload
    def propagate_with_stm(self, state: QuasiNonsingOE, dt: float) -> numpy.ndarray[numpy.float64[6, 6]]:
        ...
    @typing.overload
    def propagate_with_stm(self, state: EquinoctialOE, dt: float) -> numpy.ndarray[numpy.float64[6, 6]]:
        ...
class NBodyDynamics(NumericalOrbitDynamics):
    def __init__(self) -> None:
        ...
    def __repr__(self) -> str:
        ...
    def add_body(self, body: ...) -> None:
        ...
    def set_primary_body(self, body: ...) -> None:
        ...
class NaifId:
    """
    Members:
    
      SOLAR_SYSTEM_BARYCENTER
    
      SSB
    
      MERCURY_BARYCENTER
    
      VENUS_BARYCENTER
    
      EMB
    
      EARTH_MOON_BARYCENTER
    
      MARS_BARYCENTER
    
      JUPITER_BARYCENTER
    
      SATURN_BARYCENTER
    
      URANUS_BARYCENTER
    
      NEPTUNE_BARYCENTER
    
      PLUTO_BARYCENTER
    
      SUN
    
      MERCURY
    
      VENUS
    
      EARTH
    
      MOON
    
      MARS
    
      PHOBOS
    
      DEIMOS
    
      JUPITER
    """
    DEIMOS: typing.ClassVar[NaifId]  # value = <NaifId.DEIMOS: 402>
    EARTH: typing.ClassVar[NaifId]  # value = <NaifId.EARTH: 399>
    EARTH_MOON_BARYCENTER: typing.ClassVar[NaifId]  # value = <NaifId.EMB: 3>
    EMB: typing.ClassVar[NaifId]  # value = <NaifId.EMB: 3>
    JUPITER: typing.ClassVar[NaifId]  # value = <NaifId.JUPITER: 599>
    JUPITER_BARYCENTER: typing.ClassVar[NaifId]  # value = <NaifId.JUPITER_BARYCENTER: 5>
    MARS: typing.ClassVar[NaifId]  # value = <NaifId.MARS: 499>
    MARS_BARYCENTER: typing.ClassVar[NaifId]  # value = <NaifId.MARS_BARYCENTER: 4>
    MERCURY: typing.ClassVar[NaifId]  # value = <NaifId.MERCURY: 199>
    MERCURY_BARYCENTER: typing.ClassVar[NaifId]  # value = <NaifId.MERCURY_BARYCENTER: 1>
    MOON: typing.ClassVar[NaifId]  # value = <NaifId.MOON: 301>
    NEPTUNE_BARYCENTER: typing.ClassVar[NaifId]  # value = <NaifId.NEPTUNE_BARYCENTER: 8>
    PHOBOS: typing.ClassVar[NaifId]  # value = <NaifId.PHOBOS: 401>
    PLUTO_BARYCENTER: typing.ClassVar[NaifId]  # value = <NaifId.PLUTO_BARYCENTER: 9>
    SATURN_BARYCENTER: typing.ClassVar[NaifId]  # value = <NaifId.SATURN_BARYCENTER: 6>
    SOLAR_SYSTEM_BARYCENTER: typing.ClassVar[NaifId]  # value = <NaifId.SOLAR_SYSTEM_BARYCENTER: 0>
    SSB: typing.ClassVar[NaifId]  # value = <NaifId.SOLAR_SYSTEM_BARYCENTER: 0>
    SUN: typing.ClassVar[NaifId]  # value = <NaifId.SUN: 10>
    URANUS_BARYCENTER: typing.ClassVar[NaifId]  # value = <NaifId.URANUS_BARYCENTER: 7>
    VENUS: typing.ClassVar[NaifId]  # value = <NaifId.VENUS: 299>
    VENUS_BARYCENTER: typing.ClassVar[NaifId]  # value = <NaifId.VENUS_BARYCENTER: 2>
    __members__: typing.ClassVar[dict[str, NaifId]]  # value = {'SOLAR_SYSTEM_BARYCENTER': <NaifId.SOLAR_SYSTEM_BARYCENTER: 0>, 'SSB': <NaifId.SOLAR_SYSTEM_BARYCENTER: 0>, 'MERCURY_BARYCENTER': <NaifId.MERCURY_BARYCENTER: 1>, 'VENUS_BARYCENTER': <NaifId.VENUS_BARYCENTER: 2>, 'EMB': <NaifId.EMB: 3>, 'EARTH_MOON_BARYCENTER': <NaifId.EMB: 3>, 'MARS_BARYCENTER': <NaifId.MARS_BARYCENTER: 4>, 'JUPITER_BARYCENTER': <NaifId.JUPITER_BARYCENTER: 5>, 'SATURN_BARYCENTER': <NaifId.SATURN_BARYCENTER: 6>, 'URANUS_BARYCENTER': <NaifId.URANUS_BARYCENTER: 7>, 'NEPTUNE_BARYCENTER': <NaifId.NEPTUNE_BARYCENTER: 8>, 'PLUTO_BARYCENTER': <NaifId.PLUTO_BARYCENTER: 9>, 'SUN': <NaifId.SUN: 10>, 'MERCURY': <NaifId.MERCURY: 199>, 'VENUS': <NaifId.VENUS: 299>, 'EARTH': <NaifId.EARTH: 399>, 'MOON': <NaifId.MOON: 301>, 'MARS': <NaifId.MARS: 499>, 'PHOBOS': <NaifId.PHOBOS: 401>, 'DEIMOS': <NaifId.DEIMOS: 402>, 'JUPITER': <NaifId.JUPITER: 599>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class NumericalOrbitDynamics:
    @typing.overload
    def propagate(self, state: OrbitState, t0: float, tf: float, dt: float = 0.0) -> None:
        ...
    @typing.overload
    def propagate(self, state: numpy.ndarray[numpy.float64[6, 1]], t0: float, tf: float, dt: float = 0.0) -> numpy.ndarray[numpy.float64[6, 1]]:
        ...
    @typing.overload
    def propagate(self, state: numpy.ndarray[numpy.float64[6, 1]], t0: float, tfs: numpy.ndarray[numpy.float64[m, 1]], dt: float = 0.0, progress: bool = False) -> numpy.ndarray[numpy.float64[m, 1]]:
        ...
    @typing.overload
    def propagate_with_stm(self, state: CartesianOrbitState, t0: float, tf: float, dt: float) -> numpy.ndarray[numpy.float64[6, 6]]:
        ...
    @typing.overload
    def propagate_with_stm(self, state: numpy.ndarray[numpy.float64[6, 1]], t0: float, tf: float, dt: float) -> tuple[numpy.ndarray[numpy.float64[6, 1]], numpy.ndarray[numpy.float64[6, 6]]]:
        ...
    def set_time_step(self, arg0: float) -> None:
        ...
class OrbitState:
    frame: Frame
    state_repres: OrbitStateRepres
    vector: numpy.ndarray[numpy.float64[6, 1]]
    def __init__(self, vector: numpy.ndarray[numpy.float64[6, 1]], frame: Frame, state_repres: OrbitStateRepres, names: ..., std: ..., std: ..., std: ..., std: ..., std: ..., units: ..., std: ..., std: ..., std: ..., std: ..., std: ...) -> None:
        ...
    def __repr__(self) -> str:
        ...
    @property
    def names(self) -> ...:
        ...
    @property
    def size(self) -> int:
        ...
    @property
    def units(self) -> ...:
        ...
class OrbitStateRepres:
    """
    Members:
    
      CARTESIAN
    
      CLASSICAL_OE
    
      QUASI_NONSINGULAR_OE
    
      EQUINOTICAL_OE
    
      SINGULAR_ROE
    
      QUASINONSINGULAR_ROE
    
      DELAUNAY_OE
    """
    CARTESIAN: typing.ClassVar[OrbitStateRepres]  # value = <OrbitStateRepres.CARTESIAN: 0>
    CLASSICAL_OE: typing.ClassVar[OrbitStateRepres]  # value = <OrbitStateRepres.CLASSICAL_OE: 1>
    DELAUNAY_OE: typing.ClassVar[OrbitStateRepres]  # value = <OrbitStateRepres.DELAUNAY_OE: 6>
    EQUINOTICAL_OE: typing.ClassVar[OrbitStateRepres]  # value = <OrbitStateRepres.EQUINOTICAL_OE: 5>
    QUASINONSINGULAR_ROE: typing.ClassVar[OrbitStateRepres]  # value = <OrbitStateRepres.QUASINONSINGULAR_ROE: 9>
    QUASI_NONSINGULAR_OE: typing.ClassVar[OrbitStateRepres]  # value = <OrbitStateRepres.QUASI_NONSINGULAR_OE: 2>
    SINGULAR_ROE: typing.ClassVar[OrbitStateRepres]  # value = <OrbitStateRepres.SINGULAR_ROE: 3>
    __members__: typing.ClassVar[dict[str, OrbitStateRepres]]  # value = {'CARTESIAN': <OrbitStateRepres.CARTESIAN: 0>, 'CLASSICAL_OE': <OrbitStateRepres.CLASSICAL_OE: 1>, 'QUASI_NONSINGULAR_OE': <OrbitStateRepres.QUASI_NONSINGULAR_OE: 2>, 'EQUINOTICAL_OE': <OrbitStateRepres.EQUINOTICAL_OE: 5>, 'SINGULAR_ROE': <OrbitStateRepres.SINGULAR_ROE: 3>, 'QUASINONSINGULAR_ROE': <OrbitStateRepres.QUASINONSINGULAR_ROE: 9>, 'DELAUNAY_OE': <OrbitStateRepres.DELAUNAY_OE: 6>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class QuasiNonsingOE(OrbitState):
    Omega: float
    a: float
    ex: float
    ey: float
    i: float
    u: float
    def __init__(self, rv: numpy.ndarray[numpy.float64[6, 1]], frame: Frame = pylupnt.Frame.GCRF) -> None:
        ...
    def __repr__(self) -> str:
        ...
class QuasiNonsingROE(OrbitState):
    ada: float
    adex: float
    adey: float
    adix: float
    adiy: float
    adl: float
    def __init__(self, arg0: numpy.ndarray[numpy.float64[6, 1]], arg1: Frame) -> None:
        ...
    def __repr__(self) -> str:
        ...
class SingularROE(OrbitState):
    adM: float
    adOmega: float
    ada: float
    ade: float
    adi: float
    adw: float
    def __init__(self, arg0: numpy.ndarray[numpy.float64[6, 1]], arg1: Frame) -> None:
        ...
    def __repr__(self) -> str:
        ...
class SpiceInterface:
    @staticmethod
    def convert_time(t_tai: float, from: str, to: str) -> float:
        ...
    @staticmethod
    def extract_pck_coeffs() -> None:
        ...
    @staticmethod
    def get_body_pos_spice(t_tai: float, obs: NaifId, target: NaifId, ref_frame: Frame, ab_correction: str) -> numpy.ndarray[numpy.float64[3, 1]]:
        ...
    @staticmethod
    @typing.overload
    def get_body_pos_vel(t_tai: float, center: NaifId, target: NaifId) -> numpy.ndarray[numpy.float64[6, 1]]:
        ...
    @staticmethod
    @typing.overload
    def get_body_pos_vel(t_tai: numpy.ndarray[numpy.float64[m, 1]], center: NaifId, target: NaifId) -> numpy.ndarray[numpy.float64[m, 6]]:
        ...
    @staticmethod
    def get_frame_conversion_mat(t_tai: float, from: Frame, to: Frame) -> numpy.ndarray[numpy.float64[6, 6]]:
        ...
    @staticmethod
    def load_spice_kernel() -> None:
        ...
    @staticmethod
    def string_to_tai(gregorian_date: str) -> float:
        ...
    @staticmethod
    def string_to_tdb(gregorian_date: str) -> float:
        ...
    @staticmethod
    def tdb_to_string_utc(t_tdb: float, precision: int) -> str:
        ...
@typing.overload
def AzElRange2Cart(aer: numpy.ndarray[numpy.float64[3, 1]], xyz_ref: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def AzElRange2Cart(aer: numpy.ndarray[numpy.float64[1, 3]], xyz_ref: numpy.ndarray[numpy.float64[1, 3]]) -> numpy.ndarray[numpy.float64[1, 3]]:
    ...
@typing.overload
def AzElRange2Cart(aer: numpy.ndarray[numpy.float64[m, 3]], xyz_ref: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def AzElRange2Cart(aer: numpy.ndarray[numpy.float64[3, 1]], xyz_ref: numpy.ndarray[numpy.float64[m, 3]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def AzElRange2Cart(aer: numpy.ndarray[numpy.float64[m, 3]], xyz_ref: numpy.ndarray[numpy.float64[m, 3]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def Cart2AzElRange(xzy: numpy.ndarray[numpy.float64[3, 1]], xyz_ref: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def Cart2AzElRange(xzy: numpy.ndarray[numpy.float64[1, 3]], xyz_ref: numpy.ndarray[numpy.float64[1, 3]]) -> numpy.ndarray[numpy.float64[1, 3]]:
    ...
@typing.overload
def Cart2AzElRange(xzy: numpy.ndarray[numpy.float64[m, 3]], xyz_ref: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def Cart2AzElRange(xzy: numpy.ndarray[numpy.float64[3, 1]], xyz_ref: numpy.ndarray[numpy.float64[m, 3]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def Cart2AzElRange(xzy: numpy.ndarray[numpy.float64[m, 3]], xyz_ref: numpy.ndarray[numpy.float64[m, 3]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def Cart2EastNorthUp(xzy: numpy.ndarray[numpy.float64[3, 1]], xyz_ref: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def Cart2EastNorthUp(xzy: numpy.ndarray[numpy.float64[1, 3]], xyz_ref: numpy.ndarray[numpy.float64[1, 3]]) -> numpy.ndarray[numpy.float64[1, 3]]:
    ...
@typing.overload
def Cart2EastNorthUp(xzy: numpy.ndarray[numpy.float64[m, 3]], xyz_ref: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def Cart2EastNorthUp(xzy: numpy.ndarray[numpy.float64[3, 1]], xyz_ref: numpy.ndarray[numpy.float64[m, 3]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def Cart2EastNorthUp(xzy: numpy.ndarray[numpy.float64[m, 3]], xyz_ref: numpy.ndarray[numpy.float64[m, 3]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def Cart2LatLonAlt(r_cart: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def Cart2LatLonAlt(r_cart: numpy.ndarray[numpy.float64[m, 3]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def Decibel2Decimal(x: float) -> float:
    """
    Convert dB to decimal
    """
@typing.overload
def Decibel2Decimal(x: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1]]:
    ...
@typing.overload
def Decimal2Decibel(x: float) -> float:
    """
    Convert decimal to dB
    """
@typing.overload
def Decimal2Decibel(x: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1]]:
    ...
def DegMinSec2Degrees(hms: numpy.ndarray[numpy.float64[3, 1]]) -> float:
    """
    Convert degrees, minutes, seconds to degrees
    """
def Degrees2DegMinSec(deg: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
    Convert degrees to degrees, minutes, seconds
    """
@typing.overload
def EastNorthUp2Cart(enu: numpy.ndarray[numpy.float64[3, 1]], xyz_ref: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def EastNorthUp2Cart(enu: numpy.ndarray[numpy.float64[1, 3]], xyz_ref: numpy.ndarray[numpy.float64[1, 3]]) -> numpy.ndarray[numpy.float64[1, 3]]:
    ...
@typing.overload
def EastNorthUp2Cart(enu: numpy.ndarray[numpy.float64[m, 3]], xyz_ref: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def EastNorthUp2Cart(enu: numpy.ndarray[numpy.float64[3, 1]], xyz_ref: numpy.ndarray[numpy.float64[m, 3]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def EastNorthUp2Cart(enu: numpy.ndarray[numpy.float64[m, 3]], xyz_ref: numpy.ndarray[numpy.float64[m, 3]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def LatLonAlt2Cart(r_geo: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def LatLonAlt2Cart(r_geo: numpy.ndarray[numpy.float64[m, 3]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def Wrap2Pi(arg0: float) -> float:
    """
    Wrap angle to [-pi, pi]
    """
@typing.overload
def Wrap2Pi(arg0: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1]]:
    """
    Wrap angle to [-pi, pi]
    """
@typing.overload
def Wrap2TwoPi(arg0: float) -> float:
    """
    Wrap angle to [0, 2pi]
    """
@typing.overload
def Wrap2TwoPi(arg0: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1]]:
    """
    Wrap angle to [0, 2pi]
    """
@typing.overload
def cartesian_to_classical(cart: CartesianOrbitState, GM: float) -> ClassicalOE:
    ...
@typing.overload
def cartesian_to_classical(cart: numpy.ndarray[numpy.float64[6, 1]], GM: float) -> numpy.ndarray[numpy.float64[6, 1]]:
    ...
@typing.overload
def cartesian_to_classical(cart: numpy.ndarray[numpy.float64[6, 1]], GM: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def cartesian_to_classical(cart: numpy.ndarray[numpy.float64[m, 6]], GM: float) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def cartesian_to_classical(cart: numpy.ndarray[numpy.float64[m, 6]], GM: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_cartesian(coe: ClassicalOE, GM: float) -> CartesianOrbitState:
    ...
@typing.overload
def classical_to_cartesian(coe: numpy.ndarray[numpy.float64[6, 1]], GM: float) -> numpy.ndarray[numpy.float64[6, 1]]:
    ...
@typing.overload
def classical_to_cartesian(coe: numpy.ndarray[numpy.float64[6, 1]], GM: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_cartesian(coe: numpy.ndarray[numpy.float64[m, 6]], GM: float) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_cartesian(coe: numpy.ndarray[numpy.float64[m, 6]], GM: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_delaunay(arg0: ClassicalOE, arg1: float) -> ...:
    ...
@typing.overload
def classical_to_delaunay(coe: numpy.ndarray[numpy.float64[6, 1]], GM: float) -> numpy.ndarray[numpy.float64[6, 1]]:
    ...
@typing.overload
def classical_to_delaunay(coe: numpy.ndarray[numpy.float64[6, 1]], GM: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_delaunay(coe: numpy.ndarray[numpy.float64[m, 6]], GM: float) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_delaunay(coe: numpy.ndarray[numpy.float64[m, 6]], GM: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_equinoctial(arg0: ClassicalOE, arg1: float) -> EquinoctialOE:
    ...
@typing.overload
def classical_to_equinoctial(coe: numpy.ndarray[numpy.float64[6, 1]], GM: float) -> numpy.ndarray[numpy.float64[6, 1]]:
    ...
@typing.overload
def classical_to_equinoctial(coe: numpy.ndarray[numpy.float64[6, 1]], GM: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_equinoctial(coe: numpy.ndarray[numpy.float64[m, 6]], GM: float) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_equinoctial(coe: numpy.ndarray[numpy.float64[m, 6]], GM: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_quasi_nonsingular(coe: ClassicalOE, GM: float) -> QuasiNonsingOE:
    ...
@typing.overload
def classical_to_quasi_nonsingular(coe: numpy.ndarray[numpy.float64[6, 1]], GM: float) -> numpy.ndarray[numpy.float64[6, 1]]:
    ...
@typing.overload
def classical_to_quasi_nonsingular(coe: numpy.ndarray[numpy.float64[6, 1]], GM: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_quasi_nonsingular(coe: numpy.ndarray[numpy.float64[m, 6]], GM: float) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_quasi_nonsingular(coe: numpy.ndarray[numpy.float64[m, 6]], GM: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def compute_occultation(epoch: float, r1: numpy.ndarray[numpy.float64[3, 1]], r2: numpy.ndarray[numpy.float64[3, 1]], cs1: Frame, cs2: Frame, bodies: list[NaifId]) -> numpy.ndarray[numpy.float64[m, 1]]:
    """
    Compute occultation between two points
    """
@typing.overload
def compute_occultation(epoch: float, r1: numpy.ndarray[numpy.float64[m, 3]], r2: numpy.ndarray[numpy.float64[m, 3]], cs1: Frame, cs2: Frame, bodies: list[NaifId]) -> numpy.ndarray[numpy.float64[m, 1]]:
    """
    Compute occultation between two points
    """
@typing.overload
def compute_occultation(epoch: numpy.ndarray[numpy.float64[m, 1]], r1: numpy.ndarray[numpy.float64[m, 3]], r2: numpy.ndarray[numpy.float64[m, 3]], cs1: Frame, cs2: Frame, bodies: list[NaifId]) -> numpy.ndarray[numpy.float64[m, 1]]:
    """
    Compute occultation between two points
    """
@typing.overload
def convert_frame(t_tai: float, rv_in: numpy.ndarray[numpy.float64[6, 1]], frame_in: Frame, frame_out: Frame) -> numpy.ndarray[numpy.float64[6, 1]]:
    """
    Convert frame
    """
@typing.overload
def convert_frame(t_tai: float, r_in: numpy.ndarray[numpy.float64[3, 1]], frame_in: Frame, frame_out: Frame) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
    Convert frame
    """
@typing.overload
def convert_frame(t_tai: float, state_in: ..., frame_out: Frame) -> ...:
    """
    Convert frame
    """
@typing.overload
def convert_frame(t_tai: numpy.ndarray[numpy.float64[m, 1]], rv_in: numpy.ndarray[numpy.float64[6, 1]], frame_in: Frame, frame_out: Frame) -> numpy.ndarray[numpy.float64[m, 6]]:
    """
    Convert frame
    """
@typing.overload
def convert_frame(t_tai: numpy.ndarray[numpy.float64[m, 1]], r_in: numpy.ndarray[numpy.float64[3, 1]], frame_in: Frame, frame_out: Frame) -> numpy.ndarray[numpy.float64[m, 3]]:
    """
    Convert frame
    """
@typing.overload
def convert_frame(t_tai: float, rv_in: numpy.ndarray[numpy.float64[m, 6]], frame_in: Frame, frame_out: Frame) -> numpy.ndarray[numpy.float64[m, 6]]:
    """
    Convert frame
    """
@typing.overload
def convert_frame(t_tai: float, r_in: numpy.ndarray[numpy.float64[m, 3]], frame_in: Frame, frame_out: Frame) -> numpy.ndarray[numpy.float64[m, 3]]:
    """
    Convert frame
    """
@typing.overload
def convert_frame(t_tai: numpy.ndarray[numpy.float64[m, 1]], rv_in: numpy.ndarray[numpy.float64[m, 6]], frame_in: Frame, frame_out: Frame) -> numpy.ndarray[numpy.float64[m, 6]]:
    """
    Convert frame
    """
@typing.overload
def convert_frame(t_tai: numpy.ndarray[numpy.float64[m, 1]], r_in: numpy.ndarray[numpy.float64[m, 3]], frame_in: Frame, frame_out: Frame) -> numpy.ndarray[numpy.float64[m, 3]]:
    """
    Convert frame
    """
@typing.overload
def convert_orbit_state(arg0: numpy.ndarray[numpy.float64[6, 1]], arg1: OrbitStateRepres, arg2: OrbitStateRepres, arg3: float) -> numpy.ndarray[numpy.float64[6, 1]]:
    ...
@typing.overload
def convert_orbit_state(arg0: numpy.ndarray[numpy.float64[6, 1]], arg1: numpy.ndarray[..., double>[6, 1]], arg2: OrbitStateRepres, arg3: OrbitStateRepres, arg4: OrbitStateRepres, arg5: float) -> numpy.ndarray[numpy.float64[6, 1]]:
    ...
@typing.overload
def delaunay_to_classical(deloe: ..., GM: float) -> ClassicalOE:
    ...
@typing.overload
def delaunay_to_classical(deloe: numpy.ndarray[numpy.float64[6, 1]], GM: float) -> numpy.ndarray[numpy.float64[6, 1]]:
    ...
@typing.overload
def delaunay_to_classical(deloe: numpy.ndarray[numpy.float64[6, 1]], GM: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def delaunay_to_classical(deloe: numpy.ndarray[numpy.float64[m, 6]], GM: float) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def delaunay_to_classical(deloe: numpy.ndarray[numpy.float64[m, 6]], GM: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def eccentric_to_mean_anomaly(E: float, e: float) -> float:
    ...
@typing.overload
def eccentric_to_mean_anomaly(E: numpy.ndarray[numpy.float64[m, 1]], e: float) -> numpy.ndarray[numpy.float64[m, 1]]:
    ...
@typing.overload
def eccentric_to_mean_anomaly(E: float, e: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1]]:
    ...
@typing.overload
def eccentric_to_mean_anomaly(E: numpy.ndarray[numpy.float64[m, 1]], e: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1]]:
    ...
@typing.overload
def eccentric_to_true_anomaly(E: float, e: float) -> float:
    ...
@typing.overload
def eccentric_to_true_anomaly(E: numpy.ndarray[numpy.float64[m, 1]], e: float) -> numpy.ndarray[numpy.float64[m, 1]]:
    ...
@typing.overload
def eccentric_to_true_anomaly(E: float, e: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1]]:
    ...
@typing.overload
def eccentric_to_true_anomaly(E: numpy.ndarray[numpy.float64[m, 1]], e: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1]]:
    ...
@typing.overload
def equinoctial_to_classical(eqoe: EquinoctialOE, GM: float) -> ClassicalOE:
    ...
@typing.overload
def equinoctial_to_classical(eqoe: numpy.ndarray[numpy.float64[6, 1]], GM: float) -> numpy.ndarray[numpy.float64[6, 1]]:
    ...
@typing.overload
def equinoctial_to_classical(eqoe: numpy.ndarray[numpy.float64[6, 1]], GM: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def equinoctial_to_classical(eqoe: numpy.ndarray[numpy.float64[m, 6]], GM: float) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def equinoctial_to_classical(eqoe: numpy.ndarray[numpy.float64[m, 6]], GM: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def mean_to_eccentric_anomaly(M: float, e: float) -> float:
    ...
@typing.overload
def mean_to_eccentric_anomaly(M: numpy.ndarray[numpy.float64[m, 1]], e: float) -> numpy.ndarray[numpy.float64[m, 1]]:
    ...
@typing.overload
def mean_to_eccentric_anomaly(M: float, e: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1]]:
    ...
@typing.overload
def mean_to_eccentric_anomaly(M: numpy.ndarray[numpy.float64[m, 1]], e: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1]]:
    ...
@typing.overload
def mean_to_true_anomaly(M: float, e: float) -> float:
    ...
@typing.overload
def mean_to_true_anomaly(M: numpy.ndarray[numpy.float64[m, 1]], e: float) -> numpy.ndarray[numpy.float64[m, 1]]:
    ...
@typing.overload
def mean_to_true_anomaly(M: float, e: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1]]:
    ...
@typing.overload
def mean_to_true_anomaly(M: numpy.ndarray[numpy.float64[m, 1]], e: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1]]:
    ...
@typing.overload
def quasi_nonsingular_to_classical(arg0: QuasiNonsingOE, arg1: float) -> ClassicalOE:
    ...
@typing.overload
def quasi_nonsingular_to_classical(qnsoe: numpy.ndarray[numpy.float64[6, 1]], GM: float) -> numpy.ndarray[numpy.float64[6, 1]]:
    ...
@typing.overload
def quasi_nonsingular_to_classical(qnsoe: numpy.ndarray[numpy.float64[6, 1]], GM: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def quasi_nonsingular_to_classical(qnsoe: numpy.ndarray[numpy.float64[m, 6]], GM: float) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def quasi_nonsingular_to_classical(qnsoe: numpy.ndarray[numpy.float64[m, 6]], GM: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def relative_quasi_nonsingular_to_classical(coe: ClassicalOE, rel_qnsoe: QuasiNonsingROE) -> ClassicalOE:
    ...
@typing.overload
def relative_quasi_nonsingular_to_classical(coe: numpy.ndarray[numpy.float64[6, 1]], rel_qnsoe: numpy.ndarray[numpy.float64[6, 1]]) -> numpy.ndarray[numpy.float64[6, 1]]:
    ...
@typing.overload
def relative_quasi_nonsingular_to_classical(coe: numpy.ndarray[numpy.float64[1, 6]], rel_qnsoe: numpy.ndarray[numpy.float64[1, 6]]) -> numpy.ndarray[numpy.float64[1, 6]]:
    ...
@typing.overload
def relative_quasi_nonsingular_to_classical(coe: numpy.ndarray[numpy.float64[m, 6]], rel_qnsoe: numpy.ndarray[numpy.float64[6, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def relative_quasi_nonsingular_to_classical(coe: numpy.ndarray[numpy.float64[6, 1]], rel_qnsoe: numpy.ndarray[numpy.float64[m, 6]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def relative_quasi_nonsingular_to_classical(coe: numpy.ndarray[numpy.float64[m, 6]], rel_qnsoe: numpy.ndarray[numpy.float64[m, 6]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def true_to_eccentric_anomaly(nu: float, e: float) -> float:
    ...
@typing.overload
def true_to_eccentric_anomaly(nu: numpy.ndarray[numpy.float64[m, 1]], e: float) -> numpy.ndarray[numpy.float64[m, 1]]:
    ...
@typing.overload
def true_to_eccentric_anomaly(nu: float, e: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1]]:
    ...
@typing.overload
def true_to_eccentric_anomaly(nu: numpy.ndarray[numpy.float64[m, 1]], e: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1]]:
    ...
@typing.overload
def true_to_mean_anomaly(f: float, e: float) -> float:
    ...
@typing.overload
def true_to_mean_anomaly(f: numpy.ndarray[numpy.float64[m, 1]], e: float) -> numpy.ndarray[numpy.float64[m, 1]]:
    ...
@typing.overload
def true_to_mean_anomaly(f: float, e: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1]]:
    ...
@typing.overload
def true_to_mean_anomaly(f: numpy.ndarray[numpy.float64[m, 1]], e: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1]]:
    ...
A1_TAI_OFFSET: float = 0.0343817
ARCSEC_DEG: float = 3600.0
ARCSEC_RAD: float = 206264.80624709636
CARTESIAN: OrbitStateRepres  # value = <OrbitStateRepres.CARTESIAN: 0>
CLASSICAL_OE: OrbitStateRepres  # value = <OrbitStateRepres.CLASSICAL_OE: 1>
DAYS_CENTURY: float = 36525.0
DAYS_SEC: float = 1.1574074074074073e-05
DAYS_WEEK: float = 7.0
DAYS_YEAR: float = 365.25
DEG: float = 57.29577951308232
DEG_ARCSEC: float = 0.0002777777777777778
DEIMOS: NaifId  # value = <NaifId.DEIMOS: 402>
DELAUNAY_OE: OrbitStateRepres  # value = <OrbitStateRepres.DELAUNAY_OE: 6>
E: float = 2.718281828459045
EARTH: NaifId  # value = <NaifId.EARTH: 399>
EARTH_MOON_BARYCENTER: NaifId  # value = <NaifId.EMB: 3>
ECEF: Frame  # value = <Frame.ITRF: 0>
ECI: Frame  # value = <Frame.ECI: 2>
EMB: NaifId  # value = <NaifId.EMB: 3>
EME: Frame  # value = <Frame.ECI: 2>
EMR: Frame  # value = <Frame.EMR: 8>
EPS: float = 1e-16
EQUINOTICAL_OE: OrbitStateRepres  # value = <OrbitStateRepres.EQUINOTICAL_OE: 5>
FOOT_M: float = 0.3048
GCRF: Frame  # value = <Frame.GCRF: 1>
GM_EARTH: float = 398600.435507
GM_MERCURY: float = 22031.868551
GM_MOON: float = 4902.800118
GM_SUN: float = 132712440041.27942
GM_VENUS: float = 324858.592
GPS: str = 'GPS'
GSE: Frame  # value = <Frame.GSE: 5>
HOURS_DAY: float = 24.0
ICRF: Frame  # value = <Frame.ICRF: 3>
INCH_M: float = 0.0254
ITRF: Frame  # value = <Frame.ITRF: 0>
JD_J2000: float = 2451545.0
JD_JAN_5_1941: float = 2430000.0
JD_MJD_OFFSET: float = 2400000.5
JD_NOV_17_1858: float = 2400000.5
JD_T0: float = 2443144.5003725
JD_TDB: str = 'JDTDB'
JD_TT: str = 'JDTDT'
JULIAN_DATE_OF_010541: int = 2430000
JUPITER: NaifId  # value = <NaifId.JUPITER: 599>
JUPITER_BARYCENTER: NaifId  # value = <NaifId.JUPITER_BARYCENTER: 5>
KM_M: float = 0.001
LBM_TO_KG: float = 0.45359237
L_B: float = 1.550505e-08
L_G: float = 6.969290134e-10
MARS: NaifId  # value = <NaifId.MARS: 499>
MARS_BARYCENTER: NaifId  # value = <NaifId.MARS_BARYCENTER: 4>
MARS_FIXED: Frame  # value = <Frame.MARS_FIXED: 13>
MERCURY: NaifId  # value = <NaifId.MERCURY: 199>
MERCURY_BARYCENTER: NaifId  # value = <NaifId.MERCURY_BARYCENTER: 1>
MILE_M: float = 1609.344
MINS_DAY: float = 1440.0
MINS_HOUR: float = 60.0
MJD_J2000: float = 51544.5
MOD: Frame  # value = <Frame.MOD: 6>
MOON: NaifId  # value = <NaifId.MOON: 301>
MOON_CI: Frame  # value = <Frame.MOON_CI: 9>
MOON_ME: Frame  # value = <Frame.MOON_ME: 11>
MOON_OP: Frame  # value = <Frame.MOON_OP: 12>
MOON_PA: Frame  # value = <Frame.MOON_PA: 10>
M_KM: float = 1000.0
NEPTUNE_BARYCENTER: NaifId  # value = <NaifId.NEPTUNE_BARYCENTER: 8>
NUM_SECS: float = 86400.0
PHOBOS: NaifId  # value = <NaifId.PHOBOS: 401>
PI: float = 3.141592653589793
PI_OVER_TWO: float = 1.5707963267948966
PLUTO_BARYCENTER: NaifId  # value = <NaifId.PLUTO_BARYCENTER: 9>
QUASINONSINGULAR_ROE: OrbitStateRepres  # value = <OrbitStateRepres.QUASINONSINGULAR_ROE: 9>
QUASI_NONSINGULAR_OE: OrbitStateRepres  # value = <OrbitStateRepres.QUASI_NONSINGULAR_OE: 2>
RAD: float = 0.017453292519943295
RAD_ARCSEC: float = 4.84813681109536e-06
R_EARTH: float = 6378.137
R_MOON: float = 1737.4
SATURN_BARYCENTER: NaifId  # value = <NaifId.SATURN_BARYCENTER: 6>
SECS_DAY: float = 86400.0
SECS_HOUR: float = 3600.0
SECS_MINUTE: float = 60.0
SER: Frame  # value = <Frame.SER: 4>
SINGULAR_ROE: OrbitStateRepres  # value = <OrbitStateRepres.SINGULAR_ROE: 3>
SLUG_TO_KG: float = 14.59390294
SOLAR_SYSTEM_BARYCENTER: NaifId  # value = <NaifId.SOLAR_SYSTEM_BARYCENTER: 0>
SSB: NaifId  # value = <NaifId.SOLAR_SYSTEM_BARYCENTER: 0>
SUN: NaifId  # value = <NaifId.SUN: 10>
TAI: str = 'TAI'
TCB: str = 'TCB'
TCG: str = 'TCG'
TDB: str = 'TDB'
TIME_OF_J2000: float = 883655990.85
TOD: Frame  # value = <Frame.TOD: 7>
TT: str = 'TT'
TT_TAI_OFFSET: float = 32.184
TWO_PI: float = 6.283185307179586
URANUS_BARYCENTER: NaifId  # value = <NaifId.URANUS_BARYCENTER: 7>
UT1: str = 'UT1'
UTC: str = 'UTC'
VENUS: NaifId  # value = <NaifId.VENUS: 299>
VENUS_BARYCENTER: NaifId  # value = <NaifId.VENUS_BARYCENTER: 2>
VENUS_FIXED: Frame  # value = <Frame.VENUS_FIXED: 14>

from __future__ import annotations
import numpy
import pylupnt
import typing
__all__ = ['A1MJD_OF_J2000', 'A1_TAI_OFFSET', 'AU', 'Body', 'C', 'C22_MOON', 'CARTESIAN', 'CLASSICAL_OE', 'CartesianOrbitState', 'CartesianTwoBodyDynamics', 'ClassicalOE', 'CoordConverter', 'DAYS_PER_JULIAN_CENTURY', 'DAYS_PER_SEC', 'DAYS_PER_YEAR', 'DEG_PER_RAD', 'DEIMOS', 'DELAUNAY_OE', 'E', 'EARTH', 'EARTH_BARYCENTER', 'EARTH_MOON_BARYCENTER', 'ECEF', 'ECI', 'EME', 'EMR', 'EQUINOTICAL_OE', 'EquinoctialOE', 'Frame', 'GCRF', 'GSE', 'ICRF', 'ITRF', 'J2_EARTH', 'J2_MOON', 'JD_JAN_5_1941', 'JD_MJD_OFFSET', 'JD_NOV_17_1858', 'JD_OF_J2000', 'JUPITER', 'JUPITER_BARYCENTER', 'KeplerianDynamics', 'MARS', 'MARS_BARYCENTER', 'ME', 'MERCURY', 'MERCURY_BARYCENTER', 'MI', 'MJD_OF_J2000', 'MOD', 'MOON', 'MU_EARTH', 'MU_MOON', 'NBodyDynamics', 'NEPTUNE_BARYCENTER', 'NaifId', 'NumericalOrbitDynamics', 'OMEGA_E_M', 'OP', 'OrbitState', 'OrbitStateRepres', 'PA', 'PHOBOS', 'PI', 'PI_OVER_TWO', 'PLUTO_BARYCENTER', 'P_SUN', 'QUASINONSINGULAR_ROE', 'QUASI_NONSINGULAR_OE', 'QuasiNonsingularOE', 'QuasiNonsingularROE', 'RAD_PER_DEG', 'R_EARTH', 'R_MOON', 'SATURN_BARYCENTER', 'SECS_PER_DAY', 'SECS_PER_HOUR', 'SECS_PER_MINUTE', 'SER', 'SINGULAR_ROE', 'SOLAR_SYSTEM_BARYCENTER', 'SUN', 'S_AU', 'SingularROE', 'SpiceInterface', 'TIME_OF_J2000', 'TOD', 'TT_TAI_OFFSET', 'TWO_PI', 'URANUS_BARYCENTER', 'VENUS', 'VENUS_BARYCENTER', 'azimuth_elevation_range_to_cartesian', 'cartesian_to_azimuth_elevation_range', 'cartesian_to_classical', 'cartesian_to_east_north_up', 'cartesian_to_geographical', 'cartesian_to_spherical', 'classical_to_cartesian', 'classical_to_delaunay', 'classical_to_equinoctial', 'classical_to_quasi_nonsingular', 'compute_occultation', 'convert_orbit_state', 'dB2decimal', 'd_E_EMB', 'd_E_M', 'd_M_EMB', 'decimal2dB', 'degrees2dms', 'delaunay_to_classical', 'dms2degrees', 'east_north_up_to_cartesian', 'eccentric_to_mean_anomaly', 'eccentric_to_true_anomaly', 'equinoctial_to_classical', 'geographical_to_cartesian', 'mean_to_eccentric_anomaly', 'mean_to_true_anomaly', 'quasi_nonsingular_to_classical', 'relative_quasi_nonsingular_to_classical', 'spherical_to_cartesian', 'true_to_eccentric_anomaly', 'true_to_mean_anomaly', 'wrapTo2Pi', 'wrapToPi']
class Body:
    @staticmethod
    def Earth(n_max: int = 0, m_max: int = 0) -> Body:
        ...
    @staticmethod
    def Moon(n_max: int = 0, m_max: int = 0) -> Body:
        ...
    def __init__(self) -> None:
        ...
    def __repr__(self) -> str:
        ...
class CartesianOrbitState(OrbitState):
    r: numpy.ndarray[numpy.float64[3, 1]]
    v: numpy.ndarray[numpy.float64[3, 1]]
    def __init__(self, rv: numpy.ndarray[numpy.float64[6, 1]], coord_sys: Frame = pylupnt.Frame.MI) -> None:
        ...
    def __repr__(self) -> str:
        ...
class CartesianTwoBodyDynamics(NumericalOrbitDynamics):
    def __init__(self, mu: float, integrator: str = 'RK4') -> None:
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
class CoordConverter:
    @staticmethod
    @typing.overload
    def convert(epoch: float, rv_in: numpy.ndarray[numpy.float64[6, 1]], coord_sys_in: Frame, coord_sys_out: Frame) -> numpy.ndarray[numpy.float64[6, 1]]:
        """
        Convert coordinate system
        """
    @staticmethod
    @typing.overload
    def convert(epoch: float, r_in: numpy.ndarray[numpy.float64[3, 1]], coord_sys_in: Frame, coord_sys_out: Frame) -> numpy.ndarray[numpy.float64[3, 1]]:
        """
        Convert coordinate system
        """
    @staticmethod
    @typing.overload
    def convert(epoch: float, state_in: ..., coord_sys_out: Frame) -> ...:
        """
        Convert coordinate system
        """
    @staticmethod
    @typing.overload
    def convert(epoch: numpy.ndarray[numpy.float64[m, 1]], rv_in: numpy.ndarray[numpy.float64[6, 1]], coord_sys_in: Frame, coord_sys_out: Frame) -> numpy.ndarray[numpy.float64[m, n]]:
        """
        Convert coordinate system
        """
    @staticmethod
    @typing.overload
    def convert(epoch: numpy.ndarray[numpy.float64[m, 1]], r_in: numpy.ndarray[numpy.float64[3, 1]], coord_sys_in: Frame, coord_sys_out: Frame) -> numpy.ndarray[numpy.float64[m, n]]:
        """
        Convert coordinate system
        """
    @staticmethod
    @typing.overload
    def convert(epoch: float, rv_in: numpy.ndarray[numpy.float64[m, 6]], coord_sys_in: Frame, coord_sys_out: Frame) -> numpy.ndarray[numpy.float64[m, 6]]:
        """
        Convert coordinate system
        """
    @staticmethod
    @typing.overload
    def convert(epoch: float, r_in: numpy.ndarray[numpy.float64[m, 3]], coord_sys_in: Frame, coord_sys_out: Frame) -> numpy.ndarray[numpy.float64[m, 3]]:
        """
        Convert coordinate system
        """
    @staticmethod
    @typing.overload
    def convert(epoch: numpy.ndarray[numpy.float64[m, 1]], rv_in: numpy.ndarray[numpy.float64[m, 6]], coord_sys_in: Frame, coord_sys_out: Frame) -> numpy.ndarray[numpy.float64[m, 6]]:
        """
        Convert coordinate system
        """
    @staticmethod
    @typing.overload
    def convert(epoch: numpy.ndarray[numpy.float64[m, 1]], r_in: numpy.ndarray[numpy.float64[m, 3]], coord_sys_in: Frame, coord_sys_out: Frame) -> numpy.ndarray[numpy.float64[m, 3]]:
        """
        Convert coordinate system
        """
class EquinoctialOE(OrbitState):
    a: float
    h: float
    k: float
    lon: float
    p: float
    q: float
    def __init__(self, arg0: numpy.ndarray[numpy.float64[6, 1]], arg1: Frame) -> None:
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
    
      MI
    
      PA
    
      ME
    
      OP
    """
    ECEF: typing.ClassVar[Frame]  # value = <Frame.ITRF: 0>
    ECI: typing.ClassVar[Frame]  # value = <Frame.GCRF: 1>
    EME: typing.ClassVar[Frame]  # value = <Frame.EME: 5>
    EMR: typing.ClassVar[Frame]  # value = <Frame.EMR: 8>
    GCRF: typing.ClassVar[Frame]  # value = <Frame.GCRF: 1>
    GSE: typing.ClassVar[Frame]  # value = <Frame.GSE: 4>
    ICRF: typing.ClassVar[Frame]  # value = <Frame.ICRF: 2>
    ITRF: typing.ClassVar[Frame]  # value = <Frame.ITRF: 0>
    ME: typing.ClassVar[Frame]  # value = <Frame.ME: 11>
    MI: typing.ClassVar[Frame]  # value = <Frame.MI: 9>
    MOD: typing.ClassVar[Frame]  # value = <Frame.MOD: 6>
    OP: typing.ClassVar[Frame]  # value = <Frame.OP: 12>
    PA: typing.ClassVar[Frame]  # value = <Frame.PA: 10>
    SER: typing.ClassVar[Frame]  # value = <Frame.SER: 3>
    TOD: typing.ClassVar[Frame]  # value = <Frame.TOD: 7>
    __members__: typing.ClassVar[dict[str, Frame]]  # value = {'ITRF': <Frame.ITRF: 0>, 'ECEF': <Frame.ITRF: 0>, 'GCRF': <Frame.GCRF: 1>, 'ECI': <Frame.GCRF: 1>, 'ICRF': <Frame.ICRF: 2>, 'SER': <Frame.SER: 3>, 'GSE': <Frame.GSE: 4>, 'EME': <Frame.EME: 5>, 'MOD': <Frame.MOD: 6>, 'TOD': <Frame.TOD: 7>, 'EMR': <Frame.EMR: 8>, 'MI': <Frame.MI: 9>, 'PA': <Frame.PA: 10>, 'ME': <Frame.ME: 11>, 'OP': <Frame.OP: 12>}
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
    def propagate(self, state: QuasiNonsingularOE, dt: float) -> None:
        ...
    @typing.overload
    def propagate(self, state: EquinoctialOE, dt: float) -> None:
        ...
    @typing.overload
    def propagate_with_stm(self, state: ClassicalOE, dt: float) -> numpy.ndarray[numpy.float64[6, 6]]:
        ...
    @typing.overload
    def propagate_with_stm(self, state: QuasiNonsingularOE, dt: float) -> numpy.ndarray[numpy.float64[6, 6]]:
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
    def set_time_step(self, dt: float) -> None:
        ...
class NaifId:
    """
    Members:
    
      SOLAR_SYSTEM_BARYCENTER
    
      MERCURY_BARYCENTER
    
      VENUS_BARYCENTER
    
      EARTH_BARYCENTER
    
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
    EARTH_BARYCENTER: typing.ClassVar[NaifId]  # value = <NaifId.EARTH_BARYCENTER: 3>
    EARTH_MOON_BARYCENTER: typing.ClassVar[NaifId]  # value = <NaifId.EARTH_BARYCENTER: 3>
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
    SUN: typing.ClassVar[NaifId]  # value = <NaifId.SUN: 10>
    URANUS_BARYCENTER: typing.ClassVar[NaifId]  # value = <NaifId.URANUS_BARYCENTER: 7>
    VENUS: typing.ClassVar[NaifId]  # value = <NaifId.VENUS: 299>
    VENUS_BARYCENTER: typing.ClassVar[NaifId]  # value = <NaifId.VENUS_BARYCENTER: 2>
    __members__: typing.ClassVar[dict[str, NaifId]]  # value = {'SOLAR_SYSTEM_BARYCENTER': <NaifId.SOLAR_SYSTEM_BARYCENTER: 0>, 'MERCURY_BARYCENTER': <NaifId.MERCURY_BARYCENTER: 1>, 'VENUS_BARYCENTER': <NaifId.VENUS_BARYCENTER: 2>, 'EARTH_BARYCENTER': <NaifId.EARTH_BARYCENTER: 3>, 'EARTH_MOON_BARYCENTER': <NaifId.EARTH_BARYCENTER: 3>, 'MARS_BARYCENTER': <NaifId.MARS_BARYCENTER: 4>, 'JUPITER_BARYCENTER': <NaifId.JUPITER_BARYCENTER: 5>, 'SATURN_BARYCENTER': <NaifId.SATURN_BARYCENTER: 6>, 'URANUS_BARYCENTER': <NaifId.URANUS_BARYCENTER: 7>, 'NEPTUNE_BARYCENTER': <NaifId.NEPTUNE_BARYCENTER: 8>, 'PLUTO_BARYCENTER': <NaifId.PLUTO_BARYCENTER: 9>, 'SUN': <NaifId.SUN: 10>, 'MERCURY': <NaifId.MERCURY: 199>, 'VENUS': <NaifId.VENUS: 299>, 'EARTH': <NaifId.EARTH: 399>, 'MOON': <NaifId.MOON: 301>, 'MARS': <NaifId.MARS: 499>, 'PHOBOS': <NaifId.PHOBOS: 401>, 'DEIMOS': <NaifId.DEIMOS: 402>, 'JUPITER': <NaifId.JUPITER: 599>}
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
    def propagate(self, state: OrbitState, t0: float, tf: float, dt: float) -> None:
        ...
    @typing.overload
    def propagate(self, state: numpy.ndarray[numpy.float64[6, 1]], t0: float, tf: float, dt: float) -> numpy.ndarray[numpy.float64[6, 1]]:
        ...
    @typing.overload
    def propagate(self, state: numpy.ndarray[numpy.float64[6, 1]], t0: float, tfs: numpy.ndarray[numpy.float64[m, 1]], progress: bool = False) -> numpy.ndarray[numpy.float64[m, n]]:
        ...
    @typing.overload
    def propagate_with_stm(self, state: CartesianOrbitState, t0: float, tf: float, dt: float) -> numpy.ndarray[numpy.float64[6, 6]]:
        ...
    @typing.overload
    def propagate_with_stm(self, state: numpy.ndarray[numpy.float64[6, 1]], t0: float, tf: float, dt: float) -> tuple[numpy.ndarray[numpy.float64[6, 1]], numpy.ndarray[numpy.float64[6, 6]]]:
        ...
class OrbitState:
    coord_sys: Frame
    state_repres: OrbitStateRepres
    vector: numpy.ndarray[numpy.float64[6, 1]]
    def __init__(self, vector: numpy.ndarray[numpy.float64[6, 1]], coord_sys: Frame, state_repres: OrbitStateRepres, names: ..., std: ..., std: ..., std: ..., std: ..., std: ..., units: ..., std: ..., std: ..., std: ..., std: ..., std: ...) -> None:
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
class QuasiNonsingularOE(OrbitState):
    Omega: float
    a: float
    ex: float
    ey: float
    i: float
    u: float
    def __init__(self, arg0: numpy.ndarray[numpy.float64[6, 1]], arg1: Frame) -> None:
        ...
    def __repr__(self) -> str:
        ...
class QuasiNonsingularROE(OrbitState):
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
    def convert_time(arg0: float, arg1: str, arg2: str) -> float:
        ...
    @staticmethod
    def extract_pck_coeffs() -> None:
        ...
    @staticmethod
    def get_body_pos(*args, **kwargs) -> numpy.ndarray[numpy.float64[3, 1]]:
        ...
    @staticmethod
    @typing.overload
    def get_body_pos_vel(arg0: float, arg1: NaifId, arg2: NaifId) -> numpy.ndarray[numpy.float64[6, 1]]:
        ...
    @staticmethod
    @typing.overload
    def get_body_pos_vel(arg0: numpy.ndarray[numpy.float64[m, 1]], arg1: NaifId, arg2: NaifId) -> numpy.ndarray[numpy.float64[m, 6]]:
        ...
    @staticmethod
    def get_frame_conversion_matrix(arg0: float, arg1: Frame, arg2: Frame) -> numpy.ndarray[numpy.float64[6, 6]]:
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
    def tdb_to_string_utc(arg0: float, arg1: int) -> str:
        ...
@typing.overload
def azimuth_elevation_range_to_cartesian(r_aer_ref: numpy.ndarray[numpy.float64[3, 1]], r_aer: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def azimuth_elevation_range_to_cartesian(r_aer_ref: numpy.ndarray[numpy.float64[1, 3]], r_aer: numpy.ndarray[numpy.float64[1, 3]]) -> numpy.ndarray[numpy.float64[1, 3]]:
    ...
@typing.overload
def azimuth_elevation_range_to_cartesian(r_aer_ref: numpy.ndarray[numpy.float64[m, 3]], r_aer: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def azimuth_elevation_range_to_cartesian(r_aer_ref: numpy.ndarray[numpy.float64[3, 1]], r_aer: numpy.ndarray[numpy.float64[m, 3]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def azimuth_elevation_range_to_cartesian(r_aer_ref: numpy.ndarray[numpy.float64[m, 3]], r_aer: numpy.ndarray[numpy.float64[m, 3]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def cartesian_to_azimuth_elevation_range(r_cart_ref: numpy.ndarray[numpy.float64[3, 1]], r_cart: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def cartesian_to_azimuth_elevation_range(r_cart_ref: numpy.ndarray[numpy.float64[1, 3]], r_cart: numpy.ndarray[numpy.float64[1, 3]]) -> numpy.ndarray[numpy.float64[1, 3]]:
    ...
@typing.overload
def cartesian_to_azimuth_elevation_range(r_cart_ref: numpy.ndarray[numpy.float64[m, 3]], r_cart: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def cartesian_to_azimuth_elevation_range(r_cart_ref: numpy.ndarray[numpy.float64[3, 1]], r_cart: numpy.ndarray[numpy.float64[m, 3]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def cartesian_to_azimuth_elevation_range(r_cart_ref: numpy.ndarray[numpy.float64[m, 3]], r_cart: numpy.ndarray[numpy.float64[m, 3]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def cartesian_to_classical(cart: CartesianOrbitState, mu: float) -> ClassicalOE:
    ...
@typing.overload
def cartesian_to_classical(cart: numpy.ndarray[numpy.float64[6, 1]], mu: float) -> numpy.ndarray[numpy.float64[6, 1]]:
    ...
@typing.overload
def cartesian_to_classical(cart: numpy.ndarray[numpy.float64[6, 1]], mu: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def cartesian_to_classical(cart: numpy.ndarray[numpy.float64[m, 6]], mu: float) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def cartesian_to_classical(cart: numpy.ndarray[numpy.float64[m, 6]], mu: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def cartesian_to_east_north_up(r_ref: numpy.ndarray[numpy.float64[3, 1]], r_cart: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def cartesian_to_east_north_up(r_ref: numpy.ndarray[numpy.float64[1, 3]], r_cart: numpy.ndarray[numpy.float64[1, 3]]) -> numpy.ndarray[numpy.float64[1, 3]]:
    ...
@typing.overload
def cartesian_to_east_north_up(r_ref: numpy.ndarray[numpy.float64[m, 3]], r_cart: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def cartesian_to_east_north_up(r_ref: numpy.ndarray[numpy.float64[3, 1]], r_cart: numpy.ndarray[numpy.float64[m, 3]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def cartesian_to_east_north_up(r_ref: numpy.ndarray[numpy.float64[m, 3]], r_cart: numpy.ndarray[numpy.float64[m, 3]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def cartesian_to_geographical(r_cart: numpy.ndarray[numpy.float64[3, 1]], radius: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def cartesian_to_geographical(r_cart: numpy.ndarray[numpy.float64[3, 1]], radius: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def cartesian_to_geographical(r_cart: numpy.ndarray[numpy.float64[m, 3]], radius: float) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def cartesian_to_geographical(r_cart: numpy.ndarray[numpy.float64[m, 3]], radius: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def cartesian_to_spherical(r_cart: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def cartesian_to_spherical(r_cart: numpy.ndarray[numpy.float64[m, 3]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def classical_to_cartesian(coe: ClassicalOE, mu: float) -> CartesianOrbitState:
    ...
@typing.overload
def classical_to_cartesian(coe: numpy.ndarray[numpy.float64[6, 1]], mu: float) -> numpy.ndarray[numpy.float64[6, 1]]:
    ...
@typing.overload
def classical_to_cartesian(coe: numpy.ndarray[numpy.float64[6, 1]], mu: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_cartesian(coe: numpy.ndarray[numpy.float64[m, 6]], mu: float) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_cartesian(coe: numpy.ndarray[numpy.float64[m, 6]], mu: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_delaunay(arg0: ClassicalOE, arg1: float) -> ...:
    ...
@typing.overload
def classical_to_delaunay(coe: numpy.ndarray[numpy.float64[6, 1]], mu: float) -> numpy.ndarray[numpy.float64[6, 1]]:
    ...
@typing.overload
def classical_to_delaunay(coe: numpy.ndarray[numpy.float64[6, 1]], mu: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_delaunay(coe: numpy.ndarray[numpy.float64[m, 6]], mu: float) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_delaunay(coe: numpy.ndarray[numpy.float64[m, 6]], mu: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_equinoctial(arg0: ClassicalOE, arg1: float) -> EquinoctialOE:
    ...
@typing.overload
def classical_to_equinoctial(coe: numpy.ndarray[numpy.float64[6, 1]], mu: float) -> numpy.ndarray[numpy.float64[6, 1]]:
    ...
@typing.overload
def classical_to_equinoctial(coe: numpy.ndarray[numpy.float64[6, 1]], mu: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_equinoctial(coe: numpy.ndarray[numpy.float64[m, 6]], mu: float) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_equinoctial(coe: numpy.ndarray[numpy.float64[m, 6]], mu: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_quasi_nonsingular(coe: ClassicalOE, mu: float) -> QuasiNonsingularOE:
    ...
@typing.overload
def classical_to_quasi_nonsingular(coe: numpy.ndarray[numpy.float64[6, 1]], mu: float) -> numpy.ndarray[numpy.float64[6, 1]]:
    ...
@typing.overload
def classical_to_quasi_nonsingular(coe: numpy.ndarray[numpy.float64[6, 1]], mu: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_quasi_nonsingular(coe: numpy.ndarray[numpy.float64[m, 6]], mu: float) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def classical_to_quasi_nonsingular(coe: numpy.ndarray[numpy.float64[m, 6]], mu: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def compute_occultation(epoch: float, r1: numpy.ndarray[numpy.float64[3, 1]], r2: numpy.ndarray[numpy.float64[3, 1]], cs1: Frame, cs2: Frame, bodies: list[NaifId]) -> numpy.ndarray[numpy.float64[m, 1]]:
    """
    Compute occultation between two points
    """
@typing.overload
def compute_occultation(epoch: float, r1: numpy.ndarray[numpy.float64[m, 3]], r2: numpy.ndarray[numpy.float64[m, 3]], cs1: Frame, cs2: Frame, bodies: list[NaifId]) -> numpy.ndarray[numpy.float64[m, n]]:
    """
    Compute occultation between two points
    """
@typing.overload
def compute_occultation(epoch: numpy.ndarray[numpy.float64[m, 1]], r1: numpy.ndarray[numpy.float64[m, 3]], r2: numpy.ndarray[numpy.float64[m, 3]], cs1: Frame, cs2: Frame, bodies: list[NaifId]) -> numpy.ndarray[numpy.float64[m, n]]:
    """
    Compute occultation between two points
    """
@typing.overload
def convert_orbit_state(arg0: numpy.ndarray[numpy.float64[6, 1]], arg1: OrbitStateRepres, arg2: OrbitStateRepres, arg3: float) -> numpy.ndarray[numpy.float64[6, 1]]:
    ...
@typing.overload
def convert_orbit_state(arg0: numpy.ndarray[numpy.float64[6, 1]], arg1: numpy.ndarray[..., double>[6, 1]], arg2: OrbitStateRepres, arg3: OrbitStateRepres, arg4: OrbitStateRepres, arg5: float) -> numpy.ndarray[numpy.float64[6, 1]]:
    ...
@typing.overload
def dB2decimal(x: float) -> float:
    """
    Convert dB to decimal
    """
@typing.overload
def dB2decimal(x: numpy.ndarray[numpy.float64[m, n]]) -> numpy.ndarray[numpy.float64[m, n]]:
    ...
@typing.overload
def decimal2dB(x: float) -> float:
    """
    Convert decimal to dB
    """
@typing.overload
def decimal2dB(x: numpy.ndarray[numpy.float64[m, n]]) -> numpy.ndarray[numpy.float64[m, n]]:
    ...
def degrees2dms(deg: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
    Convert degrees to degrees, minutes, seconds
    """
@typing.overload
def delaunay_to_classical(deloe: ..., mu: float) -> ClassicalOE:
    ...
@typing.overload
def delaunay_to_classical(deloe: numpy.ndarray[numpy.float64[6, 1]], mu: float) -> numpy.ndarray[numpy.float64[6, 1]]:
    ...
@typing.overload
def delaunay_to_classical(deloe: numpy.ndarray[numpy.float64[6, 1]], mu: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def delaunay_to_classical(deloe: numpy.ndarray[numpy.float64[m, 6]], mu: float) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def delaunay_to_classical(deloe: numpy.ndarray[numpy.float64[m, 6]], mu: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
def dms2degrees(hms: numpy.ndarray[numpy.float64[3, 1]]) -> float:
    """
    Convert degrees, minutes, seconds to degrees
    """
@typing.overload
def east_north_up_to_cartesian(r_ref: numpy.ndarray[numpy.float64[3, 1]], r_enu: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def east_north_up_to_cartesian(r_ref: numpy.ndarray[numpy.float64[1, 3]], r_enu: numpy.ndarray[numpy.float64[1, 3]]) -> numpy.ndarray[numpy.float64[1, 3]]:
    ...
@typing.overload
def east_north_up_to_cartesian(r_ref: numpy.ndarray[numpy.float64[m, 3]], r_enu: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def east_north_up_to_cartesian(r_ref: numpy.ndarray[numpy.float64[3, 1]], r_enu: numpy.ndarray[numpy.float64[m, 3]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def east_north_up_to_cartesian(r_ref: numpy.ndarray[numpy.float64[m, 3]], r_enu: numpy.ndarray[numpy.float64[m, 3]]) -> numpy.ndarray[numpy.float64[m, 3]]:
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
def equinoctial_to_classical(eqoe: EquinoctialOE, mu: float) -> ClassicalOE:
    ...
@typing.overload
def equinoctial_to_classical(eqoe: numpy.ndarray[numpy.float64[6, 1]], mu: float) -> numpy.ndarray[numpy.float64[6, 1]]:
    ...
@typing.overload
def equinoctial_to_classical(eqoe: numpy.ndarray[numpy.float64[6, 1]], mu: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def equinoctial_to_classical(eqoe: numpy.ndarray[numpy.float64[m, 6]], mu: float) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def equinoctial_to_classical(eqoe: numpy.ndarray[numpy.float64[m, 6]], mu: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def geographical_to_cartesian(r_geo: numpy.ndarray[numpy.float64[3, 1]], radius: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def geographical_to_cartesian(r_geo: numpy.ndarray[numpy.float64[3, 1]], radius: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def geographical_to_cartesian(r_geo: numpy.ndarray[numpy.float64[m, 3]], radius: float) -> numpy.ndarray[numpy.float64[m, 3]]:
    ...
@typing.overload
def geographical_to_cartesian(r_geo: numpy.ndarray[numpy.float64[m, 3]], radius: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 3]]:
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
def quasi_nonsingular_to_classical(arg0: QuasiNonsingularOE, arg1: float) -> ClassicalOE:
    ...
@typing.overload
def quasi_nonsingular_to_classical(qnsoe: numpy.ndarray[numpy.float64[6, 1]], mu: float) -> numpy.ndarray[numpy.float64[6, 1]]:
    ...
@typing.overload
def quasi_nonsingular_to_classical(qnsoe: numpy.ndarray[numpy.float64[6, 1]], mu: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def quasi_nonsingular_to_classical(qnsoe: numpy.ndarray[numpy.float64[m, 6]], mu: float) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def quasi_nonsingular_to_classical(qnsoe: numpy.ndarray[numpy.float64[m, 6]], mu: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 6]]:
    ...
@typing.overload
def relative_quasi_nonsingular_to_classical(coe: ClassicalOE, rel_qnsoe: QuasiNonsingularROE) -> ClassicalOE:
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
def spherical_to_cartesian(r_sph: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def spherical_to_cartesian(r_sph: numpy.ndarray[numpy.float64[m, 3]]) -> numpy.ndarray[numpy.float64[m, 3]]:
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
@typing.overload
def wrapTo2Pi(arg0: float) -> float:
    """
    Wrap angle to [0, 2pi]
    """
@typing.overload
def wrapTo2Pi(arg0: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1]]:
    """
    Wrap angle to [0, 2pi]
    """
@typing.overload
def wrapToPi(arg0: float) -> float:
    """
    Wrap angle to [-pi, pi]
    """
@typing.overload
def wrapToPi(arg0: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1]]:
    """
    Wrap angle to [-pi, pi]
    """
A1MJD_OF_J2000: float = 21545.0
A1_TAI_OFFSET: float = 0.0343817
AU: float = 149597970.0
C: float = 299792.458
C22_MOON: float = 3.470983013194e-05
CARTESIAN: OrbitStateRepres  # value = <OrbitStateRepres.CARTESIAN: 0>
CLASSICAL_OE: OrbitStateRepres  # value = <OrbitStateRepres.CLASSICAL_OE: 1>
DAYS_PER_JULIAN_CENTURY: float = 36525.0
DAYS_PER_SEC: float = 1.1574074074074073e-05
DAYS_PER_YEAR: float = 365.25
DEG_PER_RAD: float = 57.29577951308232
DEIMOS: NaifId  # value = <NaifId.DEIMOS: 402>
DELAUNAY_OE: OrbitStateRepres  # value = <OrbitStateRepres.DELAUNAY_OE: 6>
E: float = 2.718281828459045
EARTH: NaifId  # value = <NaifId.EARTH: 399>
EARTH_BARYCENTER: NaifId  # value = <NaifId.EARTH_BARYCENTER: 3>
EARTH_MOON_BARYCENTER: NaifId  # value = <NaifId.EARTH_BARYCENTER: 3>
ECEF: Frame  # value = <Frame.ITRF: 0>
ECI: Frame  # value = <Frame.GCRF: 1>
EME: Frame  # value = <Frame.EME: 5>
EMR: Frame  # value = <Frame.EMR: 8>
EQUINOTICAL_OE: OrbitStateRepres  # value = <OrbitStateRepres.EQUINOTICAL_OE: 5>
GCRF: Frame  # value = <Frame.GCRF: 1>
GSE: Frame  # value = <Frame.GSE: 4>
ICRF: Frame  # value = <Frame.ICRF: 2>
ITRF: Frame  # value = <Frame.ITRF: 0>
J2_EARTH: float = 0.00108262668
J2_MOON: float = 9.09427845027e-05
JD_JAN_5_1941: float = 2430000.0
JD_MJD_OFFSET: float = 2400000.5
JD_NOV_17_1858: float = 2400000.5
JD_OF_J2000: float = 2451545.0
JUPITER: NaifId  # value = <NaifId.JUPITER: 599>
JUPITER_BARYCENTER: NaifId  # value = <NaifId.JUPITER_BARYCENTER: 5>
MARS: NaifId  # value = <NaifId.MARS: 499>
MARS_BARYCENTER: NaifId  # value = <NaifId.MARS_BARYCENTER: 4>
ME: Frame  # value = <Frame.ME: 11>
MERCURY: NaifId  # value = <NaifId.MERCURY: 199>
MERCURY_BARYCENTER: NaifId  # value = <NaifId.MERCURY_BARYCENTER: 1>
MI: Frame  # value = <Frame.MI: 9>
MJD_OF_J2000: float = 21545.0
MOD: Frame  # value = <Frame.MOD: 6>
MOON: NaifId  # value = <NaifId.MOON: 301>
MU_EARTH: float = 398600.4418
MU_MOON: float = 4902.800066
NEPTUNE_BARYCENTER: NaifId  # value = <NaifId.NEPTUNE_BARYCENTER: 8>
OMEGA_E_M: float = 2.6617e-06
OP: Frame  # value = <Frame.OP: 12>
PA: Frame  # value = <Frame.PA: 10>
PHOBOS: NaifId  # value = <NaifId.PHOBOS: 401>
PI: float = 3.141592653589793
PI_OVER_TWO: float = 1.5707963267948966
PLUTO_BARYCENTER: NaifId  # value = <NaifId.PLUTO_BARYCENTER: 9>
P_SUN: float = 4.529800412790905e-09
QUASINONSINGULAR_ROE: OrbitStateRepres  # value = <OrbitStateRepres.QUASINONSINGULAR_ROE: 9>
QUASI_NONSINGULAR_OE: OrbitStateRepres  # value = <OrbitStateRepres.QUASI_NONSINGULAR_OE: 2>
RAD_PER_DEG: float = 0.017453292519943295
R_EARTH: float = 6378.137
R_MOON: float = 1737.4
SATURN_BARYCENTER: NaifId  # value = <NaifId.SATURN_BARYCENTER: 6>
SECS_PER_DAY: float = 86400.0
SECS_PER_HOUR: float = 3600.0
SECS_PER_MINUTE: float = 60.0
SER: Frame  # value = <Frame.SER: 3>
SINGULAR_ROE: OrbitStateRepres  # value = <OrbitStateRepres.SINGULAR_ROE: 3>
SOLAR_SYSTEM_BARYCENTER: NaifId  # value = <NaifId.SOLAR_SYSTEM_BARYCENTER: 0>
SUN: NaifId  # value = <NaifId.SUN: 10>
S_AU: float = 0.0013579999999999998
TIME_OF_J2000: float = 883655990.85
TOD: Frame  # value = <Frame.TOD: 7>
TT_TAI_OFFSET: float = 32.184
TWO_PI: float = 6.283185307179586
URANUS_BARYCENTER: NaifId  # value = <NaifId.URANUS_BARYCENTER: 7>
VENUS: NaifId  # value = <NaifId.VENUS: 299>
VENUS_BARYCENTER: NaifId  # value = <NaifId.VENUS_BARYCENTER: 2>
d_E_EMB: float = 4671.0
d_E_M: float = 384400.0
d_M_EMB: float = 379729.0

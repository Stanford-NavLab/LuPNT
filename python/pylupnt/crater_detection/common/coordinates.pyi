from __future__ import annotations
from astropy.coordinates.funcs import cartesian_to_spherical
from astropy.coordinates.funcs import spherical_to_cartesian
import numpy as np
import numpy
from numpy import linalg as LA
from pylupnt.crater_detection.common import constants as const
from scipy.spatial.transform._rotation import Rotation
__all__ = ['ENU_system', 'LA', 'OrbitingBodyBase', 'Rotation', 'cartesian_to_spherical', 'const', 'nadir_attitude', 'np', 'spherical_to_cartesian', 'suborbital_coords']
class OrbitingBodyBase:
    """
    
        Base class implementing all positional and orientation attributes + methods.
        
    """
    q = ...
    quaternion = ...
    @classmethod
    def from_coordinates(cls, lat, long, height, attitude = None, Rbody = 1737.1, convert_to_radians = False):
        ...
    def __init__(self, position = None, attitude = None, primary_body_radius = 1737.1):
        ...
    def point_nadir(self):
        ...
    def rotate(self, axis: str, angle: float, degrees: bool = True, reset_first: bool = False):
        ...
    def set_coordinates(self, lat, long, height = None, point_nadir = False, convert_to_radians = True):
        ...
    def set_random_position(self, min_height = 150, max_height = 400, height = None):
        ...
    def suborbital_position(self):
        ...
    @property
    def T(self) -> numpy.ndarray:
        ...
    @T.setter
    def T(self, attitude: typing.Union[numpy.ndarray, scipy.spatial.transform._rotation.Rotation]):
        """
        
                Sets instance's attitude
        
                Parameters
                ----------
                attitude : np.ndarray, Rotation
                    Orientation / attitude matrix (3x3) or scipy.spatial.transform.Rotation
                
        """
    @property
    def attitude(self) -> numpy.ndarray:
        ...
    @attitude.setter
    def attitude(self, attitude: typing.Union[numpy.ndarray, scipy.spatial.transform._rotation.Rotation]):
        """
        
                Sets instance's attitude
        
                Parameters
                ----------
                attitude : np.ndarray, Rotation
                    Orientation / attitude matrix (3x3) or scipy.spatial.transform.Rotation
                
        """
    @property
    def coordinates(self):
        ...
    @property
    def height(self):
        ...
    @height.setter
    def height(self, height):
        """
        
                Adjusts radial height without changing angular position.
        
                Parameters
                ----------
                height: int, float
                    Height to set to in km.
                
        """
    @property
    def latitude(self):
        ...
    @property
    def longitude(self):
        ...
    @property
    def position(self) -> numpy.ndarray:
        ...
    @position.setter
    def position(self, position: numpy.ndarray):
        """
        
                Sets instance's position in Cartesian space.
        
                If set to None, a random position above the moon will be generated between 150 and 400 km height.
        
                Parameters
                ----------
                position : np.ndarray
                    3x1 position vector of camera.
                
        """
    @property
    def r(self) -> numpy.ndarray:
        ...
    @r.setter
    def r(self, position: numpy.ndarray):
        """
        
                Sets instance's position in Cartesian space.
        
                If set to None, a random position above the moon will be generated between 150 and 400 km height.
        
                Parameters
                ----------
                position : np.ndarray
                    3x1 position vector of camera.
                
        """
def ENU_system(r):
    """
    Return local East-North-Up (ENU) coordinate system for point defined by p.
    
        Using the coordinate system defined using:
    
        .. math::
            \\mathbf{u}_i = \\mathbf{p}^{(c)}_{M_i}/||\\mathbf{p}^{(c)}_{M_i}||
    
            \\mathbf{e}_i = cross(\\mathbf{k}, \\mathbf{u}_i )/|| cross(\\mathbf{k}, \\mathbf{u}_i) ||
    
            \\mathbf{n}_i = cross(\\mathbf{u}_i, \\mathbf{e}_i)/|| cross(\\mathbf{u}_i, \\mathbf{e}_i) ||
    
        with
    
        .. math::
            \\mathbf{k} = [0 & 0 & 1]^T
    
        and :math:`p_{Mi}` is the selenographic 3D cartesian coordinate derived from latitude & longitude.
    
        Parameters
        ----------
        r : np.ndarray
            (Nx)3x1 vector that defines origin.
    
        Returns
        -------
        e_i, n_i, u_i : np.ndarray
            Normalized i, j, k components of coordinate system.
    
        
    """
def nadir_attitude(r):
    """
    Return nadir-pointing (z-axis) coordinate system for point defined by r in world reference frame. X- and
        Y-components are defined by East and South respectively.
    
        Parameters
        ----------
        r : np.ndarray
            (Nx)3x1 vector that defines origin.
    
        Returns
        -------
        e_i, n_i, d_i : np.ndarray
            Normalized i, j, k components of coordinate system.
    
        
    """
def suborbital_coords(r, R_body = 1737.1):
    """
    Return coordinates directly below orbital position.
    
        Parameters
        ----------
        r : np.ndarray
            Position above body (e.g. Moon)
        R_body : np.ndarray
            Radius of body in km, defaults to const.RMOON
    
        Returns
        -------
        np.ndarray
            Suborbital coordinates
        
    """

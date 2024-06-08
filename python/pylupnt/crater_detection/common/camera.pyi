from __future__ import annotations
from astropy.coordinates.funcs import spherical_to_cartesian
from collections.abc import Iterable
import numpy
import numpy as np
from numpy import linalg as LA
from pylupnt.crater_detection.common import constants as const
import pylupnt.crater_detection.common.coordinates
from pylupnt.crater_detection.common.coordinates import OrbitingBodyBase
__all__ = ['Camera', 'Iterable', 'LA', 'OrbitingBodyBase', 'camera_matrix', 'const', 'np', 'projection_matrix', 'spherical_to_cartesian']
class Camera(pylupnt.crater_detection.common.coordinates.OrbitingBodyBase):
    """
    
        Camera data class with associated state attributes & functions.
        
    """
    @classmethod
    def from_coordinates(cls, lat, long, height, fov = 45, resolution = (256, 256), attitude = None, Rbody = 1737.1, convert_to_radians = False):
        ...
    def __init__(self, fov = 45, resolution = (256, 256), **kwargs):
        ...
    @property
    def K(self) -> numpy.ndarray:
        ...
    @property
    def P(self) -> numpy.ndarray:
        ...
    @property
    def camera_matrix(self) -> numpy.ndarray:
        ...
    @property
    def fov(self) -> typing.Tuple:
        ...
    @fov.setter
    def fov(self, fov):
        """
        
                Set instance's Field-of-View in radians.
        
                Parameters
                ----------
                fov: int, float, Iterable
                    Field-of-View angle (radians), if type is Iterable it will be interpreted as (fov_x, fov_y)
                
        """
    @property
    def projection_matrix(self) -> numpy.ndarray:
        ...
    @property
    def resolution(self) -> typing.Tuple:
        ...
    @resolution.setter
    def resolution(self, resolution):
        """
        
                Set instance's resolution in pixels.
        
                Parameters
                ----------
                resolution : int, Iterable
                    Image resolution, if type is Iterable it will be interpreted as (res_x, res_y)
                
        """
def camera_matrix(fov = 45, resolution = (256, 256), alpha = 0):
    """
    Returns camera matrix [1] from Field-of-View, skew, and offset.
    
        Parameters
        ----------
        fov : float, Iterable
            Field-of-View angle (degrees), if type is Iterable it will be interpreted as (fov_x, fov_y)
        resolution : float, Iterable
            X- and Y-resolution of the image in pixels
        alpha : float
            Camera skew angle.
    
        Returns
        -------
        np.ndarray
            3x3 camera matrix
    
        References
        ----------
        .. [1] https://www.cs.ucf.edu/~mtappen/cap5415/lecs/lec19.pdf
        
    """
def projection_matrix(K, T_CM, r_M):
    """
    Return Projection matrix [1] according to:
    
        .. math:: \\mathbf{P}_C = \\mathbf{K} [ \\mathbf{T^C_M} & -r_C]
    
        Parameters
        ----------
        K : np.ndarray
            3x3 camera matrix
        T_CM : np.ndarray
            3x3 attitude matrix of camera in selenographic frame.
        r_M : np.ndarray
            3x1 camera position in world reference frame
    
        Returns
        -------
        np.ndarray
            3x4 projection matrix
    
        References
        ----------
        .. [1] https://www.cs.ucf.edu/~mtappen/cap5415/lecs/lec19.pdf
    
        See Also
        --------
        camera_matrix
    
        
    """

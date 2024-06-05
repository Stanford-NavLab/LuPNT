from __future__ import annotations
from astropy.coordinates.funcs import spherical_to_cartesian
import cv2 as cv2
from matplotlib.collections import EllipseCollection
from matplotlib import pyplot as plt
from numba.core.decorators import njit
import numpy
import numpy as np
from numpy import linalg as LA
import pylupnt.crater_detection.common.camera
from pylupnt.crater_detection.common.camera import Camera
from pylupnt.crater_detection.common.camera import camera_matrix
from pylupnt.crater_detection.common.camera import projection_matrix
from pylupnt.crater_detection.common import constants as const
from pylupnt.crater_detection.common.coordinates import ENU_system
from pylupnt.crater_detection.common.robbins import extract_robbins_dataset
from pylupnt.crater_detection.common.robbins import load_craters
from scipy.spatial.distance import cdist
import torch as torch
__all__ = ['Camera', 'ConicProjector', 'ENU_system', 'EllipseCollection', 'LA', 'MaskGenerator', 'camera_matrix', 'cdist', 'conic_center', 'conic_center_numba', 'conic_matrix', 'const', 'crater_camera_homography', 'cv2', 'ellipse_angle', 'ellipse_axes', 'extract_robbins_dataset', 'generate_mask', 'load_craters', 'matrix_adjugate', 'njit', 'np', 'plot_conics', 'plt', 'project_crater_centers', 'project_crater_conics', 'projection_matrix', 'scale_det', 'spherical_to_cartesian', 'torch']
class ConicProjector(pylupnt.crater_detection.common.camera.Camera):
    def generate_mask(self, A_craters = None, C_craters = None, r_craters = None, **kwargs):
        ...
    def plot(self, A_craters = None, C_craters = None, r_craters = None, **kwargs):
        ...
    def project_crater_centers(self, r_craters):
        ...
    def project_crater_conics(self, C_craters, r_craters):
        ...
class MaskGenerator(ConicProjector):
    @classmethod
    def from_robbins_dataset(cls, file_path = 'data/lunar_crater_database_robbins_2018.csv', diamlims = (4, 100), ellipse_limit = 1.3, arc_lims = 0.0, axis_threshold = (5, 100), filled = False, instancing = True, mask_thickness = 1, position = None, resolution = (256, 256), fov = 45, primary_body_radius = 1737.1, **load_crater_kwargs):
        ...
    def __init__(self, r_craters_catalogue: numpy.ndarray, C_craters_catalogue: numpy.ndarray, axis_threshold = (5, 100), filled = False, instancing = True, mask_thickness = 1, mask_margin = 0, **kwargs):
        ...
    def _visible(self):
        ...
    def craters_in_image(self, margin = None):
        ...
    def generate_mask(self, **kwargs):
        ...
    def plot(self, *args, **kwargs):
        ...
    def visible_catalogue_craters(self, margin = None):
        ...
def conic_center(A):
    ...
def conic_center_numba(*args, **kwargs):
    ...
def conic_matrix(a, b, psi, x = 0, y = 0):
    """
    Returns matrix representation for crater derived from ellipse parameters
    
        Parameters
        ----------
        a: np.ndarray, torch.Tensor, int, float
            Semi-major ellipse axis
        b: np.ndarray, torch.Tensor, int, float
            Semi-minor ellipse axis
        psi: np.ndarray, torch.Tensor, int, float
            Ellipse angle (radians)
        x: np.ndarray, torch.Tensor, int, float
            X-position in 2D cartesian coordinate system (coplanar)
        y: np.ndarray, torch.Tensor, int, float
            Y-position in 2D cartesian coordinate system (coplanar)
    
        Returns
        -------
        np.ndarray, torch.Tensor
            Array of ellipse matrices
        
    """
def crater_camera_homography(r_craters, P_MC):
    """
    Calculate homography between crater-plane and camera reference frame.
    
        .. math:: \\mathbf{H}_{C_i} =  ^\\mathcal{M}\\mathbf{P}_\\mathcal{C_{craters}} [[H_{M_i}], [k^T]]
    
        Parameters
        ----------
        r_craters : np.ndarray
            (Nx)3x1 position vector of craters.
        P_MC : np.ndarray
            (Nx)3x4 projection matrix from selenographic frame to camera pixel frame.
    
        Returns
        -------
            (Nx)3x3 homography matrix
        
    """
def ellipse_angle(A):
    ...
def ellipse_axes(A):
    ...
def generate_mask(A_craters, resolution = (256, 256), filled = False, instancing = False, thickness = 1):
    ...
def matrix_adjugate(matrix):
    """
    Return adjugate matrix [1].
    
        Parameters
        ----------
        matrix : np.ndarray
            Input matrix
    
        Returns
        -------
        np.ndarray
            Adjugate of input matrix
    
        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Adjugate_matrix
        
    """
def plot_conics(A_craters: typing.Union[numpy.ndarray, torch.Tensor], resolution = (256, 256), figsize = (15, 15), plot_centers = False, ax = None, rim_color = 'r', alpha = 1.0):
    ...
def project_crater_centers(r_craters, fov, resolution, T_CM, r_M):
    """
    Project crater centers into digital pixel frame.
    
        Parameters
        ----------
        r_craters : np.ndarray
            Nx3x1 position vector of craters.
        fov : int, float, Iterable
            Field-of-View angle (radians), if type is Iterable it will be interpreted as (fov_x, fov_y)
        resolution : int, Iterable
            Image resolution, if type is Iterable it will be interpreted as (res_x, res_y)
        T_CM : np.ndarray
            3x3 matrix representing camera attitude in world reference frame
        r_M : np.ndarray
            3x1 position vector of camera
    
        Returns
        -------
        np.ndarray
            Nx2x1 2D positions of craters in pixel frame
        
    """
def project_crater_conics(C_craters, r_craters, fov, resolution, T_CM, r_M):
    """
    Project crater conics into digital pixel frame. See pages 17 - 25 from [1] for methodology.
    
        Parameters
        ----------
        C_craters : np.ndarray
            Nx3x3 array of crater conics
        r_craters : np.ndarray
            Nx3x1 position vector of craters.
        fov : float, Iterable
            Field-of-View angle (radians), if type is Iterable it will be interpreted as (fov_x, fov_y)
        resolution : int, Iterable
            Image resolution, if type is Iterable it will be interpreted as (res_x, res_y)
        T_CM : np.ndarray
            3x3 matrix representing camera attitude in world reference frame
        r_M : np.ndarray
            3x1 position vector of camera
    
        Returns
        -------
        np.ndarray
            Nx3x3 Homography matrix H_Ci
    
        References
        ----------
        .. [1] Christian, J. A., Derksen, H., & Watkins, R. (2020). Lunar Crater Identification in Digital Images. https://arxiv.org/abs/2009.01228
        
    """
def scale_det(matrix):
    """
    Rescale matrix such that det(A) = 1.
    
        Parameters
        ----------
        matrix: np.ndarray, torch.Tensor
            Matrix input
        Returns
        -------
        np.ndarray
            Normalised matrix.
        
    """

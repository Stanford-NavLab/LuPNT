from __future__ import annotations
from numba.core.decorators import njit
import numpy as np
import numpy
from numpy import linalg as LA
from pylupnt.crater_detection.common.conics import ConicProjector
from pylupnt.crater_detection.common.conics import ellipse_axes
from pylupnt.crater_detection.common import constants as const
from pylupnt.crater_detection.common.coordinates import ENU_system
from pylupnt.crater_detection.detection.metrics import gaussian_angle_distance
from scipy.optimize._lsq.least_squares import least_squares
from sklearn.linear_model._ransac import RANSACRegressor
import warnings as warnings
__all__ = ['ConicProjector', 'ENU_system', 'HJacobian', 'Hx', 'LA', 'PositionRegressor', 'RANSACRegressor', 'const', 'ellipse_axes', 'gaussian_angle_distance', 'least_squares', 'njit', 'np', 'systems_dynamics_matrix', 'warnings']
class PositionRegressor:
    def __init__(self, sigma_pix = 6, **ransac_kwargs):
        ...
    def fit(self, A_query, C_query, r_query, attitude, camera_matrix, reprojection = True):
        ...
    def ransac_match(self):
        ...
    def reprojection_match(self):
        ...
    @property
    def num_inliers(self):
        ...
    @property
    def num_verified(self):
        ...
def HJacobian(*args) -> numpy.ndarray:
    ...
def Hx(x_state: numpy.ndarray) -> numpy.ndarray:
    ...
def _model_validator(min_alt = 30, max_alt = 500, primary_body_radius = 1737.1):
    ...
def systems_dynamics_matrix(*args, **kwargs):
    ...

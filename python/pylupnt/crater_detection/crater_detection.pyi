from PIL import Image
from __future__ import annotations
import cv2 as cv2
import numpy
import numpy as np
import os as os
import pylupnt as pnt
from pylupnt.crater_detection.common.conics import conic_center
from pylupnt.crater_detection.common.conics import conic_matrix
from pylupnt.crater_detection.common.conics import ellipse_angle
from pylupnt.crater_detection.common.conics import ellipse_axes
from pylupnt.crater_detection.common.conics import plot_conics
from pylupnt.crater_detection.common.conics import project_crater_centers
from pylupnt.crater_detection.common.conics import scale_det
from pylupnt.crater_detection.common.coordinates import ENU_system
from pylupnt.crater_detection.common.coordinates import nadir_attitude
from pylupnt.crater_detection.common.robbins import extract_robbins_dataset
from pylupnt.crater_detection.common.robbins import load_craters
from scipy.ndimage._morphology import binary_dilation
__all__ = ['ENU_system', 'Image', 'a_crater', 'b_crater', 'binary_dilation', 'conic_center', 'conic_matrix', 'crater_detection', 'cv2', 'ellipse_angle', 'ellipse_axes', 'extract_robbins_dataset', 'horizon_detection', 'horizon_matching', 'id_crater', 'lat_crater', 'load_craters', 'lon_crater', 'nadir_attitude', 'np', 'os', 'plot_conics', 'pnt', 'project_crater_centers', 'psi_crater', 'scale_det', 'trace_lines_and_record']
def crater_detection(r_cam_pa, R_pa2cam, r_sun_pa, fov_cam, res_cam, seed = None):
    ...
def horizon_detection(img, R_pa2cam, r_sun_pa, K, N_lines = 200, threshold = 110):
    ...
def horizon_matching(u: numpy.ndarray, K_inv: numpy.ndarray, R_cam2pa: numpy.ndarray, a: float, b: float = None, c: float = None):
    """
    
        Christian-Robinson OPNAV algorithm
    
        Args:
            u (np.ndarray):
        
    """
def trace_lines_and_record(img, u_illum, N_lines, threshold, return_lines = False):
    ...
__warningregistry__: dict = {'version': 194}
a_crater: numpy.ndarray  # value = array([289.44  , 233.763 , 234.321 , ...,  72.8571,  30.0635,  72.1062])
b_crater: numpy.ndarray  # value = array([245.786 , 225.266 , 212.619 , ...,  69.2105,  25.1541,  61.5883])
id_crater: numpy.ndarray  # value = array(['00-1-000001', '00-1-000020', '00-1-000023', ..., '10-2-014982',...
lat_crater: numpy.ndarray  # value = array([ 0.77507105,  0.79423127, -0.39516   , ..., -1.41853729,...
lon_crater: numpy.ndarray  # value = array([5.74293609, 2.66817209, 6.17034977, ..., 0.64088316, 5.34720014,...
psi_crater: numpy.ndarray  # value = array([2.21662051, 2.96021804, 2.57024167, ..., 1.07212609, 3.0668402 ,...

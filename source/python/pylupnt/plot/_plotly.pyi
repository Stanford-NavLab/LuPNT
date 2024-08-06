from PIL import Image
import PIL.JpegImagePlugin
from __future__ import annotations
import numpy as np
import numpy
import os as os
from plotly import express as px
from plotly import graph_objs as go
import plotly.graph_objs._figure
import plotly.graph_objs._mesh3d
from plotly import io as pio
import pylupnt as pnt
from pylupnt import utils
from sklearn.cluster._kmeans import KMeans
from sklearn.utils._indexing import shuffle
__all__ = ['EARTH_SURFACE', 'IMAGES', 'Image', 'KMeans', 'MOON_SURFACE', 'PLOTLY_COLORS', 'RADII', 'axis_dict', 'create_sphere_meshgrid', 'get_body_trace', 'go', 'image2zvals', 'mesh_data', 'np', 'os', 'pio', 'plot_3d_arrow', 'plot_constellation', 'plot_frame', 'pnt', 'px', 'regular_tri', 'set_equal_aspect_ratio', 'set_view', 'shuffle', 'utils']
def create_sphere_meshgrid(rows: int, cols: int, radius: float = 1.0) -> tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """
    
        Create a sphere meshgrid
    
        Args:
            rows (int): number of rows
            cols (int): number of columns
        Returns:
            tuple[np.ndarray, np.ndarray, np.ndarray]: x, y, z coordinates of the sphere
        
    """
def get_body_trace(body = ..., size_factor: int = 5, R_b2frame: numpy.ndarray = None, r_b2s_pa: numpy.ndarray = None, r_body: numpy.ndarray = None, alpha: float = 0.2) -> plotly.graph_objs._mesh3d.Mesh3d:
    """
    
        Get the moon trace
    
        Args:
            body (int): celestial body
            size_factor (int): size factor
            R_b2frame (np.ndarray): rotation matrix from body to frame
            r_b2s_pa (np.ndarray): vector from body to sun in the frame
            alpha (float): light intensity
        
    """
def image2zvals(img: numpy.ndarray, n_colors: int = 64, n_training_pixels: int = 800, rngs: int = 123) -> tuple[numpy.ndarray, list]:
    """
    
        Image color quantization
    
        Args:
            img (np.ndarray): image array
            n_colors (int): number of colors for color quantization
            n_training_pixels (int): number of image pixels to fit a KMeans instance to them
            rngs (int): random seed
        Returns:
            tuple[np.ndarray, list]: z_values for the heatmap representation, and a plotly colorscale
        
    """
def mesh_data(img, n_colors = 32, n_training_pixels = 800):
    ...
def plot_3d_arrow(fig: plotly.graph_objs._figure.Figure, origin: numpy.ndarray, direction: numpy.ndarray, length_scale: float = 1, tip_scale: float = 1, color: typing.Union[str, list[str]] = 'black', width: float = 2) -> None:
    """
    
        Add 3D arrow to a plotly figure
    
        Args:
            fig (go.Figure): plotly figure
            origin (np.ndarray): origin of the arrow
            direction (np.ndarray): direction of the arrow
            length_scale (float): scale factor for the length of the arrow
            tip_scale (float): scale factor for the tip of the arrow
            color (str): color of the arrow
            width (float): width of the arrow
        
    """
def plot_constellation(rv: numpy.ndarray, t: int = 0, marker_size = 4, fig = None, **kwargs) -> plotly.graph_objs._figure.Figure:
    ...
def plot_frame(fig: plotly.graph_objs._figure.Figure, origin: numpy.ndarray, rotation: numpy.ndarray, length_scale: float = 1, tip_scale: float = 1, width: float = 2) -> None:
    ...
def regular_tri(rows, cols):
    ...
def set_equal_aspect_ratio(fig: plotly.graph_objs._figure.Figure):
    ...
def set_view(fig: plotly.graph_objs._figure.Figure, azimuth: float, elevation: float, zoom: float = 1.0):
    ...
EARTH_SURFACE: PIL.JpegImagePlugin.JpegImageFile  # value = <PIL.JpegImagePlugin.JpegImageFile image mode=RGB size=1024x512 at 0x141577490>
IMAGES: dict  # value = {<NaifId.MOON: 301>: array([[[183, 180, 173],...
MOON_SURFACE: PIL.JpegImagePlugin.JpegImageFile  # value = <PIL.JpegImagePlugin.JpegImageFile image mode=RGB size=1024x512 at 0x142E68210>
PLOTLY_COLORS: list = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD', '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF']
RADII: dict  # value = {<NaifId.MOON: 301>: 1737.4, <NaifId.EARTH: 399>: 6378.137}
axis_dict: dict = {'mirror': True, 'ticks': 'outside', 'showline': True, 'showgrid': True, 'automargin': True}

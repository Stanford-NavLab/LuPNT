from PIL import Image
import PIL.JpegImagePlugin
from __future__ import annotations
from matplotlib import colors as mcolors
from matplotlib import pyplot as plt
import numpy as np
import os as os
from plotly import express as px
from plotly import graph_objs as go
from plotly import io as pio
import pylupnt as pnt
from pylupnt.plot._mpl import Plot3D
from pylupnt.plot._plotly import create_sphere_meshgrid
from pylupnt.plot._plotly import get_body_trace
from pylupnt.plot._plotly import image2zvals
from pylupnt.plot._plotly import mesh_data
from pylupnt.plot._plotly import plot_3d_arrow
from pylupnt.plot._plotly import plot_constellation
from pylupnt.plot._plotly import plot_frame
from pylupnt.plot._plotly import regular_tri
from pylupnt.plot._plotly import set_equal_aspect_ratio
from pylupnt.plot._plotly import set_view
from pylupnt import utils
from sklearn.cluster._kmeans import KMeans
from sklearn.utils._indexing import shuffle
from . import _mpl
from . import _plotly
__all__ = ['COLORS', 'EARTH_SURFACE', 'IMAGES', 'Image', 'KMeans', 'MOON_SURFACE', 'PLOTLY_COLORS', 'Plot3D', 'RADII', 'axis_dict', 'create_sphere_meshgrid', 'get_body_trace', 'go', 'image2zvals', 'mcolors', 'mesh_data', 'np', 'os', 'pio', 'plot_3d_arrow', 'plot_constellation', 'plot_data', 'plot_frame', 'plt', 'pnt', 'px', 'regular_tri', 'set_equal_aspect_ratio', 'set_view', 'shuffle', 'utils']
COLORS: list = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
EARTH_SURFACE: PIL.JpegImagePlugin.JpegImageFile  # value = <PIL.JpegImagePlugin.JpegImageFile image mode=RGB size=1024x512 at 0x16B4F5D90>
IMAGES: dict  # value = {<NaifId.MOON: 301>: array([[[183, 180, 173],...
MOON_SURFACE: PIL.JpegImagePlugin.JpegImageFile  # value = <PIL.JpegImagePlugin.JpegImageFile image mode=RGB size=1024x512 at 0x16AC088D0>
PLOTLY_COLORS: list = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD', '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF']
RADII: dict  # value = {<NaifId.MOON: 301>: 1737.4, <NaifId.EARTH: 399>: 6378.137}
axis_dict: dict = {'mirror': True, 'ticks': 'outside', 'showline': True, 'showgrid': True, 'automargin': True}
plot_data: dict  # value = {<NaifId.EARTH: 399>: {'filename': 'earth_surface.jpg', 'RE': 6378.137, 'lim': 25000.0, 'brightness': 3}, <NaifId.MOON: 301>: {'filename': 'moon_surface.jpeg', 'RE': 1737.1, 'lim': 10000.0, 'brightness': 1.5}}

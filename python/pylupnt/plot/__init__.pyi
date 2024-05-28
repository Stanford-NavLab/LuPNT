from PIL import Image
import PIL.JpegImagePlugin
from __future__ import annotations
from matplotlib import colors as mcolors
from matplotlib import pyplot as plt
import numpy as np
import numpy
import os as os
from plotly import express as px
from plotly import graph_objs as go
from plotly import io as pio
from pylupnt.plot._mpl import Plot3D
from pylupnt.plot._plotly import create_sphere_meshgrid
from pylupnt.plot._plotly import get_moon_trace
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
__all__ = ['COLORS', 'Image', 'KMeans', 'MOON_SURFACE', 'Plot3D', 'axis_dict', 'create_sphere_meshgrid', 'get_moon_trace', 'go', 'image2zvals', 'img', 'mcolors', 'mesh_data', 'np', 'os', 'pio', 'plot_3d_arrow', 'plot_constellation', 'plot_data', 'plot_frame', 'plt', 'px', 'regular_tri', 'set_equal_aspect_ratio', 'set_view', 'shuffle', 'utils']
COLORS: list = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
MOON_SURFACE: PIL.JpegImagePlugin.JpegImageFile  # value = <PIL.JpegImagePlugin.JpegImageFile image mode=RGB size=1024x512 at 0x14E6DC050>
axis_dict: dict = {'mirror': True, 'ticks': 'outside', 'showline': True, 'showgrid': True, 'automargin': True}
img: numpy.ndarray  # value = array([[[183, 180, 173],...
plot_data: dict  # value = {<NaifId.EARTH: 399>: {'filename': 'earth_surface.jpg', 'RE': 6378.137, 'lim': 25000.0, 'brightness': 3}, <NaifId.MOON: 301>: {'filename': 'moon_surface.jpeg', 'RE': 1737.1, 'lim': 10000.0, 'brightness': 1.5}}

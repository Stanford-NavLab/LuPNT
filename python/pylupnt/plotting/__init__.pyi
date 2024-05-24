from PIL import Image
from __future__ import annotations
from matplotlib import colors as mcolors
from matplotlib import pyplot as plt
import numpy as np
import numpy
import os as os
from plotly import express as px
from plotly import graph_objs as go
from plotly import io as pio
from pylupnt.plotting._mpl import Plot3D
from pylupnt.plotting._plotly import add_3d_arrow
from pylupnt.plotting._plotly import create_sphere_meshgrid
from pylupnt.plotting._plotly import get_moon_trace
from pylupnt.plotting._plotly import image2zvals
from pylupnt.plotting._plotly import mesh_data
from pylupnt.plotting._plotly import plot_constellation
from pylupnt.plotting._plotly import regular_tri
from pylupnt.plotting._plotly import set_equal_aspect_ratio
from pylupnt.plotting._plotly import set_view
from pylupnt import utils
from sklearn.cluster._kmeans import KMeans
from sklearn.utils._indexing import shuffle
from . import _mpl
from . import _plotly
__all__ = ['COLORS', 'Image', 'KMeans', 'Plot3D', 'RE', 'add_3d_arrow', 'axis_dict', 'create_sphere_meshgrid', 'get_moon_trace', 'go', 'image2zvals', 'img', 'mcolors', 'mesh_data', 'np', 'os', 'pio', 'plot_constellation', 'plot_data', 'plt', 'px', 'regular_tri', 'set_equal_aspect_ratio', 'set_view', 'shuffle', 'utils']
COLORS: list = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
RE: float = 1737.1
axis_dict: dict = {'mirror': True, 'ticks': 'outside', 'showline': True, 'showgrid': True, 'automargin': True}
img: numpy.ndarray  # value = array([[[183, 180, 173],...
plot_data: dict  # value = {<NaifId.EARTH: 399>: {'filename': 'earth_surface.jpg', 'RE': 6378.137, 'lim': 25000.0, 'brightness': 3}, <NaifId.MOON: 301>: {'filename': 'moon_surface.jpeg', 'RE': 1737.1, 'lim': 10000.0, 'brightness': 1.5}}

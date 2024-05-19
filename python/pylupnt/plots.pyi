from PIL import Image
from __future__ import annotations
from matplotlib import colors as mcolors
from matplotlib import pyplot as plt
import numpy as np
import os as os
from pylupnt import pylupnt_pybind as _pnt
from pylupnt import utils as u
__all__ = ['COLORS', 'Image', 'Plot3D', 'mcolors', 'np', 'os', 'plot_data', 'plt', 'u']
class Plot3D:
    def __init__(self, azim = -60, elev = 30, figsize = (10, 10)):
        ...
    def check_occultation(self, data):
        ...
    def plot(self, data, *args, mask = False, **kwargs):
        """
        
                Plot Cartesian coordinates
                
        """
    def plot_surface(self, name, offset = ..., adjust_axis = True, limit = None, scale = 3):
        ...
    def plot_visible(self, azimuth, elev):
        ...
    def rotate(self, event):
        ...
    def scatter(self, data, mask = False, *args, **kwargs):
        """
        
                Plot Cartesian coordinates
                
        """
    def set_labelpad(self, padx: int, pady: int, padz: int) -> None:
        ...
    def set_labels(self, x: str, y: str, z: str) -> None:
        ...
    def set_lims(self, xlims: tuple, ylims: tuple, zlims: tuple, equal = True) -> None:
        ...
    def set_pane_color(self, color: tuple) -> None:
        ...
    def set_tick_multiplier(self, factor: int) -> None:
        ...
    def set_tickpad(self, pad: int) -> None:
        ...
    def set_ticks(self, x: list, y: list, z: list) -> None:
        ...
COLORS: list = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
plot_data: dict  # value = {<NaifId.EARTH: 399>: {'filename': 'earth_surface.jpg', 'RE': 6378.137, 'lim': 25000.0, 'brightness': 3}, <NaifId.MOON: 301>: {'filename': 'moon_surface.jpeg', 'RE': 1737.1, 'lim': 10000.0, 'brightness': 1.5}}

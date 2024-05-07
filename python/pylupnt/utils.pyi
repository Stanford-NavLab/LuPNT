from __future__ import annotations
from matplotlib.gridspec import GridSpec
from matplotlib import pyplot as plt
import numpy as np
import os as os
import pandas as pd
import pickle as pickle
from time import time
__all__ = ['GridSpec', 'arr_to_mat_idx', 'basepath', 'dump_pickle', 'format_element', 'i_to_arr_idxs', 'load_data', 'load_pickle', 'mat_to_arr_idx', 'np', 'os', 'pd', 'pickle', 'plot_RTN', 'plt', 'print_aligned', 'set_axes_equal', 'time', 'timed', 'timer_func', 'wrapToPi']
def arr_to_mat_idx(k, n):
    ...
def dump_pickle(obj, path):
    ...
def format_element(x, fmt = '{}'):
    ...
def i_to_arr_idxs(i, N):
    ...
def load_data(directory):
    ...
def load_pickle(path):
    ...
def mat_to_arr_idx(i, j, n):
    ...
def plot_RTN(rv_RTN, labels = None, legend_text = None, init = False, final = True, center = True):
    ...
def print_aligned(matrix):
    ...
def set_axes_equal(ax):
    """
    Make axes of 3D plot have equal scale so that spheres appear as spheres,
        cubes as cubes, etc..  This is one possible solution to Matplotlib's
        ax.set_aspect('equal') and ax.axis('equal') not working for 3D.
    
        Input
          ax: a matplotlib axis, e.g., as output from plt.gca().
        
    """
def timed(func, *args, **kwargs):
    ...
def timer_func(func):
    ...
def wrapToPi(x):
    """
    
        Wrap angle in radians to [-pi pi]
    
        Args:
            x (float): angle in radians
        Returns:
            x (float): angle in radians wrapped to [-pi pi]
        
    """
basepath: str = '/Users/guillemcv/Development/NavLab/LuPNT/data'

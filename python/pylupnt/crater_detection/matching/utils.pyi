from __future__ import annotations
import networkx as nx
from numba.core.decorators import njit
import numpy as np
import numpy
from numpy import linalg as LA
__all__ = ['LA', 'all_clockwise', 'cw_or_ccw', 'cyclic_permutations', 'enhanced_pattern_shifting', 'eps_array', 'get_cliques_by_length', 'is_clockwise', 'is_colinear', 'latlong2cartesian', 'njit', 'np', 'np_swap_columns', 'nx', 'shift_nd', 'triad_splice']
def all_clockwise(x_triads_, y_triads_):
    ...
def cw_or_ccw(x_triads, y_triads):
    ...
def cyclic_permutations(it, step = 1):
    """
    Returns cyclic permutations for iterable.
    
        Parameters
        ----------
        it : iterable object
        step : int, optional
    
        Yields
        -------
        Cyclic permutation of it
        
    """
def enhanced_pattern_shifting(*args, **kwargs):
    """
    Generator function returning next crater triad according to Enhanced Pattern Shifting Method [1].
    
        Parameters
        ----------
        n : int
            Number of detected instances.
        start_n: int
            Iteration to start from, useful for batch processing of triads.
    
        Returns
        -------
        i, j, k : int
    
        References
        ----------
        .. [1] Arnas, D., Fialho, M. A. A., & Mortari, D. (2017). Fast and robust kernel generators for star trackers. Acta Astronautica, 134 (August 2016), 291–302. https://doi.org/10.1016/j.actaastro.2017.02.016
        
    """
def eps_array(*args, **kwargs):
    ...
def get_cliques_by_length(G, length_clique):
    """
    Return the list of all cliques in an undirected graph G with length
        equal to length_clique.
    """
def is_clockwise(x_triads, y_triads):
    """
    Returns boolean array which tells whether the three points in 2D plane given by x_triads & y_triads are
        oriented clockwise. https://en.wikipedia.org/wiki/Curve_orientation#Orientation_of_a_simple_polygon
    
        Parameters
        ----------
        x_triads, y_triads : np.ndarray
            Array of 2D coordinates for triangles in a plane
    
        Returns
        -------
        np.ndarray
        
    """
def is_colinear(x_triads_, y_triads_):
    ...
def latlong2cartesian(lat, long, alt = 0, rad = 1737.1):
    """
    
        Calculate Cartesian coordinates from latitude + longitude information
        
    """
def np_swap_columns(arr):
    ...
def shift_nd(arr: numpy.ndarray, shift: numpy.ndarray):
    ...
def triad_splice(arr, triads):
    ...

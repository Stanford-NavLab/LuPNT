from __future__ import annotations
import numpy
import numpy as np
__all__ = ['arr_to_mat_idx', 'cross_norm', 'i_to_arr_idxs', 'mat_to_arr_idx', 'np', 'wrapToPi']
def arr_to_mat_idx(k: int, n: int) -> tuple:
    """
    
        Convert array index to matrix indices
    
        Args:
            k (int): array index
            n (int): number of rows
        Returns:
            i (int): row index
            j (int): column index
        
    """
def cross_norm(a: numpy.ndarray, b: numpy.ndarray) -> numpy.ndarray:
    """
    
        Compute the norm of the cross product of two vectors
    
        Args:
            a (np.ndarray): first vector
            b (np.ndarray): second vector
        Returns:
            np.ndarray: norm of the cross product
        
    """
def i_to_arr_idxs(i: int, N: int) -> list:
    """
    
        Convert matrix index to array indices
    
        Args:
            i (int): row index
            N (int): number of rows
        Returns:
            arr1 (list): array indices for upper triangular part of the matrix
            arr2 (list): array indices for lower triangular part of the matrix
        
    """
def mat_to_arr_idx(i: int, j: int, n: int) -> int:
    """
    
        Convert matrix indices to array index
    
        Args:
            i (int): row index
            j (int): column index
            n (int): number of rows
        Returns:
            k (int): array index
        
    """
def wrapToPi(x: typing.Union[float, numpy.ndarray]) -> typing.Union[float, numpy.ndarray]:
    """
    
        Wrap angle in radians to [-pi pi]
    
        Args:
            x (float): angle in radians
        Returns:
            x (float): angle in radians wrapped to [-pi pi]
        
    """

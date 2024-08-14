import numpy as np
from typing import Union


def wrap2Pi(x: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Wrap angle in radians to [-pi pi]

    Args:
        x (float): angle in radians
    Returns:
        x (float): angle in radians wrapped to [-pi pi]
    """
    x = np.mod(x + np.pi, 2 * np.pi) - np.pi
    return x


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
    assert i != j, "Diagonal elements are not stored"
    assert j > i, "This function is only valid for j > i"
    k = (n * (n - 1) / 2) - (n - i) * ((n - i) - 1) / 2 + j - i - 1
    return int(k)


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
    i = n - 2 - np.floor(np.sqrt(-8 * k + 4 * n * (n - 1) - 7) / 2.0 - 0.5)
    j = k + i + 1 - n * (n - 1) / 2 + (n - i) * ((n - i) - 1) / 2
    return int(i), int(j)


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
    arr1 = [mat_to_arr_idx(i, j, N) for j in range(i + 1, N)]
    arr2 = [mat_to_arr_idx(j, i, N) for j in range(0, i)]
    return sorted(arr1 + arr2)


def cross_norm(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Compute the norm of the cross product of two vectors

    Args:
        a (np.ndarray): first vector
        b (np.ndarray): second vector
    Returns:
        np.ndarray: norm of the cross product
    """
    cross = np.cross(a, b, axis=-1)
    return cross / np.linalg.norm(cross, axis=-1)[..., np.newaxis]

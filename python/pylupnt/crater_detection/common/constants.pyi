from __future__ import annotations
import pathlib
from pathlib import Path
__all__ = ['ARC_LIMS', 'AXIS_THRESHOLD', 'CAMERA_FOV', 'CAMERA_RESOLUTION', 'DB_CAM_ALTITUDE', 'DIAMLIMS', 'FILLED', 'GENERATION_KWARGS', 'INSTANCING', 'KERNEL_ROOT', 'MASK_THICKNESS', 'MAX_ELLIPTICITY', 'MAX_SOL_INCIDENCE', 'MIN_SOL_INCIDENCE', 'Path', 'RANDOMIZED_ORIENTATION', 'RMOON', 'SAVE_CRATERS', 'SPICE_BASE_URL', 'TRIAD_RADIUS']
ARC_LIMS: float = 0.0
AXIS_THRESHOLD: tuple = (5, 100)
CAMERA_FOV: int = 45
CAMERA_RESOLUTION: tuple = (256, 256)
DB_CAM_ALTITUDE: int = 300
DIAMLIMS: tuple = (4, 100)
FILLED: bool = True
GENERATION_KWARGS: dict = {'axis_threshold': (5, 100), 'resolution': (256, 256), 'fov': 45, 'min_sol_incidence': 10, 'max_sol_incidence': 80, 'filled': True, 'ellipse_limit': 1.3, 'arc_lims': 0.0, 'diamlims': (4, 100), 'instancing': True, 'randomized_orientation': True, 'mask_thickness': 1, 'save_craters': True}
INSTANCING: bool = True
KERNEL_ROOT: pathlib.PosixPath  # value = PosixPath('data/spice_kernels')
MASK_THICKNESS: int = 1
MAX_ELLIPTICITY: float = 1.3
MAX_SOL_INCIDENCE: int = 80
MIN_SOL_INCIDENCE: int = 10
RANDOMIZED_ORIENTATION: bool = True
RMOON: float = 1737.1
SAVE_CRATERS: bool = True
SPICE_BASE_URL: str = 'https://naif.jpl.nasa.gov/pub/naif/'
TRIAD_RADIUS: int = 200

from __future__ import annotations
import bpy as bpy
import bpy_types
import mathutils as mathutils
import numpy
import numpy as np
import os as os
from pylupnt import utils
from scipy.spatial.transform._rotation import Rotation as R
__all__ = ['ALBEDO', 'BODY', 'BODY_SECONDARY', 'CAMERA', 'CYCLES', 'EEVEE', 'ENGINE', 'R', 'R_ocv2ogl', 'SCALE_BU', 'SCENE', 'SUN', 'SUN_ENERGY', 'bpy', 'mathutils', 'np', 'os', 'render', 'update_scene', 'utils']
def render(r_c_pa, R_pa2c, r_s_pa, filepath, frame = 'OpenCV'):
    ...
def update_scene():
    ...
ALBEDO: float = 0.169
BODY: bpy_types.Object  # value = bpy.data.objects['Moon']
BODY_SECONDARY = None
CAMERA: bpy_types.Object  # value = bpy.data.objects['Camera']
CYCLES: str = 'CYCLES'
EEVEE: str = 'BLENDER_EEVEE'
ENGINE: str = 'BLENDER_EEVEE'
R_ocv2ogl: numpy.ndarray  # value = array([[ 1.,  0.,  0.],...
SCALE_BU: float = 0.001
SCENE: dict = {'fov': 60, 'resx': 1024, 'resy': 1024, 'encoding': 8, 'rendSamples': 64, 'viewSamples': 4, 'scattering': 0, 'labelDepth': 0, 'labelID': 0, 'labelSlopes': 0, 'viewtransform': 'Filmic', 'filmexposure': 1}
SUN: bpy_types.Object  # value = bpy.data.objects['Sun']
SUN_ENERGY: int = 30
_SCENE_OLD = None
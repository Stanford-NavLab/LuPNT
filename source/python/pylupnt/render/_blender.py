import logging

try:
    import bpy
except ImportError:
    pass
import os
import mathutils
import numpy as np
from .. import _pylupnt as _pnt
from scipy.spatial.transform import Rotation as R
from .. import utils


class Blender:

    EEVEE = "BLENDER_EEVEE"
    CYCLES = "CYCLES"

    ENGINE = EEVEE
    ALBEDO = 0.169
    SUN_ENERGY = 30
    SCALE_BU = 1e-3
    R_ocv2ogl = np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]], dtype=float)

    CAMERA = None
    SUN = None
    BODY = None
    BODY_SECONDARY = None

    SCENE = dict(
        fov=60,  # [deg]
        resx=1024,
        resy=1024,
        encoding=8,
        rendSamples=64,
        viewSamples=4,
        scattering=0,
        labelDepth=0,
        labelID=0,
        labelSlopes=0,
        viewtransform="Filmic",
        filmexposure=1,
    )
    _SCENE_OLD = None

    def __init__(self):
        bpy.ops.wm.open_mainfile(
            filepath=os.path.join(utils.LUPNT_DATA_PATH, "surface", "Moon.blend")
        )
        self.CAMERA = bpy.data.objects["Camera"]
        self.SUN = bpy.data.objects["Sun"]
        self.BODY = bpy.data.objects["Moon"]
        self.BODY_SECONDARY = None
        self.update_scene()

    def update_scene(self):
        # CAMERA properties
        self.CAMERA.data.type = "PERSP"
        self.CAMERA.data.lens_unit = "FOV"
        self.CAMERA.data.angle = self.SCENE["fov"] * np.pi / 180
        self.CAMERA.data.clip_start = 0.1  # [m]
        self.CAMERA.data.clip_end = 10000  # [m]
        if self.ENGINE == "CYCLES":
            bpy.context.scene.cycles.film_exposure = self.SCENE["filmexposure"]
        bpy.context.scene.view_settings.view_transform = self.SCENE["viewtransform"]
        bpy.context.scene.render.pixel_aspect_x = 1
        bpy.context.scene.render.pixel_aspect_y = 1
        bpy.context.scene.render.resolution_x = self.SCENE["resx"]  # CAM resolution (x)
        bpy.context.scene.render.resolution_y = self.SCENE["resy"]  # CAM resolution (y)
        bpy.context.scene.render.image_settings.color_mode = "BW"
        bpy.context.scene.render.image_settings.color_depth = str(
            self.SCENE["encoding"]
        )
        if self.ENGINE == "CYCLES":
            bpy.context.scene.cycles.diffuse_bounces = 0

        # SUN properties
        self.SUN.data.type = "SUN"
        self.SUN.data.energy = self.SUN_ENERGY  # To perform quantitative analysis
        self.SUN.data.angle = 0.53 * np.pi / 180

        # WORLD properties
        # black = (0, 0, 0, 1)
        # bpy.data.worlds["World"].node_tree.nodes["Background"].inputs[0].default_value = black

        # RENDERING ENGINE properties
        bpy.context.scene.render.engine = self.ENGINE
        if self.ENGINE == "CYCLES":
            bpy.context.scene.cycles.device = "GPU"
            bpy.context.scene.cycles.samples = self.SCENE["rendSamples"]
            bpy.context.scene.cycles.preview_samples = self.SCENE["viewSamples"]

        self._SCENE_OLD = self.SCENE.copy()

    def render(self, r_c_pa, R_pa2c, r_s_pa, filepath, frame="OpenCV"):
        if self.SCENE != self._SCENE_OLD:
            self.update_scene()

        SUN_DISTANCE = 3e6  # [km]

        if frame == "OpenCV":
            R_pa2ogl = self.R_ocv2ogl @ R_pa2c
        elif frame == "OpenGL":
            R_pa2ogl = R_pa2c
        q_c_pa = _pnt.rot2quat(R_pa2ogl.T)
        q_c_pa = q_c_pa[[3, 0, 1, 2]]
        r_m_pa = np.zeros(3)
        q_m_pa = np.array([1, 0, 0, 0])
        r_s_pa = np.array(r_s_pa) / np.linalg.norm(r_s_pa)

        bpy.context.scene.frame_current = 0

        self.BODY.rotation_mode = "QUATERNION"
        self.BODY.location = r_m_pa * self.SCALE_BU
        self.BODY.rotation_quaternion = q_m_pa

        self.CAMERA.rotation_mode = "QUATERNION"
        self.CAMERA.location = r_c_pa * self.SCALE_BU
        self.CAMERA.rotation_quaternion = q_c_pa

        self.SUN.rotation_mode = "QUATERNION"
        self.SUN.location = r_s_pa * SUN_DISTANCE * self.SCALE_BU
        self.SUN.rotation_quaternion = mathutils.Vector(r_s_pa).to_track_quat("Z", "Y")

        bpy.context.view_layer.update()
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        bpy.context.scene.render.filepath = os.path.join(filepath)
        bpy.ops.render.render(write_still=1)

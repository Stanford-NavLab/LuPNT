import bpy
import os
import mathutils
import numpy as np
from scipy.spatial.transform import Rotation as R
from .. import utils

bpy.ops.wm.open_mainfile(
    filepath=os.path.join(utils.LUPNT_DATA_PATH, "surface", "Moon.blend")
)

CYCLES = "CYCLES"
EEVEE = "BLENDER_EEVEE"

ENGINE = EEVEE
CAMERA = bpy.data.objects["Camera"]
SUN = bpy.data.objects["Sun"]
ALBEDO = 0.169  # TBD
SUN_ENERGY = 30  # TBD
BODY = bpy.data.objects["Moon"]
BODY_SECONDARY = None
SCALE_BU = 1e-3  # [km] to [Blender Units]
R_ocv2ogl = np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]], dtype=float)

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


def update_scene():
    # CAMERA properties
    CAMERA.data.type = "PERSP"
    CAMERA.data.lens_unit = "FOV"
    CAMERA.data.angle = SCENE["fov"] * np.pi / 180
    CAMERA.data.clip_start = 0.1  # [m]
    CAMERA.data.clip_end = 10000  # [m]
    if ENGINE == "CYCLES":
        bpy.context.scene.cycles.film_exposure = SCENE["filmexposure"]
    bpy.context.scene.view_settings.view_transform = SCENE["viewtransform"]
    bpy.context.scene.render.pixel_aspect_x = 1
    bpy.context.scene.render.pixel_aspect_y = 1
    bpy.context.scene.render.resolution_x = SCENE["resx"]  # CAM resolution (x)
    bpy.context.scene.render.resolution_y = SCENE["resy"]  # CAM resolution (y)
    bpy.context.scene.render.image_settings.color_mode = "BW"
    bpy.context.scene.render.image_settings.color_depth = str(SCENE["encoding"])
    if ENGINE == "CYCLES":
        bpy.context.scene.cycles.diffuse_bounces = 0

    # SUN properties
    SUN.data.type = "SUN"
    SUN.data.energy = SUN_ENERGY  # To perform quantitative analysis
    SUN.data.angle = 0.53 * np.pi / 180

    # WORLD properties
    # black = (0, 0, 0, 1)
    # bpy.data.worlds["World"].node_tree.nodes["Background"].inputs[0].default_value = black

    # RENDERING ENGINE properties
    bpy.context.scene.render.engine = ENGINE
    if ENGINE == "CYCLES":
        bpy.context.scene.cycles.device = "GPU"
        bpy.context.scene.cycles.samples = SCENE["rendSamples"]
        bpy.context.scene.cycles.preview_samples = SCENE["viewSamples"]

    _SCENE_OLD = SCENE.copy()


update_scene()


def render(r_c_pa, R_pa2c, r_s_pa, filepath, frame="OpenCV"):
    if SCENE != _SCENE_OLD:
        update_scene()

    SUN_DISTANCE = 3e6  # [km]

    if frame == "OpenCV":
        R_pa2ogl = R_ocv2ogl @ R_pa2c
    elif frame == "OpenGL":
        R_pa2ogl = R_pa2c
    q_c_pa = R.from_matrix(R_pa2ogl.T).as_quat()
    q_c_pa = q_c_pa[[3, 0, 1, 2]]
    r_m_pa = np.zeros(3)
    q_m_pa = np.array([1, 0, 0, 0])
    r_s_pa = np.array(r_s_pa) / np.linalg.norm(r_s_pa)

    bpy.context.scene.frame_current = 0

    BODY.rotation_mode = "QUATERNION"
    BODY.location = r_m_pa * SCALE_BU
    BODY.rotation_quaternion = q_m_pa

    CAMERA.rotation_mode = "QUATERNION"
    CAMERA.location = r_c_pa * SCALE_BU
    CAMERA.rotation_quaternion = q_c_pa

    SUN.rotation_mode = "QUATERNION"
    SUN.location = r_s_pa * SUN_DISTANCE * SCALE_BU
    SUN.rotation_quaternion = mathutils.Vector(r_s_pa).to_track_quat("Z", "Y")

    bpy.context.view_layer.update()
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    bpy.context.scene.render.filepath = os.path.join(filepath)
    bpy.ops.render.render(write_still=1)

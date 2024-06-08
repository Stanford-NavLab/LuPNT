import os
import numpy as np
import cv2
from scipy.ndimage.morphology import binary_dilation
import pylupnt as pnt
from PIL import Image
from .common.robbins import (
    load_craters,
    extract_robbins_dataset,
)

from .common.conics import (
    plot_conics,
    ellipse_axes,
    scale_det,
    conic_center,
    ellipse_angle,
    conic_matrix,
    project_crater_centers,
)

from pylupnt.crater_detection.common.coordinates import ENU_system, nadir_attitude

lat_crater, lon_crater, a_crater, b_crater, psi_crater, id_crater = (
    extract_robbins_dataset(
        load_craters(
            pnt.find_file("lunar_crater_database_robbins_2018.csv"),
            diamlims=[20, 250],
            ellipse_limit=1.4,
        )
    )
)


def crater_detection(
    r_cam_pa,
    R_pa2cam,
    r_sun_pa,
    fov_cam,
    res_cam,
    seed=None,
    alpha_sun=-0.05,
    alpha_cam=0.5,
):
    if seed is not None:
        np.random.seed(seed)
    r_M = r_cam_pa[:, None]
    R_cam2pa = R_pa2cam.T

    r_geo = np.vstack((lat_crater, lon_crater, np.zeros_like(lat_crater))).T
    p_Mi = pnt.geographical_to_cartesian(r_geo, pnt.R_MOON)[:, :, None]

    # Visibility
    height = np.linalg.norm(r_M) - pnt.R_MOON
    dist = height / np.cos(np.radians(fov_cam / 2) * np.sqrt(2) * 1.5)
    vis_crater = (np.sqrt(np.sum(np.square(p_Mi - r_M), axis=1)) < dist).ravel()
    # Sun angle
    e1 = p_Mi.squeeze() / np.linalg.norm(p_Mi.squeeze(), axis=1)[:, None]
    e2 = r_sun_pa / np.linalg.norm(r_sun_pa)
    vis_crater &= np.einsum("ij,j->i", e1, e2) > alpha_sun
    # Spacecraft angle
    e1 = p_Mi.squeeze() / np.linalg.norm(p_Mi.squeeze(), axis=1)[:, None]
    e2 = r_cam_pa / np.linalg.norm(r_cam_pa)
    vis_crater &= np.einsum("ij,j->i", e1, e2) > alpha_cam
    # Camera frame
    r_center = project_crater_centers(p_Mi[vis_crater], fov_cam, res_cam, R_cam2pa, r_M)
    margin = 0.0
    vis_crater2 = np.logical_and.reduce(
        (
            r_center[:, 0] > margin * res_cam[0],
            r_center[:, 0] < (1 - margin) * res_cam[0],
            r_center[:, 1] > margin * res_cam[1],
            r_center[:, 1] < (1 - margin) * res_cam[1],
        )
    )
    r_center = r_center[vis_crater2]
    vis_crater[vis_crater] &= vis_crater2
    p_Mi = p_Mi[vis_crater]

    T_EM = np.concatenate(ENU_system(p_Mi), axis=-1)
    S = np.concatenate((np.identity(2), np.zeros((1, 2))), axis=0)
    H_Mi = np.concatenate((T_EM @ S, p_Mi), axis=-1)

    f_x = (res_cam[0] / 2) / np.tan(np.radians(fov_cam) / 2)
    f_y = (res_cam[1] / 2) / np.tan(np.radians(fov_cam) / 2)
    K = np.array(
        [
            [f_x, 0, res_cam[0] / 2],
            [0, f_y, res_cam[1] / 2],
            [0, 0, 1],
        ]
    )

    k = np.array([0.0, 0.0, 1.0])[:, None]
    r_C = R_pa2cam @ r_M
    # P_MC = K @ np.concatenate((T_MC, -r_C), axis=1)
    P_MC = K @ R_pa2cam @ np.concatenate((np.identity(3), -r_M), axis=1)
    H_Ci = P_MC @ np.concatenate(
        (H_Mi, np.tile(k.T[None, ...], (len(H_Mi), 1, 1))), axis=1
    )

    C_i = conic_matrix(
        a_crater[vis_crater], b_crater[vis_crater], psi_crater[vis_crater]
    )
    A_i = np.linalg.inv(H_Ci).transpose((0, 2, 1)) @ C_i @ np.linalg.inv(H_Ci)

    sigma_pix = 1
    A_craters = A_i

    n_det = len(A_craters)
    major_det, minor_det = ellipse_axes(A_craters)
    psi_det = ellipse_angle(A_craters)
    r_craters_det = conic_center(A_craters)

    major_det += np.random.uniform(-sigma_pix, sigma_pix, size=n_det)
    minor_det += np.random.uniform(-sigma_pix, sigma_pix, size=n_det)
    psi_det += np.random.normal(scale=(20 / 180) * np.pi, size=n_det)
    r_craters_det += np.random.normal(-sigma_pix, sigma_pix, size=r_craters_det.shape)

    A_detections = conic_matrix(major_det, minor_det, psi_det, *r_craters_det.T)

    T_EM = np.concatenate(ENU_system(p_Mi), axis=-1)
    T_ME = np.linalg.inv(T_EM)

    B_craters = R_cam2pa @ K.T @ A_detections @ K @ R_pa2cam
    A_lstsq = (S.T @ T_ME @ B_craters).reshape(-1, 3)
    b_lstsq = (S.T @ T_ME @ B_craters @ p_Mi).reshape(-1)
    r_lstsq = np.linalg.lstsq(A_lstsq, b_lstsq, rcond=None)[0]
    return r_lstsq


def trace_lines_and_record(img, u_illum, N_lines, threshold, return_lines=False):
    img_array = np.array(img)
    height, width = img_array.shape
    mask = np.zeros_like(img_array)
    mask_lines = np.zeros_like(img_array)

    # Calculate starting points along the borders
    start_points = []
    start_points.extend(
        [(0, int(y)) for y in np.linspace(0, height, N_lines // 2, endpoint=False)]
    )
    start_points.extend(
        [(int(x), 0) for x in np.linspace(0, width, N_lines // 2, endpoint=False)]
    )
    for start in start_points:
        x, y = start
        while 0 <= x < width and 0 <= y < height:
            if img_array[y, x] > threshold:
                mask[y, x] = 1
                break
            else:
                mask_lines[y, x] = 1
            x += int(np.round(u_illum[0]))
            y += int(np.round(u_illum[1]))

    if return_lines:
        return mask, mask_lines

    return mask


def horizon_matching(
    u: np.ndarray,
    K_inv: np.ndarray,
    R_cam2pa: np.ndarray,
    a: float,
    b: float = None,
    c: float = None,
):
    if b is None:
        b = a
    if c is None:
        c = a
    D = np.diag([1 / a, 1 / b, 1 / c])
    R = D @ R_cam2pa @ K_inv
    H = u @ R.T
    H /= np.linalg.norm(H, axis=-1)[:, None]
    n = np.linalg.lstsq(H, np.ones((len(H), 1)), rcond=None)[0]
    R_pa2cam = R_cam2pa.T
    D_inv = np.diag([a, b, c])
    rC = (n.T @ n - 1) ** -0.5 * R_pa2cam @ D_inv @ n
    rP = R_cam2pa @ rC
    return np.ravel(rP)


def horizon_detection(img, R_pa2cam, r_sun_pa, K, N_lines=200, threshold=110):
    R_cam2pa = R_pa2cam.T

    # Sun direction
    u_illum = -np.array([[1, 0, 0], [0, 1, 0]]) @ R_pa2cam @ r_sun_pa
    u_illum /= np.linalg.norm(u_illum)
    u_illum

    mask, mask_lines = trace_lines_and_record(
        img, u_illum, N_lines, threshold, return_lines=True
    )

    # Canny edge detection
    import cv2

    edges = cv2.Canny(np.array(img), 100, 200)

    # Dilate the mask
    structure = np.ones((3, 3))
    dilated_mask = binary_dilation(mask, structure=structure).astype(mask.dtype)
    final_mask = dilated_mask * edges

    uv_horizon = np.ones((np.count_nonzero(final_mask), 3))
    uv_horizon[:, :2] = np.vstack(np.where(final_mask)).T
    r_cam_pa_horizon = horizon_matching(
        uv_horizon, np.linalg.inv(K), R_cam2pa, pnt.R_MOON
    )
    return r_cam_pa_horizon

import os
from typing import List
from matplotlib.patches import Ellipse
from matplotlib.collections import LineCollection
import matplotlib
from matplotlib import cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import numpy as np
import open3d as o3d
import plotly.graph_objects as go
import torch

from pylupnt.math import ensure_torch
from pylupnt import Logger

COLORBAR_FRAC = 0.027


def get_colormap(colormap: str, num_labels: int):
    """common colormaps: tab20, jet"""
    colors = matplotlib.colormaps[colormap](np.linspace(0, 1, num_labels))
    return mcolors.ListedColormap(colors)


def get_GnYlRd_colormap():
    """Get a green-yellow-red colormap (good for errors/quality).

    Green = low values (good), Yellow = medium, Red = high values (bad)
    """
    colors = ["green", "yellow", "red"]
    return mcolors.LinearSegmentedColormap.from_list("GnYlRd", colors)


def plot_colored_line_by_value(
    x, y, values, cmap="GnYlRd", vmin=None, vmax=None, colorbar=True, cbar_label=None, **kwargs
):
    """Plot a line colored by associated values (e.g., errors).

    Args:
        x, y: Coordinates of the line
        values: Array of values to map to colors (same length as x, y)
        cmap: Colormap name or object (default: 'gnrlrd' green-yellow-red)
        vmin, vmax: Value range for color mapping. If None, uses min/max of values
        colorbar: If True, add a colorbar
        cbar_label: Label for the colorbar
        **kwargs: Additional arguments passed to LineCollection

    Returns:
        LineCollection object
    """
    if cmap == "GnYlRd":
        cmap = get_GnYlRd_colormap()

    if vmin is None:
        vmin = np.nanmin(values)
    if vmax is None:
        vmax = np.nanmax(values)

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    norm = plt.Normalize(vmin, vmax)
    lc = LineCollection(segments, cmap=cmap, norm=norm, **kwargs)
    lc.set_array(values[:-1])  # Use values for segment colors
    lc.set_linewidth(kwargs.get("linewidth", 2))

    plt.gca().add_collection(lc)

    if colorbar:
        cbar = plt.colorbar(lc, ax=plt.gca())
        if cbar_label:
            cbar.set_label(cbar_label)

    return lc


def use_latex():
    plt.rcParams.update(
        {
            "text.usetex": True,
            "font.family": "serif",
            "font.serif": ["Computer Modern Roman"],
        }
    )


def sample_colors(colormap: str, num_labels: int):
    """common colormaps: tab20, jet"""
    return np.array(
        [matplotlib.colormaps[colormap](i / max(num_labels - 1, 1))[:3] for i in range(num_labels)]
    )


class RemoteDisplay:
    """
    A minimal context manager to configure the environment for Open3D
    visualization in a remote Jupyter notebook.

    It sets the DISPLAY variable and temporarily disables Open3D's
    WebVisualizer detection.
    """

    def __init__(self, display: int = 1):
        self.display = display
        self.original_display = None
        self.original_jupyter_url = None

    def __enter__(self):
        if "DISPLAY" in os.environ:
            self.original_display = int(os.environ.pop("DISPLAY")[1:])
        if "JUPYTER_SERVER_URL" in os.environ:
            self.original_jupyter_url = os.environ.pop("JUPYTER_SERVER_URL")
        os.environ["DISPLAY"] = f":{self.display}"
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.original_jupyter_url is not None:
            os.environ["JUPYTER_SERVER_URL"] = self.original_jupyter_url
        if self.original_display is not None:
            os.environ["DISPLAY"] = f":{self.original_display}"
        elif "DISPLAY" in os.environ:
            del os.environ["DISPLAY"]
        return False


def plot_mapper_cameras(dataset, mapper, cam_name, idxs, cols=2):
    rows = int(np.ceil(len(idxs) / cols))
    fig, axs = plt.subplots(rows, 2 * cols, figsize=(2.5 * 2 * cols, 2.5 * rows))
    axs = axs.flatten()

    for i, idx in enumerate(idxs):
        data = dataset.get_data(idx)
        cam_data = data["cameras"][cam_name]

        wTb = data["pose"]
        bTc = cam_data["body_T_cam"]
        wTc = wTb @ bTc
        rgb_gt = ensure_torch(cam_data["rgb"], channels_first=False)
        intrinsics = cam_data["intrinsics"]

        outputs = mapper.render(wTc, intrinsics)
        rgb = outputs["rgb"][0]
        depth = outputs["depth"][0]
        axs[i * 2].imshow(rgb_gt.cpu())
        axs[i * 2 + 1].imshow(rgb.cpu())
        # axs[i * 3 + 2].imshow(depth.cpu())
        axs[i * 2].set_axis_off()
        axs[i * 2 + 1].set_axis_off()
    for i in range(len(idxs) * 3, len(axs)):
        axs[i].remove()
    plt.tight_layout()
    return fig, axs


def plot_grid(
    inputs, titles=None, repeat_title=False, figsize=(10, 10), axis=False, show=False
) -> tuple:
    n, m = len(inputs[0]), len(inputs)
    fig, axs = plt.subplots(n, m, figsize=figsize)
    if n == 1:
        axs = np.array([axs])
    if m == 1:
        axs = np.array([axs])
    for i in range(n):
        for j in range(m):
            axs[i, j].imshow(inputs[j][i])
            if titles is not None and (i == 0 or repeat_title):
                axs[i, j].set_title(titles[j])
            if not axis:
                axs[i, j].set_axis_off()
    plt.tight_layout()
    if show:
        plt.show()
    return fig, axs


def plot_image_grid(images, titles=None, figsize=(10, 6), cmap="gray", show_axes=False):
    # 1D list -> single row; 2D list -> rows-first grid
    is_1d = isinstance(images[0], np.ndarray)
    grid = [images] if is_1d else images
    rows, cols = (1, len(images)) if is_1d else (len(images), len(images[0]))

    fig, axs = plt.subplots(rows, cols, figsize=figsize)

    # normalize axs to 2D array
    if rows == 1 and cols == 1:
        axs = np.array([[axs]])
    elif rows == 1:
        axs = axs[None, :]
    elif cols == 1:
        axs = axs[:, None]

    for r in range(rows):
        for c in range(cols):
            ax = axs[r, c]
            img = grid[r][c]
            ax.imshow(img, cmap=cmap if img.ndim == 2 else None)
            if titles is not None:
                ax.set_title(titles[c] if is_1d else titles[r][c])
            if not show_axes:
                ax.set_axis_off()

    fig.tight_layout()
    return fig, axs


def plot_field(depth, ax=None, cbar=True, cmap="viridis", vmin=None, vmax=None) -> tuple:
    depth_np = depth if isinstance(depth, np.ndarray) else depth.cpu().numpy()
    if ax is None:
        ax = plt.gca()
    im = ax.imshow(depth_np, cmap=cmap, vmin=vmin, vmax=vmax)
    if cbar:
        cb = ax.figure.colorbar(im, ax=ax, fraction=COLORBAR_FRAC)
    return im, cb


def plot_horizons(img: np.ndarray, horizons: list[np.ndarray], ax=None, colors=None) -> tuple:
    if ax is None:
        ax = plt.gca()

    # Create a color image to support multiple colors
    horizon_img = np.zeros((*img.shape[:2], 3), dtype=np.uint8)

    # Define colors for each horizon (RGB values)
    if colors is None:
        colors = [
            [255, 0, 0],  # Red
            [0, 255, 0],  # Green
            [0, 0, 255],  # Blue
            [255, 255, 0],  # Yellow
            [255, 0, 255],  # Magenta
            [0, 255, 255],  # Cyan
            [255, 128, 0],  # Orange
            [128, 0, 255],  # Purple
        ]

    for i, horizon in enumerate(horizons):
        color = colors[i % len(colors)]  # Cycle through colors if more horizons than colors
        for c, r in enumerate(horizon):
            if 0 <= r < horizon_img.shape[0] and 0 <= c < horizon_img.shape[1]:
                horizon_img[r, c] = color

    ax.imshow(horizon_img)


def plot_path_3d(
    path: np.ndarray,
    fig=None,
    color="red",
    markersize=3,
    width=3,
    markers=True,
    **kwargs,
) -> go.Figure:
    """Plot a 3D path with optional markers."""
    if fig is None:
        fig = go.Figure()
    if markers:
        fig.add_scatter3d(
            x=path[:, 0],
            y=path[:, 1],
            z=path[:, 2],
            mode="markers+lines",
            marker=dict(size=markersize, color=color),
            line=dict(color=color, width=width),
            hovertext=np.arange(len(path)),
            **kwargs,
        )
    else:
        fig.add_scatter3d(
            x=path[:, 0],
            y=path[:, 1],
            z=path[:, 2],
            mode="lines",
            line=dict(color=color, width=width),
            hovertext=np.arange(len(path)),
            **kwargs,
        )
    return fig


def plot_path_2d(path: np.ndarray, fig=None, color="blue", linewidth=2, **kwargs):
    """Plot a 2D path."""
    if fig is None:
        fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=path[:, 0],
            y=path[:, 1],
            mode="lines",
            line=dict(color=color, width=linewidth),
            **kwargs,
        )
    )
    fig.update_layout(
        width=900,
        height=900,
        xaxis=dict(scaleanchor="y"),
        xaxis_title="X (m)",
        yaxis_title="Y (m)",
        font=dict(size=25),
    )
    return fig


def plot_heatmap(
    data: np.ndarray, fig=None, colorscale="Viridis", no_axes=False, scale_title="Z (m)", **kwargs
) -> go.Figure:
    """Plot a 2D heatmap."""
    if fig is None:
        fig = go.Figure()

    if data.ndim == 3:
        fig.add_trace(
            go.Heatmap(
                x=data[0, :, 0],
                y=data[:, 0, 1],
                z=data[:, :, 2],
                colorscale=colorscale,
                colorbar=dict(
                    lenmode="fraction",
                    len=0.8,
                    title=scale_title,
                ),
                **kwargs,
            )
        )

    else:
        fig.add_trace(go.Heatmap(z=data, colorscale=colorscale, **kwargs))

    fig.update_xaxes(
        visible=not no_axes,
        showgrid=False,
        zeroline=False,
        constrain="domain",  # no padding from axis
        scaleanchor="y",
    )
    fig.update_yaxes(visible=not no_axes, showgrid=False, zeroline=False, constrain="domain")

    if no_axes:
        fig.update_layout(xaxis=dict(visible=False), yaxis=dict(visible=False))
    fig.update_layout(
        width=800,
        height=800,
        xaxis=dict(scaleanchor="y"),
        xaxis_title="X (m)",
        yaxis_title="Y (m)",
    )
    return fig


def pose_trace(
    pose: np.ndarray | tuple,
    name: str = "",
    style: str = "solid",
    width: int = 5,
    length: float = 1.0,
) -> list:
    """Create a plotly trace to visualize a pose

    RGB vectors are used to represent the X, Y, Z axes of the rotation matrix

    Parameters
    ----------
    pose : 4x4 np.ndarray or tuple

    Returns
    -------
    traces : list
        List of plotly traces for displaying the pose

    """
    # If pose is tuple
    if isinstance(pose, tuple):
        # Unpack pose into rotation matrix R and translation vector t
        R, t = pose
        t = np.array(t)  # Ensure t is a numpy array
    else:  # If pose is a 4x4 matrix
        R = pose[:3, :3]
        t = pose[:3, 3]

    # Define arrow colors for each axis (RGB)
    colors = ["red", "green", "blue"]

    # Define the unit vectors from the columns of R
    axis_vectors = [R[:, 0], R[:, 1], R[:, 2]]

    # Create traces for each axis (X, Y, Z)
    traces = []
    for i, vec in enumerate(axis_vectors):
        arrow_start = t
        arrow_end = t + vec * length  # Arrow points in the direction of the column of R

        # Create an arrow trace for the axis
        if name == "":
            trace = go.Scatter3d(
                x=[arrow_start[0], arrow_end[0]],
                y=[arrow_start[1], arrow_end[1]],
                z=[arrow_start[2], arrow_end[2]],
                mode="lines",
                marker=dict(size=4),
                line=dict(color=colors[i], width=width, dash=style),
                showlegend=False,
            )
        else:
            AXES_NAMES = ["X", "Y", "Z"]
            trace = go.Scatter3d(
                x=[arrow_start[0], arrow_end[0]],
                y=[arrow_start[1], arrow_end[1]],
                z=[arrow_start[2], arrow_end[2]],
                mode="lines",
                marker=dict(size=4),
                line=dict(color=colors[i], width=width, dash=style),
                name=name + f"_{AXES_NAMES[i]}",
                showlegend=True,
            )
        traces.append(trace)

    return traces


def pose_traces(pose_list: list | np.ndarray, **kwargs) -> list:
    """Create traces for a sequence of poses

    Parameters
    ----------
    pose_list : list of tuples
        List of poses, where each pose is a tuple (R, t)

    Returns
    -------
    all_traces : list
        List of plotly traces for displaying all poses

    """
    all_traces = []

    for pose in pose_list:
        traces = pose_trace(pose, **kwargs)
        all_traces.extend(traces)

    return all_traces


def plot_poses(poses: list | np.ndarray, fig=None, no_axes=False, **kwargs) -> go.Figure:
    """poses is a list of 4x4 arrays or an Nx4x4 array"""
    if fig is None:
        fig = go.Figure()
    if len(poses) > 100 and not no_axes:
        Logger.warning(
            "Plotting axes for more than 100 poses is not recommended. Defaulting to no_axes=True"
        )
        no_axes = True
    if no_axes:
        positions = np.array([pose[:3, 3] for pose in poses])
        fig = plot_path_3d(positions, fig=fig, **kwargs)
    else:
        fig.add_traces(pose_traces(poses, **kwargs))
    fig.update_layout(
        width=1600,
        height=900,
        scene_aspectmode="data",
        scene=dict(xaxis_title="X (m)", yaxis_title="Y (m)", zaxis_title="Z (m)"),
    )
    return fig


def plot_box2d(xlims: np.ndarray, ylims: np.ndarray, center: np.ndarray = None, **kwargs):
    if center is None:
        center = np.zeros(2)
    plt.plot(
        [xlims[0], xlims[1], xlims[1], xlims[0], xlims[0]],
        [ylims[0], ylims[0], ylims[1], ylims[1], ylims[0]],
        **kwargs,
    )


def plot_mesh(
    fig: go.Figure,
    mesh: o3d.geometry.TriangleMesh,
    lighting=None,
    light_position=None,
    colorscale=None,
    resolution=None,
):
    if colorscale is None:
        colorscale = [0, "rgb(153, 153, 153)"], [1.0, "rgb(160,160,160)"]
    if light_position is None:
        light_position = dict(x=1000, y=500, z=2000)
    if lighting is None:
        lighting = dict(
            ambient=0.02,
            diffuse=0.8,
            specular=0.15,
            roughness=0.5,
            fresnel=0.2,
            facenormalsepsilon=1e-15,
            vertexnormalsepsilon=1e-15,
        )

    if resolution is not None:
        mesh = mesh.simplify_vertex_clustering(resolution)

    vertices = np.asarray(mesh.vertices)
    triangles = np.asarray(mesh.triangles)

    mesh_plotly = go.Mesh3d(
        x=vertices[:, 0],
        y=vertices[:, 1],
        z=vertices[:, 2],
        i=triangles[:, 0],
        j=triangles[:, 1],
        k=triangles[:, 2],
        flatshading=True,
        colorscale=colorscale,
        intensity=vertices[:, 0],
        lighting=lighting,
        lightposition=light_position,
        showlegend=False,
        showscale=False,
    )
    fig.add_trace(mesh_plotly)
    return


def plot_colored_line(
    x,
    y,
    colors_or_cmap,
    label=None,
    linestyle="-",
    linewidth=2,
    colorbar=False,
    cbar_label=None,
    vmin=None,
    vmax=None,
):
    """Plot a line with colors varying along its length.

    Args:
        x, y: Coordinates of the line
        colors_or_cmap: Either a numpy array of colors for each segment, a colormap name string,
                       or a colormap object
        label: Optional label for legend
        linestyle: Line style
        linewidth: Width of the line
        colorbar: If True, add a colorbar
        cbar_label: Label for the colorbar
        vmin, vmax: Color scale limits (only used with colormap, not explicit colors)

    Returns:
        LineCollection object (can be used to add colorbar later)
    """
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Check if colors_or_cmap is an explicit color array or a colormap name
    if (
        isinstance(colors_or_cmap, (np.ndarray, list))
        and len(colors_or_cmap) > 0
        and isinstance(colors_or_cmap[0], (tuple, list, np.ndarray))
    ):
        # Explicit colors provided (array of RGBA tuples)
        colors = np.array(colors_or_cmap)
        lc = LineCollection(segments, colors=colors[:-1], linestyle=linestyle, linewidth=linewidth)
    else:
        # Colormap name or object provided
        if vmin is None:
            vmin = 0
        if vmax is None:
            vmax = len(x) - 1
        norm = plt.Normalize(vmin, vmax)
        lc = LineCollection(
            segments, cmap=colors_or_cmap, norm=norm, linestyle=linestyle, linewidth=linewidth
        )
        lc.set_array(np.arange(len(x)))

    plt.gca().add_collection(lc)

    if label:
        if (
            isinstance(colors_or_cmap, (np.ndarray, list))
            and len(colors_or_cmap) > 0
            and isinstance(colors_or_cmap[0], (tuple, list, np.ndarray))
        ):
            # Use first color for legend
            plt.plot([], [], color=colors[0], label=label, linestyle=linestyle)
        else:
            cmap_obj = (
                plt.get_cmap(colors_or_cmap) if isinstance(colors_or_cmap, str) else colors_or_cmap
            )
            plt.plot([], [], color=cmap_obj(0.7), label=label, linestyle=linestyle)

    if colorbar:
        cbar = plt.colorbar(lc, ax=plt.gca())
        if cbar_label:
            cbar.set_label(cbar_label)

    return lc


def plot_meshes(fig: go.Figure, meshes: List[o3d.geometry.TriangleMesh], **kwargs):
    for mesh in meshes:
        plot_mesh(fig, mesh, **kwargs)


def plot_covariances(xs, Ps, **kwargs):
    ax = plt.gca()
    for x, P in zip(xs, Ps):
        confidence_ellipse(x, P, ax, **kwargs)


def confidence_ellipse(x, P, ax, n_std=3.0, edgecolor="black", facecolor="none", **kwargs):
    if isinstance(P, torch.Tensor):
        P = P.detach().cpu().numpy()

    cov = P[:2, :2]
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensional dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse(
        (0, 0),
        width=ell_radius_x * 2,
        height=ell_radius_y * 2,
        edgecolor=edgecolor,
        facecolor=facecolor,
        **kwargs,
    )

    # Calculating the standard deviation of x from the squareroot of the variance and multiplying with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    scale_y = np.sqrt(cov[1, 1]) * n_std

    transf = transforms.Affine2D().rotate_deg(45).scale(scale_x, scale_y).translate(x[0], x[1])

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


def plot_depth(depth, depth_gt, name, max_error=0.1, figsize=(8, 8), show=True):
    fig, axs = plt.subplots(2, 2, figsize=figsize)
    depth_np = depth.cpu().numpy() if isinstance(depth, torch.Tensor) else depth
    depth_gt_np = depth_gt.cpu().numpy() if isinstance(depth_gt, torch.Tensor) else depth_gt

    plot_field(depth, ax=axs[0, 0])
    axs[0, 0].set_title(f"{name} [m]")
    plot_field(depth_gt, ax=axs[1, 0])
    axs[1, 0].set_title("Ground truth [m]")

    error = np.abs(depth_np - depth_gt_np)
    plot_field(error, ax=axs[0, 1], cmap="jet")
    axs[0, 1].set_title("Error [m]")
    plot_field(np.where(error < max_error, error, np.nan), ax=axs[1, 1], cmap="jet")
    axs[1, 1].set_title(f"Error below {max_error:.2f} m")

    idxs = np.isfinite(error)
    mae = np.mean(error[idxs])
    absrel = np.mean(error[idxs] / depth_gt_np[idxs])
    d25 = np.mean(
        np.maximum(depth_gt_np[idxs] / depth_np[idxs], depth_np[idxs] / depth_gt_np[idxs]) < 1.25
    )
    print(f"MAE: {mae:.2f} m")
    print(f"AbsRel: {absrel:.2f}")
    print(f"d25: {np.mean(d25):.2f}")

    plt.tight_layout()
    if show:
        plt.show()
    else:
        return fig


def plot_surface(
    grid: np.ndarray, fig=None, colorscale="Viridis", no_axes=False, showscale=True, **kwargs
) -> go.Figure:
    """
    grid is NxNx3 array representing the coordinates and elevation data for a surface plot.

    """
    if fig is None:
        fig = go.Figure()
    # Downsample large grids to 500x500 for plotting
    M, N = grid.shape[:2]
    if M > 500 or N > 500:
        downsample_factor_M = max(1, M // 500)
        downsample_factor_N = max(1, N // 500)
        grid = grid[::downsample_factor_M, ::downsample_factor_N, :]
        print(f"Downsampling grid from {(M, N)} to {grid.shape[:2]} for plotting")
    fig.add_trace(
        go.Surface(
            x=grid[:, :, 0],
            y=grid[:, :, 1],
            z=grid[:, :, 2],
            colorscale=colorscale,
            showscale=showscale,
            **kwargs,
        )
    )
    fig.update_layout(
        width=1600,
        height=900,
        scene_aspectmode="data",
        scene=dict(xaxis_title="X (m)", yaxis_title="Y (m)", zaxis_title="Z (m)"),
    )
    if no_axes:
        fig.update_layout(
            scene=dict(
                xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False)
            )
        )
    return fig


def plot_3d_points(
    points: np.ndarray, fig=None, color="blue", markersize=3, name=None, symbol="circle"
) -> go.Figure:
    """Plot 3D points."""
    if fig is None:
        fig = go.Figure()
    N = len(points)
    MAX_POINTS = 1e5
    if N > MAX_POINTS:
        downsample_factor = int(N / MAX_POINTS)
        points = points[::downsample_factor]
        print(f"Downsampling points from {N} to {len(points)} for plotting")
    fig.add_trace(
        go.Scatter3d(
            x=points[:, 0],
            y=points[:, 1],
            z=points[:, 2],
            mode="markers",
            marker=dict(size=markersize, color=color, symbol=symbol),
            name=name,
        )
    )
    fig.update_layout(
        width=1200,
        height=900,
        scene_aspectmode="data",
        scene=dict(xaxis_title="X (m)", yaxis_title="Y (m)", zaxis_title="Z (m)"),
    )
    return fig


def plot_loop_closures(trajectory, loop_closures: list, fig=None, **kwargs):
    """Plot loop closures

    trajectory: sequence of poses or positions
    loop_closures: list of tuples

    """
    if fig is None:
        fig = go.Figure()
    for i, j in loop_closures:
        fig.add_trace(
            go.Scatter3d(
                x=[trajectory[i][0, 3], trajectory[j][0, 3]],
                y=[trajectory[i][1, 3], trajectory[j][1, 3]],
                z=[trajectory[i][2, 3], trajectory[j][2, 3]],
                mode="markers+lines",
                marker=dict(color="red", size=5),
                line=dict(color="red", width=5),
                name=f"LC {i}-{j}",
                **kwargs,
            )
        )
    fig.update_layout(width=1600, height=900, scene_aspectmode="data")
    return fig


def plot_rock_height_maps(
    height_map,
    height_map_gt,
    rock_map,
    rock_map_gt,
    xlims,
    ylims,
    gt_color="orange",
    pred_color="red",
    log_scale=False,
    poses=None,
    show=True,
    plot_height_error=True,
):
    height_error = np.abs(height_map - np.where(rock_map_gt, np.nan, height_map_gt)) * 100
    mean = np.nanmean(height_error)
    std = np.nanstd(height_error)
    height_error = np.where(height_error > mean + 3 * std, np.nan, height_error)

    fig = plt.figure()
    kwargs = dict(origin="lower", extent=[xlims[0], xlims[1], ylims[0], ylims[1]])
    plt.imshow(
        np.where(rock_map_gt, 0, np.nan),
        cmap=plt.cm.colors.ListedColormap([gt_color]),
        **kwargs,
    )
    plt.imshow(
        np.where(rock_map, 0.9, np.nan),
        cmap=plt.cm.colors.ListedColormap([pred_color]),
        **kwargs,
    )
    if plot_height_error:
        if log_scale:
            plt.imshow(np.log10(height_error), **kwargs)
            min_exp = int(np.floor(np.log10(np.nanmin(height_error[height_error > 0]))))
            max_exp = int(np.ceil(np.log10(np.nanmax(height_error))))
            ticks = np.linspace(min_exp, max_exp, max_exp - min_exp + 1)
            cbar = plt.colorbar(
                label="Height error [cm]",
                fraction=0.03,
                ticks=ticks,
                spacing="proportional",  # Add spacing between ticks
            )
            cbar.ax.set_yticklabels([f"$10^{{{int(t)}}}$" for t in ticks])
        else:
            plt.imshow(height_error, **kwargs)
            plt.colorbar(label="Height error [cm]", fraction=0.03)
    plt.plot([], [], "o", c=gt_color, label="Ground truth")
    plt.plot([], [], "o", c=pred_color, label="Predicted")
    if poses is not None:
        xt = np.vstack([rt[:2, 3] for rt in poses])
        plt.plot(xt[:, 0], xt[:, 1], "k-", label="Trajectory")
    plt.xlabel("X [m]")
    plt.ylabel("Y [m]")
    plt.legend(framealpha=1.0, loc="upper left")
    if show:
        plt.show()
    else:
        return fig


def plot_camera_frustums(poses: list | np.ndarray, intrinsics: dict, fig=None, **kwargs):
    if fig is None:
        fig = go.Figure()
    for pose, intrinsics in zip(poses, intrinsics):
        fig = plot_camera_frustum(pose, intrinsics, fig=fig, **kwargs)
    return fig


def plot_camera_frustum(
    pose: np.ndarray,
    intrinsics: dict,
    fig=None,
    frustum_length: float = 1.0,
    color: str = "black",
    opacity: float = 0.3,
    name: str = "Camera",
    show_axes: bool = True,
) -> go.Figure:
    """
    Plot a camera frustum given pose and camera intrinsics.

    Parameters
    ----------
    pose : np.ndarray
        4x4 transformation matrix representing camera pose
    intrinsics : dict
        Camera intrinsics dictionary with keys: 'fx', 'fy', 'cx', 'cy', 'W', 'H'
    fig : go.Figure, optional
        Existing plotly figure to add frustum to
    frustum_length : float, optional
        Length of the frustum pyramid, by default 1.0
    color : str, optional
        Color of the frustum, by default "blue"
    opacity : float, optional
        Opacity of the frustum faces, by default 0.3
    name : str, optional
        Name for the camera in the legend, by default "Camera"
    show_axes : bool, optional
        Whether to show camera coordinate axes, by default True

    Returns
    -------
    go.Figure
        Plotly figure with camera frustum
    """
    if fig is None:
        fig = go.Figure()

    # Extract rotation matrix and translation vector
    R = pose[:3, :3]
    t = pose[:3, 3]

    # Calculate field of view from intrinsics
    fx = intrinsics["fx"]
    fy = intrinsics["fy"]
    W = intrinsics["W"]
    H = intrinsics["H"]

    fov_x = 2 * np.arctan(W / (2 * fx))
    fov_y = 2 * np.arctan(H / (2 * fy))

    # Calculate frustum corners in camera frame
    # Camera convention: X=forward, Y=left, Z=up
    half_width = frustum_length * np.tan(fov_x / 2)
    half_height = frustum_length * np.tan(fov_y / 2)

    # Define frustum corners in camera frame (camera at origin, looking down +X)
    corners_cam = np.array(
        [
            [0, 0, 0],  # Camera center
            [frustum_length, -half_width, -half_height],  # Bottom-left
            [frustum_length, half_width, -half_height],  # Bottom-right
            [frustum_length, half_width, half_height],  # Top-right
            [frustum_length, -half_width, half_height],  # Top-left
        ]
    )

    # Transform corners to world frame
    corners_world = (R @ corners_cam.T).T + t

    # Define frustum edges as line segments
    edges = [
        [0, 1],  # Camera center to bottom-left
        [0, 2],  # Camera center to bottom-right
        [0, 3],  # Camera center to top-right
        [0, 4],  # Camera center to top-left
        [1, 2],  # Bottom edge
        [2, 3],  # Right edge
        [3, 4],  # Top edge
        [4, 1],  # Left edge
    ]

    # Create line traces for frustum edges
    for edge in edges:
        start_point = corners_world[edge[0]]
        end_point = corners_world[edge[1]]

        fig.add_trace(
            go.Scatter3d(
                x=[start_point[0], end_point[0]],
                y=[start_point[1], end_point[1]],
                z=[start_point[2], end_point[2]],
                mode="lines",
                line=dict(color=color, width=2),
                showlegend=False,
                name=name if edge == edges[0] else None,  # Only show legend for first trace
            )
        )

    # Add camera coordinate axes if requested
    if show_axes:
        axis_length = frustum_length * 0.5
        axes_traces = pose_trace(pose, name=name, length=axis_length)
        fig.add_traces(axes_traces)

    return fig


def plot_loss_images(loss_images, batch, model_outputs, output_path):
    for j in range(len(loss_images["rgb_l1"])):
        fig, axs = plt.subplots(2, 5, figsize=(15, 5))
        axs = axs.flatten()
        i = 0

        axs[i].set_title("RGB (Ground Truth)")
        axs[i].imshow(batch["rgb"][j])
        i += 1

        axs[i].set_title("Depth (RAFT Stereo)")
        im = axs[i].imshow(batch["depth"][j])
        i += 1

        axs[i].set_title("Label (U-Net++)")
        axs[i].imshow(batch["label"][j])
        i += 1

        axs[i].set_title("RGB L1 Loss")
        axs[i].imshow(loss_images["rgb_l1"][j])
        i += 1

        axs[i].set_title("Empty Loss")
        im = axs[i].imshow(loss_images["empty"][j])
        fig.colorbar(im, ax=axs[i], fraction=0.046, pad=0.04)
        i += 1

        axs[i].set_title("RGB (GSPlat)")
        axs[i].imshow(model_outputs["rgb"][j].detach().cpu())
        i += 1

        axs[i].set_title("Depth (GSPlat)")
        im = axs[i].imshow(model_outputs["depth"][j].detach().cpu())
        fig.colorbar(im, ax=axs[i], fraction=0.046, pad=0.04)
        i += 1

        axs[i].set_title("Accumulation (GSPlat)")
        im = axs[i].imshow(model_outputs["accumulation"][j].detach().cpu())
        fig.colorbar(im, ax=axs[i], fraction=0.046, pad=0.04)
        i += 1

        axs[i].set_title("Depth L1 Loss")
        im = axs[i].imshow(loss_images["depth_l1"][j])
        fig.colorbar(im, ax=axs[i], fraction=0.046, pad=0.04)
        i += 1

        axs[i].set_title("Full Loss")
        im = axs[i].imshow(loss_images["full"][j])
        fig.colorbar(im, ax=axs[i], fraction=0.046, pad=0.04)
        i += 1

        plt.tight_layout()
        plt.savefig(output_path / f"loss_image{j:04d}.pdf")
        plt.show()


def plot_semantic_image(label_img):
    # from pylupnt.applications.segmentation import LABEL_COLORS_DICT
    # HACK: for current unreal data
    LABEL_COLORS_DICT = {
        0: [1.0, 1.0, 1.0],  # ground
        2: [1.0, 0.0, 1.0],  # rover
        3: [1.0, 0.251, 0.251],  # rock
        4: [1.0, 0.843, 0.0],  # lander
    }

    color_img = np.zeros((label_img.shape[0], label_img.shape[1], 3))
    for label, color in LABEL_COLORS_DICT.items():
        color_img[label_img == label] = color
    return color_img


def plot_keypoints(keypoints, labels=None, color="green", markersize=5):
    LABEL_COLORS_DICT = {
        0: [0.80, 0.80, 0.80],  # Ground – light gray, readable but subtle
        2: [0.00, 0.60, 1.00],  # Rover – bright cyan-blue, distinct and tech-like
        3: [1.00, 0.35, 0.00],  # Rock – vivid orange, pops against gray terrain
        4: [1.00, 0.90, 0.20],  # Lander – warm yellow-gold, clearly visible in shadow
    }
    if labels is not None:
        color = [LABEL_COLORS_DICT[label] for label in labels]
    ax = plt.gca()
    ax.scatter(keypoints[:, 0], keypoints[:, 1], c=color, s=markersize, marker="o")

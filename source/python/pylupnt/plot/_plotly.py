import os
import numpy as np
from PIL import Image
from sklearn.cluster import KMeans
from sklearn.utils import shuffle
from typing import Union
import plotly.graph_objs as go
from plotly.express.colors import sample_colorscale
import plotly.express as px
import plotly.io as pio
import pylupnt as pnt
from .. import utils

axis_dict = dict(
    mirror=True,
    ticks="outside",
    showline=True,
    showgrid=True,
    automargin=True,
)
pio.templates["lupnt"] = go.layout.Template(
    layout=go.Layout(
        font=dict(family="Computer Modern", size=18),
        xaxis=axis_dict.copy(),
        yaxis=axis_dict.copy(),
        colorway=px.colors.qualitative.D3,
        margin=dict(l=100, r=10, t=10, b=10),
    )
)
plotly_colors = px.colors.qualitative.D3
pio.templates.default = "plotly_white+lupnt"

MOON_SURFACE = Image.open(
    os.path.join(utils.LUPNT_DATA_PATH, "topo", "moon_surface.jpeg")
)
EARTH_SURFACE = Image.open(
    os.path.join(utils.LUPNT_DATA_PATH, "topo", "earth_surface.jpg")
)
IMAGES = {
    pnt.MOON: np.asarray(MOON_SURFACE),
    pnt.EARTH: np.asarray(EARTH_SURFACE),
}
RADII = {
    pnt.MOON: pnt.R_MOON,
    pnt.EARTH: pnt.R_EARTH,
}


def set_view(fig: go.Figure, azimuth: float, elevation: float, zoom: float = 1.0):
    eye = np.zeros(3) * 0.5
    eye[0] += np.cos(np.radians(azimuth)) * np.cos(np.radians(elevation))
    eye[1] += np.sin(np.radians(azimuth)) * np.cos(np.radians(elevation))
    eye[2] += np.sin(np.radians(elevation))
    eye *= zoom

    fig.update_layout(
        scene=dict(
            camera=dict(
                eye=dict(x=eye[0], y=eye[1], z=eye[2]),
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0),
            )
        )
    )


def set_equal_aspect_ratio(fig: go.Figure):
    lims = np.zeros((3, 2))
    # Iterate over the axis ranges
    for i, ax in enumerate(["x", "y", "z"]):
        lims[i] = fig.layout.scene[f"{ax}axis"].range
    # Get the max range
    max_range = np.max(lims[:, 1] - lims[:, 0])
    ratio = (lims[:, 1] - lims[:, 0]) / max_range
    # Set aspect ratio
    fig.update_layout(
        scene=dict(
            aspectmode="manual",
            aspectratio=dict(x=ratio[0], y=ratio[1], z=ratio[2]),
        )
    )


def create_sphere_meshgrid(
    rows: int, cols: int, radius: float = 1.0
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Create a sphere meshgrid

    Args:
        rows (int): number of rows
        cols (int): number of columns
    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray]: x, y, z coordinates of the sphere
    """
    u, v = np.meshgrid(
        np.linspace(0, 2 * np.pi, cols) - np.pi, np.linspace(0, np.pi, rows) - np.pi / 2
    )
    return (
        radius * np.cos(u) * np.cos(v),
        radius * np.sin(u) * np.cos(v),
        radius * np.sin(v),
    )


def plot_body(
    fig: go.Figure,
    body: pnt.NaifId,
    size_factor: int = 5,
    R_b2frame: np.ndarray = None,
    r_b2s: np.ndarray = None,
    r_body: np.ndarray = None,
    alpha: float = 0.2,
    scale: float = 3,
) -> go.Figure:
    """
    Plot a celestial body

    Args:
        body (int): celestial body
        size_factor (int): size factor
        R_b2frame (np.ndarray): rotation matrix from body to frame
        r_b2s (np.ndarray): vector from body to sun in the frame
        alpha (float): light intensity
    """
    img = IMAGES[body]
    radius = RADII[body]
    if r_body is None:
        r_body = np.zeros(3)
    reduced_img = img[::size_factor, ::size_factor]
    reduced_img = reduced_img[:, ::-1]
    r, c, _ = reduced_img.shape
    x, y, z = create_sphere_meshgrid(r, c, radius)
    xyz = np.array([x.flatten(), y.flatten(), z.flatten()]).T
    if r_b2s is not None:
        light = xyz @ r_b2s / np.linalg.norm(xyz, axis=1) / np.linalg.norm(r_b2s)
        reduced_img = reduced_img * (alpha + (1 - alpha) * light.reshape(r, c, 1))

    if R_b2frame is not None:
        xyz = xyz @ R_b2frame.T

    I, J, K, tri_color_intensity, pl_colorscale = mesh_data(
        reduced_img, n_colors=32, n_training_pixels=10000
    )

    r_body = r_body / 10**scale
    xyz = xyz / 10**scale

    fig.add_trace(
        go.Mesh3d(
            x=xyz[:, 0] + r_body[0],
            y=xyz[:, 1] + r_body[1],
            z=xyz[:, 2] + r_body[2],
            i=I,
            j=J,
            k=K,
            intensity=tri_color_intensity,
            intensitymode="cell",
            colorscale=pl_colorscale,
            showscale=False,
        )
    )
    return fig


def plot_frame(
    fig: go.Figure,
    origin: np.ndarray,
    rotation: np.ndarray,
    length_scale: float = 1,
    tip_scale: float = 1,
    width: float = 2,
) -> None:
    plot_arrow3(fig, origin, rotation[0], length_scale, tip_scale, "red", width)
    plot_arrow3(fig, origin, rotation[1], length_scale, tip_scale, "green", width)
    plot_arrow3(fig, origin, rotation[2], length_scale, tip_scale, "blue", width)


def plot_arrow3(
    fig: go.Figure,
    origin: np.ndarray = np.zeros(3),
    direction: np.ndarray = np.array([1, 0, 0]),
    length_scale: float = 1,
    tip_scale: float = 1,
    color: Union[str, list[str]] = "black",
    width: float = 2,
    scale: float = 3,
) -> go.Figure:
    """
    Add 3D arrow to a plotly figure

    Args:
        fig (go.Figure): plotly figure
        origin (np.ndarray): origin of the arrow
        direction (np.ndarray): direction of the arrow
        length_scale (float): scale factor for the length of the arrow
        tip_scale (float): scale factor for the tip of the arrow
        color (str): color of the arrow
        width (float): width of the arrow
    """
    assert fig is not None, "Please provide a plotly figure"
    if origin.ndim == 1:
        origin = origin[np.newaxis, :]
    if direction.ndim == 1:
        direction = direction[np.newaxis, :]
    if isinstance(color, str):
        color = [color] * origin.shape[0]

    origin = origin / 10**scale
    direction = direction / 10**scale
    for i in range(origin.shape[0]):
        fig.add_trace(
            go.Scatter3d(
                x=[origin[i, 0], origin[i, 0] + direction[i, 0] * length_scale],
                y=[origin[i, 1], origin[i, 1] + direction[i, 1] * length_scale],
                z=[origin[i, 2], origin[i, 2] + direction[i, 2] * length_scale],
                mode="lines",
                line=dict(color=color[i], width=width),
                showlegend=False,
            )
        )
        fig.add_trace(
            go.Cone(
                x=[origin[i, 0] + direction[i, 0] * length_scale],
                y=[origin[i, 1] + direction[i, 1] * length_scale],
                z=[origin[i, 2] + direction[i, 2] * length_scale],
                u=[direction[i, 0] * tip_scale],
                v=[direction[i, 1] * tip_scale],
                w=[direction[i, 2] * tip_scale],
                showscale=False,
                colorscale=[[0, color[i]], [1, color[i]]],
                showlegend=False,
            )
        )
    return fig


def plot_orbits(
    fig: go.Figure,
    rv: np.ndarray,
    t: int = None,
    marker_size: float = 4,
    scale: float = 3,
    color: Union[str, list[str]] = plotly_colors,
    **kwargs,
) -> go.Figure:
    rv = rv / 10**scale

    if rv.ndim == 2:
        rv = rv[np.newaxis, :, :]

    N_sat = rv.shape[0]
    for i in range(N_sat):
        fig.add_scatter3d(
            **dict(x=rv[i, :, 0], y=rv[i, :, 1], z=rv[i, :, 2]),
            mode="lines",
            line=dict(
                color=color[i % len(color)] if type(color) == list else color, width=3
            ),
            name="sat_lines",
            showlegend=False,
        )
    if t is not None:
        fig.add_scatter3d(
            **dict(x=rv[:, t, 0], y=rv[:, t, 1], z=rv[:, t, 2]),
            mode="markers",
            marker=dict(
                color=color, size=marker_size, line=dict(color=color, width=0.5)
            ),
            name="sat_markers",
            showlegend=False,
        )

    # Dummy trace
    c = 220
    axis_dict = dict(
        linecolor=f"rgb({c},{c},{c})",
        gridcolor=f"rgb({c},{c},{c})",
        linewidth=1.5,
        gridwidth=1.5,
        showline=False,
        mirror=True,
    )
    fig.update_layout(
        scene=dict(
            xaxis=dict(title="X [10<sup>3</sup> km]", **axis_dict),
            yaxis=dict(title="Y [10<sup>3</sup> km]", **axis_dict),
            zaxis=dict(title="Z [10<sup>3</sup> km]", **axis_dict),
            # xaxis=dict(title="X [10<sup>3</sup> km]", **axis_dict),
            # yaxis=dict(title="Y [10<sup>3</sup> km]", **axis_dict),
            # zaxis=dict(title="Z [10<sup>3</sup> km]", **axis_dict),
        ),
        font=dict(size=12, family="serif"),
        margin=dict(l=10, r=10, t=10, b=10),
        height=600,
        width=600,
        legend=dict(
            x=0.9,
            y=0.85,
            xanchor="right",
            yanchor="top",
            bgcolor="rgba(255, 255, 255, 0.5)",
        ),
    )
    fig.update_layout(scene_aspectmode="data")
    fig.update_layout(**kwargs)
    xticks = fig.layout.scene.xaxis.tickvals
    yticks = fig.layout.scene.yaxis.tickvals
    zticks = fig.layout.scene.zaxis.tickvals
    if xticks is not None and yticks is not None and zticks is not None:
        fig.update_layout(
            scene=dict(
                xaxis_ticktext=[f"{x/1e3:.0f}" for x in xticks],
                yaxis_ticktext=[f"{y/1e3:.0f}" for y in yticks],
                zaxis_ticktext=[f"{z/1e3:.0f}" for z in zticks],
            )
        )

    return fig


def image2zvals(
    img: np.ndarray, n_colors: int = 64, n_training_pixels: int = 800, rngs: int = 123
) -> tuple[np.ndarray, list]:
    """
    Image color quantization

    Args:
        img (np.ndarray): image array
        n_colors (int): number of colors for color quantization
        n_training_pixels (int): number of image pixels to fit a KMeans instance to them
        rngs (int): random seed
    Returns:
        tuple[np.ndarray, list]: z_values for the heatmap representation, and a plotly colorscale
    """
    # Image color quantization
    # img - np.ndarray of shape (m, n, 3) or (m, n, 4)
    # n_colors: int,  number of colors for color quantization
    # n_training_pixels: int, the number of image pixels to fit a KMeans instance to them
    # returns the array of z_values for the heatmap representation, and a plotly colorscale

    if img.ndim != 3:
        raise ValueError(
            f"Your image does not appear to  be a color image. It's shape is  {img.shape}"
        )
    rows, cols, d = img.shape
    if d < 3:
        raise ValueError(
            f"A color image should have the shape (m, n, d), d=3 or 4. Your  d = {d}"
        )

    range0 = img[:, :, 0].max() - img[:, :, 0].min()
    if range0 > 1:  # normalize the img values
        img = np.clip(img.astype(float) / 255, 0, 1)

    observations = img[:, :, :3].reshape(rows * cols, 3)
    training_pixels = shuffle(observations, random_state=rngs)[:n_training_pixels]
    model = KMeans(n_clusters=n_colors, random_state=rngs, n_init="auto").fit(
        training_pixels
    )

    codebook = model.cluster_centers_
    indices = model.predict(observations)
    z_vals = indices.astype(float) / (
        n_colors - 1
    )  # normalization (i.e. map indices to  [0,1])
    z_vals = z_vals.reshape(rows, cols)
    # define the Plotly colorscale with n_colors entries
    scale = np.linspace(0, 1, n_colors)
    colors = (codebook * 255).astype(np.uint8)
    pl_colorscale = [
        [sv, f"rgb{tuple(int(c) for c in color)}"] for sv, color in zip(scale, colors)
    ]

    # Reshape z_vals  to  img.shape[:2]
    return z_vals.reshape(rows, cols), pl_colorscale


def regular_tri(rows, cols):
    # define triangles for a np.meshgrid(np.linspace(a, b, cols), np.linspace(c,d, rows))
    triangles = []
    for i in range(rows - 1):
        for j in range(cols - 1):
            k = j + i * cols
            triangles.extend([[k, k + cols, k + 1 + cols], [k, k + 1 + cols, k + 1]])
    return np.array(triangles)


def mesh_data(img, n_colors=32, n_training_pixels=800):
    rows, cols, _ = img.shape
    z_data, pl_colorscale = image2zvals(
        img, n_colors=n_colors, n_training_pixels=n_training_pixels
    )
    triangles = regular_tri(rows, cols)
    I, J, K = triangles.T
    zc = z_data.flatten()[triangles]
    tri_color_intensity = np.array(
        [zc[k][2] if k % 2 else zc[k][1] for k in range(len(zc))]
    )
    return I, J, K, tri_color_intensity, pl_colorscale


def scatter(
    fig: go.Figure,
    xyz: np.ndarray,
    mode: str = "markers",
    marker_size: float = 4,
    color: Union[str, list[str]] = plotly_colors,
    scale: float = 3,
    **kwargs,
) -> go.Figure:
    """
    Create a 3D scatter plot

    Args:
        x (np.ndarray): x coordinates
        y (np.ndarray): y coordinates
        z (np.ndarray): z coordinates
        mode (str): plot mode
        marker_size (float): marker size
        color (str): marker color
    """
    if xyz.ndim == 2:
        xyz = xyz[np.newaxis, :, :]
    xyz = xyz / 10**scale
    N = xyz.shape[0]
    for i in range(N):
        fig.add_scatter3d(
            **dict(x=xyz[i, :, 0], y=xyz[i, :, 1], z=xyz[i, :, 2]),
            mode=mode,
            marker=dict(
                color=color[i % len(color)] if type(color) == list else color,
                size=marker_size,
            ),
            showlegend=False,
        )
    return fig

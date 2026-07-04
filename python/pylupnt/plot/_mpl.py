import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "cmr10"
plt.rcParams["font.sans-serif"] = "cmss10"
plt.rcParams["font.monospace"] = "cmtt10"
plt.rcParams["axes.formatter.use_mathtext"] = True
from PIL import Image
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

# Removed utils import - functionality moved to other modules

from .. import _pylupnt as _pnt

###patch start###
from mpl_toolkits.mplot3d.axis3d import Axis

if not hasattr(Axis, "_get_coord_info_old"):

    def _get_coord_info_new(self, renderer):
        mins, maxs, centers, deltas, tc, highs = self._get_coord_info_old(renderer)
        mins += deltas / 4
        maxs -= deltas / 4
        return mins, maxs, centers, deltas, tc, highs

    Axis._get_coord_info_old = Axis._get_coord_info
    Axis._get_coord_info = _get_coord_info_new
###patch end###

COLORS = list(mcolors.TABLEAU_COLORS.keys())

# plt.rc("text", usetex=True)
plt.rc("font", family="serif")

plot_data = {
    _pnt.EARTH: {
        "filename": "earth_surface.jpg",
        "RE": 6378.137,
        "lim": 25e3,
        "brightness": 3,
    },
    _pnt.MOON: {
        "filename": "moon_surface.jpeg",
        "RE": 1737.1,
        "lim": 10e3,
        "brightness": 1.5,
    },
}


class Plot3D:
    fig: plt.Figure
    ax: plt.Axes
    name: str
    scatters: list
    plots: list

    def __init__(self, azim=-60, elev=30, figsize=(10, 10)):
        self.fig = plt.figure(figsize=figsize)
        self.ax = self.fig.add_subplot(
            111, projection="3d", proj_type="ortho", computed_zorder=False
        )
        self.ax.view_init(azim=azim, elev=elev)
        self.azim = self.ax.azim
        self.elev = self.ax.elev
        self.scatters = []
        self.plots = []

    def plot_surface(self, name, offset=np.array([0, 0, 0]), adjust_axis=True, limit=None, scale=6):
        self.name = name

        img_data = plot_data[name]
        from pylupnt.core.pylupnt_utils import LUPNT_DATA_PATH

        image_file = os.path.join(LUPNT_DATA_PATH, "topo", img_data["filename"])
        img = Image.open(image_file)
        img = np.array(img.resize([int(d / scale) for d in img.size])) / 256.0
        img = np.minimum(img * img_data["brightness"], 1)
        # img = np.roll(img, int(img.shape[0] * -180 / 360), axis=1)
        img = img[::-1, :]
        img = img[:, ::-1]

        lons = np.linspace(-180, 180, img.shape[1]) * np.pi / 180
        lats = np.linspace(-90, 90, img.shape[0])[::-1] * np.pi / 180

        RE = plot_data[name]["RE"]
        x = np.outer(np.cos(lons), np.cos(lats)).T * RE + offset[0]
        y = np.outer(np.sin(lons), np.cos(lats)).T * RE + offset[1]
        z = np.outer(np.ones(np.size(lons)), np.sin(lats)).T * RE + offset[2]

        self.ax.plot_surface(x, y, z, rstride=4, cstride=4, facecolors=img, zorder=-1)

        if adjust_axis:
            self.ax.axis("scaled")

            lim = plot_data[name]["lim"] if limit is None else limit
            ticks = (-lim, -lim / 2, 0, lim / 2, lim)
            lims = (-lim, lim)
            self.set_ticks(ticks, ticks, ticks)
            self.set_tickpad(0)
            self.set_labels("X [m]", "Y [m]", "Z [m]")
            self.set_labelpad(0, 0, 0)
            self.set_lims(lims, lims, lims, equal=True)

        # self.fig.canvas.mpl_connect("motion_notify_event", self.rotate)

    def rotate(self, event):
        if event.inaxes == self.ax:
            self.plot_visible(self.ax.azim, self.ax.elev)

    def plot_visible(self, azimuth, elev):
        # transform viewing angle to normal vector in data coordinates
        a = azimuth * np.pi / 180.0 - np.pi
        e = elev * np.pi / 180.0 - np.pi / 2.0
        for points, data in zip(self.points, self.data):
            X = np.array([np.sin(e) * np.cos(a), np.sin(e) * np.sin(a), np.cos(e)])
            # concatenate coordinates
            Z = data
            # calculate dot product
            # the points where this is positive are to be shown
            proj = np.dot(Z, X)
            RE = plot_data[self.name]["RE"]
            cond = np.logical_or(
                proj >= 0,
                np.linalg.norm(Z - X * proj.reshape(-1, 1), axis=1) > RE,
            )
            # filter points by the above condition
            x_c = data[cond, 0]
            y_c = data[cond, 1]
            z_c = data[cond, 2]
            # set the new data points
            points.set_data(x_c, y_c)
            points.set_3d_properties(z_c, zdir="z")
        self.fig.canvas.draw_idle()

    def check_occultation(self, data):
        a = self.ax.azim * np.pi / 180.0 - np.pi
        e = self.ax.elev * np.pi / 180.0 - np.pi / 2.0
        view = np.array([np.sin(e) * np.cos(a), np.sin(e) * np.sin(a), np.cos(e)])
        proj = np.dot(data, view)
        RE = plot_data[self.name]["RE"]
        cond = np.logical_or(
            proj <= 0,
            np.linalg.norm(data - view * proj.reshape(-1, 1), axis=1) >= RE,
        )
        alphas = (proj - np.min(proj)) / (np.max(proj) - np.min(proj))
        return cond, alphas

    def scatter(self, data, mask=False, *args, **kwargs):
        """
        Plot Cartesian coordinates
        """
        if mask:
            cond, _ = self.check_occultation(data)
            data[np.logical_not(cond), :] = [np.nan, np.nan, np.nan]
        self.ax.scatter(data[:, 0], data[:, 1], data[:, 2], *args, zorder=1, **kwargs)

    def plot(self, data, *args, mask=False, **kwargs):
        """
        Plot Cartesian coordinates
        """
        if len(data.shape) == 3:
            n_data = data.shape[0]
            for i in range(n_data):
                if mask:
                    cond, _ = self.check_occultation(data[i])
                    data[i, np.logical_not(cond), :] = [np.nan, np.nan, np.nan]
                self.ax.plot(
                    data[i, :, 0],
                    data[i, :, 1],
                    data[i, :, 2],
                    *args,
                    zorder=0,
                    **kwargs,
                )
        else:
            if mask:
                cond, _ = self.check_occultation(data)
                data[np.logical_not(cond), :] = [np.nan, np.nan, np.nan]
            self.ax.plot(data[:, 0], data[:, 1], data[:, 2], *args, zorder=0, **kwargs)

    def set_labels(self, x: str, y: str, z: str) -> None:
        self.ax.set_xlabel(x)
        self.ax.set_ylabel(y)
        self.ax.set_zlabel(z)

    def set_pane_color(self, color: tuple) -> None:
        self.ax.xaxis.set_pane_color(color)
        self.ax.yaxis.set_pane_color(color)
        self.ax.zaxis.set_pane_color(color)

    def set_labelpad(self, padx: int, pady: int, padz: int) -> None:
        self.ax.xaxis.labelpad = padx
        self.ax.yaxis.labelpad = pady
        self.ax.zaxis.labelpad = padz

    def set_tickpad(self, pad: int) -> None:
        self.ax.tick_params(axis="both", which="major", pad=0)

    def set_tick_multiplier(self, factor: int) -> None:
        self.ax.set_xticklabels([f"{int(x * factor)}" for x in self.ax.get_xticks()])
        self.ax.set_yticklabels([f"{int(y * factor)}" for y in self.ax.get_yticks()])
        self.ax.set_zticklabels([f"{int(z * factor)}" for z in self.ax.get_zticks()])

    def set_ticks(self, x: list, y: list, z: list) -> None:
        self.ax.set_xticks(x)
        self.ax.set_yticks(y)
        self.ax.set_zticks(z)

    def set_lims(self, xlims: tuple, ylims: tuple, zlims: tuple, equal=True) -> None:

        self.ax.set_xlim(xlims)
        self.ax.set_ylim(ylims)
        self.ax.set_zlim(zlims)

        if equal:
            self.ax.set_box_aspect([xlims[1] - xlims[0], ylims[1] - ylims[0], zlims[1] - zlims[0]])


def plot_antenna_gain_patter_2D(
    ax: plt.Axes = None,
    antenna: _pnt.Antenna = None,
    phi: np.ndarray = None,
    theta: np.ndarray = None,
    show_max: bool = True,
    name: str = None,
):
    if ax is not None:
        plt.sca(ax)
    if phi is None:
        phi = np.linspace(-180, 180, 500)  # [deg]
    if theta is None:
        theta = np.linspace(0, 90, 4)  # [deg]
    for az in theta:
        gain = antenna.compute_gain(az * _pnt.RAD, phi * _pnt.RAD)
        plt.plot(
            _pnt.DEG * _pnt.wrap_to_pi(phi * _pnt.RAD),
            gain,
            label=f"$\\varphi = {az:.0f}^\\circ$",
        )
    if name is None:
        plt.title(rf"{antenna.name}")
    else:
        plt.title(rf"{name}")
    plt.xlabel("Boresite Angle $\\theta$ [deg]")
    plt.ylabel("Gain $G$ [dB]")
    plt.text(
        0.98,
        0.95,
        f"$G_{{\\max}} = {antenna.get_gain_matrix().max():.2f}$ dB",
        transform=plt.gca().transAxes,
        ha="right",
        va="top",
    )
    plt.grid()
    # Legend with titl
    if len(theta) > 1:
        plt.legend(title="Azimuth $\\theta$")


def plot_2d_ellipse(mu, cov, n_std=3.0, facecolor="none", **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The Axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensional dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse(
        (0, 0),
        width=ell_radius_x * 2,
        height=ell_radius_y * 2,
        facecolor=facecolor,
        **kwargs,
    )

    # Calculating the standard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = mu[0]

    # calculating the standard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = mu[1]

    transf = transforms.Affine2D().rotate_deg(45).scale(scale_x, scale_y).translate(mean_x, mean_y)

    ax = plt.gca()
    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)

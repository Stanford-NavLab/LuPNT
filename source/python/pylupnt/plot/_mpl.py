import os
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from .. import utils

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

    def plot_surface(
        self, name, offset=np.array([0, 0, 0]), adjust_axis=True, limit=None, scale=3
    ):
        self.name = name

        img_data = plot_data[name]
        image_file = os.path.join(utils.LUPNT_DATA_PATH, "topo", img_data["filename"])
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
            self.set_labels("X [km]", "Y [km]", "Z [km]")
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
            self.ax.set_box_aspect(
                [xlims[1] - xlims[0], ylims[1] - ylims[0], zlims[1] - zlims[0]]
            )

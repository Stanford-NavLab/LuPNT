import os
from PIL import Image
import numpy as np
from . import utils as u
import matplotlib.pyplot as plt

import matplotlib.colors as mcolors

COLORS = list(mcolors.TABLEAU_COLORS.keys())

plt.rc("text", usetex=True)
plt.rc("font", family="serif")
os.chdir(os.path.dirname(os.path.abspath(__file__)))

plot_data = {
    "EARTH": {
        "filename": "earth_surface.jpg",
        "RE": 6378.137,
        "lim": 25e3,
        "scale": 1,
        "brightness": 3,
    },
    "MOON": {
        "filename": "moon_surface.jpeg",
        "RE": 1737.1,
        "lim": 10e3,
        "scale": 15,
        "brightness": 1.5,
    },
}

images_dict = {}
for k, v in plot_data.items():
    image_file = os.path.join(u.basepath, "topo", v["filename"])
    img = Image.open(image_file)
    img = np.array(img.resize([int(d / v["scale"]) for d in img.size])) / 256.0
    img = np.minimum(img * v["brightness"], 1)
    # img = np.roll(img, int(img.shape[0] * -180 / 360), axis=1)
    img = img[::-1, :]
    img = img[:, ::-1]
    plot_data[k]["img"] = img


class Plot3D:
    fig = None
    ax = None
    points = []
    data = []
    name = None

    def __init__(self, azim=-60, elev=30):
        self.fig = plt.figure(figsize=(10, 10))
        self.ax = self.fig.add_subplot(111, projection="3d", computed_zorder=False)
        self.ax.view_init(azim=azim, elev=elev)
        self.azim = self.ax.azim
        self.elev = self.ax.elev

    def plot_surface(
        self, name, offset=np.array([0, 0, 0]), adjust_axis=True, limit=None
    ):
        self.name = name
        img = plot_data[name]["img"]
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
            self.ax.set_xticks(ticks)
            self.ax.set_yticks(ticks)
            self.ax.set_zticks(ticks)

            self.ax.tick_params(axis="z", which="major", pad=10)
            self.ax.set_xlabel("X [km]")
            self.ax.set_zlabel("Z [km]", labelpad=10)
            self.ax.set_ylabel("Y [km]", labelpad=5)

            self.ax.set_xlim([-lim, lim])
            self.ax.set_ylim([-lim, lim])
            self.ax.set_zlim([-lim, lim])

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
            RE = plot_data[name]["RE"]
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

    def scatter(self, data, *args, **kwargs):
        """
        Plot Cartesian coordinates
        """
        cond, _ = self.check_occultation(data)
        data[np.logical_not(cond), :] = [np.nan, np.nan, np.nan]
        self.ax.scatter(data[:, 0], data[:, 1], data[:, 2], *args, zorder=1, **kwargs)
        # (points,) = self.ax.plot([], [], [], "ro", zorder=0)
        # self.points.append(points)
        # self.data.append(data)

    def plot(self, data, *args, **kwargs):
        """
        Plot Cartesian coordinates
        """
        if len(data.shape) == 3:
            n_data = data.shape[0]
            for i in range(n_data):
                cond, _ = self.check_occultation(data[i])
                data[i, np.logical_not(cond), :] = [np.nan, np.nan, np.nan]
                self.ax.plot(
                    data[i, :, 0],
                    data[i, :, 1],
                    data[i, :, 2],
                    *args,
                    zorder=0,
                    **kwargs
                )
        else:
            cond, _ = self.check_occultation(data)
            data[np.logical_not(cond), :] = [np.nan, np.nan, np.nan]
            self.ax.plot(data[:, 0], data[:, 1], data[:, 2], *args, zorder=0, **kwargs)
            # self.ax.plot(
            #     data[i, cond, 0], data[i, cond, 1], data[i, cond, 2], *args, **kwargs
            # )
        # for i in range(n_data):
        #     (points,) = self.ax.plot([], [], [], zorder=0)
        #     self.points.append(points)
        #     self.data.append(data[i])

    def label_axis(self, x="X [km]", y="Y [km]", z="Z [km]"):
        self.ax.set_xlabel(x)
        self.ax.set_ylabel(y)
        self.ax.set_zlabel(z)

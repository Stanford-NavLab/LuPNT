import pylupnt as pnt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

try:
    from .gve import get_rtn_matrix
    from .orbit_utils import wrapToPi
except ImportError:
    from gve import get_rtn_matrix
    from orbit_utils import wrapToPi


class Ephemeris:

    def __init__(self, body=pnt.MOON):
        self.body = body
        if body == pnt.MOON:
            self.GM = pnt.GM_MOON
            self.omega_b = 2.6616957278e-6
            self.inertial = pnt.MOON_CI
            self.bodyfixed = pnt.MOON_PA
            self.T_sidereal = 27.321661 * 86400  # seconds
        elif body == pnt.EARTH:
            self.GM = pnt.GM_EARTH
            self.omega_b = 7.292115e-5
            self.inertial = pnt.GCRF
            self.bodyfixed = pnt.ITRF
            self.T_sidereal = 86164.0905  # seconds
        else:
            # error
            print("Unsupported central body")
            return

    def set_omega_b(self, t_tai):
        if self.body == pnt.MOON:
            angles = pnt.get_lunar_orientation_angles(t_tai)
            self.omega_b = angles[5]
            print(f"set omega_b to: {self.omega_b}")
        else:
            # error
            print("Unsupported central body")
            return

    def ephem2dict(self, ephem, scale=None):
        dict = {}

        for i, x in enumerate(ephem):
            key = self.keys_list[i]
            if scale is None:
                dict[key] = x
            else:
                if scale[i] == 0:
                    dict[key] = x
                else:
                    # quantize based on the scale
                    scale_i = int(scale[i])
                    # print(f"Key: {key}, Scale: {scale_i}")
                    x_int = np.round(x * (2**scale_i)).astype(int)
                    dict[key] = x_int.astype(float) / (1 << scale_i)

        return dict

    def dict2ephem(self, ephem_dict):
        ephem = np.zeros(self.n_params)

        for key, value in ephem_dict.items():
            if key in self.idx_keys.keys():
                idx = self.idx_keys[key]
                ephem[idx] = value
            else:
                print("Key not found in ephem_dict: ", key)
                return None

        return ephem

    def get_index(self, key):
        if key in self.idx_keys.keys():
            return self.idx_keys[key]
        else:
            print("Key not found in ephem_dict: ", key)
            return None

    def replace_ephem(self, ephem_base, replace_keys, replace_values):
        """
        Replace the parameters in the ephemeris

        ephem_base: np.array of base ephemeris
        keys: list of keys to replace
        values: np.array of values to replace
        """
        ephem_keys = self.ephem2dict(ephem_base).keys()
        ephem_array = np.zeros(self.n_params)

        for key in ephem_keys:
            idx = self.get_index(key)
            if key in replace_keys:  # replace the value
                ephem_array[idx] = replace_values[replace_keys.index(key)]
            else:  # keep the original value
                ephem_array[idx] = ephem_base[idx]

        return ephem_array

    def plot_fit_error(
        self,
        t_data,
        rvf_data,
        rvfw_data,
        ephem,
        plot_diff=False,
        plot_init=False,
        print_stats=True,
        use_grad_for_velfit=False,
        axes=None,
        title_txt=None,
        use_rtn=True,
        moving_average_size=60,
        ylim_pos=10,
        ylim_vel=20,
        t_data_fit=None,
        in_kms=False,
        plot_azel=False,
        figname_error=None,
        figname_azel=None,
        azel_xlim=90,
        mask_southpole=False,
        plot_prctile=False,
        plot_azel_legend=False,
        tworows_azel=False,
        plot_azel_range_doppler=False,
        figname_azel_range_doppler=None,
    ):
        """
        Plot the fit of the ephemeris

        Args:
            t_data: array of times
            rvf_data: array of position and velocity vectors in fixed frame
            ephem: ephemeris parameters
            axes: axes to plot the data
        """
        ms = 1  # marker size
        lw = 1  # line width

        t_ref, t_fit, x0 = self.init_guess(t_data, rvf_data)

        # trim off both the edge points
        print("Eval Time Interval [min]: ", (t_data[-1] - t_data[0]) / 60)

        if axes is None:
            fig, axes = plt.subplots(2, 4, figsize=(16, 6))

        # compute fitted trajectory and errors --------------------------------
        lent = rvf_data.shape[0]
        fit_y0 = self.ephem2cart(t_data, x0, compute_velocity=True, return_params=False)
        fit_y = self.ephem2cart(t_data, ephem, compute_velocity=True, return_params=False)

        # use gradient to compute velocity
        if use_grad_for_velfit:
            for i in range(3):
                fit_y0[:, 3 + i] = np.gradient(fit_y0[:, i], t_data, axis=0)
                fit_y[:, 3 + i] = np.gradient(fit_y[:, i], t_data, axis=0)

        # plot errors ---------------------------------------------------------
        diff_y0 = fit_y0 - rvfw_data
        diff_y = fit_y - rvfw_data

        if in_kms:
            diff_y0[:, :3] = diff_y0[:, :3] / 1000  # convert to km
            diff_y[:, :3] = diff_y[:, :3] / 1000

            diff_y0[:, 3:] = diff_y0[:, 3:]  # convert to m/s
            diff_y[:, 3:] = diff_y[:, 3:]

        # convert to rtn
        vel_rtn_true = np.zeros((lent, 3))
        if use_rtn:
            print("Convert to rtn")
            for ti in range(lent):
                M0 = get_rtn_matrix(rvf_data[ti])
                diff_y0[ti, :3] = np.dot(M0, diff_y0[ti, :3])
                diff_y[ti, :3] = np.dot(M0, diff_y[ti, :3])
                diff_y0[ti, 3:] = np.dot(M0, diff_y0[ti, 3:])
                diff_y[ti, 3:] = np.dot(M0, diff_y[ti, 3:])
                vel_rtn_true[ti, :] = np.dot(M0, rvfw_data[ti, 3:])

            # if use_Hz:
            #     for i in range(3):
            #         fc = 2492e6  # Hz (carrier frequency)
            #         diff_y0[:, 3 + i] = - diff_y0[:, 3 + i] * fc / pnt.C
            #         diff_y[:, 3 + i] = - diff_y[:, 3 + i] * fc / pnt.C

        t_data = (t_data - t_data[0]) / 60
        if t_data_fit is not None:
            t_data_fit = (t_data_fit - t_data_fit[0]) / 60

        # Compute Elevation at South Pole
        user_pos = np.array([0, 0, -pnt.R_MOON]).reshape(1, 3)
        user_pos_tile = np.tile(user_pos, (rvf_data.shape[0], 1))
        r_u2s = np.linalg.norm(rvf_data[:, :3] - user_pos_tile, axis=1)
        r_u2m = -user_pos_tile
        theta_mus = np.arccos(
            np.sum(rvf_data[:, :3] * r_u2m, axis=1)
            / (np.linalg.norm(rvf_data[:, :3], axis=1) * np.linalg.norm(r_u2m, axis=1))
        )
        idx_visible = theta_mus > np.deg2rad(95)  # above horizon

        if plot_diff:
            if in_kms:
                if use_rtn:
                    labels = [
                        "Radial Pos Err [km]",
                        "Tangential Pos Err [km]",
                        "Normal Pos Err [km]",
                        "3D Pos Err [km]",
                    ]
                else:
                    labels = [
                        r"$\Delta x$ [km]",
                        r"$\Delta y$ [km]",
                        r"$\Delta z$ [km]",
                        r"$\Delta pos$ [km]",
                    ]
            else:
                if use_rtn:
                    labels = [
                        r"$\Delta R$ [m]",
                        r"$\Delta T$ [m]",
                        r"$\Delta N$ [m]",
                        r"$\Delta pos$ [m]",
                    ]
                else:
                    labels = [
                        r"$\Delta x$ [m]",
                        r"$\Delta y$ [m]",
                        r"$\Delta z$ [m]",
                        r"$\Delta pos$ [m]",
                    ]
        else:
            labels = ["$x$ [km]", "$y$ [km]", "$z$ [km]", r"$\Delta pos$ [m]"]

        masize = moving_average_size

        # positions ------------------------------------------------------------------------------------------
        for i, ax in enumerate(axes[0, :]):
            if i == 3:
                if plot_init:
                    ax.plot(
                        t_data,
                        np.linalg.norm(diff_y0[:, :3], axis=1),
                        "r--",
                        markersize=ms,
                        linewidth=lw,
                        label="True - Initial",
                    )
                # moving average
                maplot = np.convolve(
                    np.linalg.norm(diff_y[:, :3], axis=1), np.ones(masize) / masize, mode="valid"
                )
                # Plot raw data
                if mask_southpole:
                    ax.plot(
                        t_data[idx_visible],
                        np.linalg.norm(diff_y[idx_visible, :3], axis=1),
                        "bo",
                        markersize=ms,
                        linewidth=lw,
                        label="Elev over 5 deg",
                    )
                    ax.plot(
                        t_data[~idx_visible],
                        np.linalg.norm(diff_y[~idx_visible, :3], axis=1),
                        "ro",
                        markersize=ms,
                        linewidth=lw,
                        label="Elev below 5 deg",
                    )
                else:
                    ax.plot(
                        t_data,
                        np.linalg.norm(diff_y[:, :3], axis=1),
                        "bo",
                        markersize=ms,
                        linewidth=lw,
                        label="Fit - True",
                    )

                if ylim_pos is not None:
                    ax.set_ylim([0, ylim_pos])
                # ax.plot(t_data[masize:lent-masize], 1000 * maplot[masize:lent-masize], 'mo-', markersize=ms,linewidth=2*lw, label='Fit - True (MA)')
            else:  # x, y, z
                if plot_diff:
                    if plot_init:
                        ax.plot(
                            t_data,
                            diff_y0[:, i],
                            "r--",
                            markersize=ms,
                            linewidth=lw,
                            label="True - Initial",
                        )
                    maplot = np.convolve(diff_y[:, i], np.ones(masize) / masize, mode="valid")

                    if mask_southpole:
                        ax.plot(
                            t_data[idx_visible],
                            diff_y[idx_visible, i],
                            "bo",
                            markersize=ms,
                            linewidth=lw,
                            label="Elev over 5 deg",
                        )
                        ax.plot(
                            t_data[~idx_visible],
                            diff_y[~idx_visible, i],
                            "ro",
                            markersize=ms,
                            linewidth=lw,
                            label="Elev below 5 deg",
                        )
                    else:
                        ax.plot(
                            t_data,
                            diff_y[:, i],
                            "bo",
                            markersize=ms,
                            linewidth=lw,
                            label="Fit - True",
                        )

                    if ylim_pos is not None:
                        ax.set_ylim([-ylim_pos, ylim_pos])
                    # ax.plot(t_data[masize:lent-masize], 1000 * maplot[masize:lent-masize], 'mo-', markersize=ms,linewidth=2*lw, label='Fit - True (MA)')
                else:
                    ax.plot(
                        t_data, rvfw_data[:, i], "k--", markersize=ms, linewidth=lw, label="True"
                    )
                    if plot_init:
                        ax.plot(
                            t_data,
                            fit_y0[:, i],
                            "r--",
                            markersize=ms,
                            linewidth=lw,
                            label="Initial",
                        )
                    ax.plot(t_data, fit_y[:, i], "b-", markersize=ms, linewidth=lw, label="Fit")
            ax.set_xlabel("Time [min]")
            ax.set_ylabel(labels[i])
            ax.legend()
            ax.grid(True)

        # velocities
        if use_grad_for_velfit:
            # remove first and last point
            t_data = t_data[1:-1]
            diff_y0 = diff_y0[1:-1, :]
            diff_y = diff_y[1:-1, :]
            rvf_data = rvf_data[1:-1, :]
            rvfw_data = rvfw_data[1:-1, :]
            fit_y0 = fit_y0[1:-1, :]
            fit_y = fit_y[1:-1, :]
            theta_mus = theta_mus[1:-1]
            idx_visible = idx_visible[1:-1]
            user_pos_tile = user_pos_tile[1:-1, :]

        if plot_diff:
            if use_rtn:
                labels = [
                    "Radial Vel Err [m/s]",
                    "Tangential Vel Err [m/s]",
                    "Normal Vel Err [m/s]",
                    "3D Vel Err [m/s]",
                ]
            else:
                labels = [
                    r"$\Delta V_x$ [mm/s]",
                    r"$\Delta V_y$ [mm/s]",
                    r"$\Delta V_z$ [mm/s]",
                    r"$\Delta vel$ [mm/s]",
                ]
        else:
            labels = ["$V_R$ [km/s]", "$V_T$ [km/s]", "$V_N$ [km/s]", r"$\Delta vel$ [mm/s]"]

        # Velocity plots
        for i, ax in enumerate(axes[1, :]):
            if i == 3:
                if plot_init:
                    ax.plot(
                        t_data,
                        np.linalg.norm(diff_y0[:, 3:], axis=1),
                        "r--",
                        markersize=ms,
                        linewidth=lw,
                        label="True - Initial",
                    )
                # moving average
                maplot = np.convolve(
                    np.linalg.norm(diff_y[:, 3:], axis=1), np.ones(masize) / masize, mode="valid"
                )
                if mask_southpole:
                    ax.plot(
                        t_data[idx_visible],
                        np.linalg.norm(diff_y[idx_visible, 3:], axis=1),
                        "bo",
                        markersize=ms,
                        linewidth=lw,
                        label="Visible",
                    )
                    ax.plot(
                        t_data[~idx_visible],
                        np.linalg.norm(diff_y[~idx_visible, 3:], axis=1),
                        "ro",
                        markersize=ms,
                        linewidth=lw,
                        label="Below Horizon",
                    )
                else:
                    ax.plot(
                        t_data,
                        np.linalg.norm(diff_y[:, 3:], axis=1),
                        "bo",
                        markersize=ms,
                        linewidth=lw,
                        label="Fit - True",
                    )
                # ax.plot(t_data[masize:lent-masize], 1e6 * maplot[masize:lent-masize], 'mo-', markersize=ms,linewidth=2*lw, label='Fit - True')
                if ylim_vel is not None:
                    ax.set_ylim([0, ylim_vel])
            else:  # x, y, z
                if plot_diff:
                    if plot_init:
                        ax.plot(
                            t_data,
                            diff_y0[:, i + 3],
                            "r--",
                            markersize=ms,
                            linewidth=lw,
                            label="True - Initial",
                        )
                    maplot = np.convolve(diff_y[:, i + 3], np.ones(masize) / masize, mode="valid")

                    if mask_southpole:
                        ax.plot(
                            t_data[idx_visible],
                            diff_y[idx_visible, i + 3],
                            "bo",
                            markersize=ms,
                            linewidth=lw,
                            label="Visible",
                        )
                        ax.plot(
                            t_data[~idx_visible],
                            diff_y[~idx_visible, i + 3],
                            "ro",
                            markersize=ms,
                            linewidth=lw,
                            label="Below Horizon",
                        )
                    else:
                        ax.plot(
                            t_data,
                            diff_y[:, i + 3],
                            "bo",
                            markersize=ms,
                            linewidth=lw,
                            label="Fit - True",
                        )
                    # ax.plot(t_data[masize:lent-masize], 1e6 * maplot[masize:lent-masize], 'mo-', markersize=ms,linewidth=2*lw, label='Fit - True (MA)')
                    if ylim_vel is not None:
                        ax.set_ylim([-ylim_vel, ylim_vel])
                else:
                    ax.plot(
                        t_data, rvfw_data[:, i + 3], "k-", markersize=ms, linewidth=lw, label="True"
                    )
                    if plot_init:
                        ax.plot(
                            t_data,
                            fit_y0[:, i + 3],
                            "r--",
                            markersize=ms,
                            linewidth=lw,
                            label="Initial",
                        )
                    ax.plot(t_data, fit_y[:, i + 3], "b-", markersize=ms, linewidth=lw, label="Fit")
            ax.set_xlabel("Time [min]")
            ax.set_ylabel(labels[i])
            # ax.set_xlim([3, 4])
            ax.legend()
            ax.grid(True)

        if title_txt is not None:
            plt.suptitle(title_txt)

        if t_data_fit is not None:
            idx_t_data_fit = np.where((t_data_fit >= t_data[0]) & (t_data_fit <= t_data[-1]))[0]
            pos_norm = np.linalg.norm(diff_y[:, :3], axis=1)
            vel_norm = np.linalg.norm(diff_y[:, 3:], axis=1)
            interp_tdata_fit_pos = interp1d(t_data, pos_norm, kind="cubic")(
                t_data_fit[idx_t_data_fit]
            )
            axes[3, 0].plot(
                t_data_fit[idx_t_data_fit], interp_tdata_fit_pos, "ro", markersize=ms * 3
            )
            interp_tdata_fit_vel = interp1d(t_data, vel_norm, kind="cubic")(
                t_data_fit[idx_t_data_fit]
            )
            axes[3, 1].plot(
                t_data_fit[idx_t_data_fit], interp_tdata_fit_vel, "ro", markersize=ms * 3
            )

        plt.tight_layout()
        if figname_error is not None:
            plt.savefig(figname_error, dpi=300)
        plt.show()

        # ----------------------------------------------------------------------
        if plot_azel:
            if tworows_azel:
                fig, axes = plt.subplots(2, 2, figsize=(10, 6))
                ax = axes.flatten()
            else:
                fig, ax = plt.subplots(1, 4, figsize=(16, 3))

            r_vect = rvfw_data[:, :3]
            lat = np.arcsin(r_vect[:, 2] / np.linalg.norm(r_vect, axis=1))
            lon = np.arctan2(r_vect[:, 1], r_vect[:, 0])

            # plot fitted az/el
            r_fit_vect = fit_y[:, :3]
            lat_fit = np.arcsin(r_fit_vect[:, 2] / np.linalg.norm(r_fit_vect, axis=1))
            lon_fit = np.arctan2(r_fit_vect[:, 1], r_fit_vect[:, 0])

            # ax1: latitude
            lat_diff = np.abs(np.rad2deg(wrapToPi(lat - lat_fit)))
            lon_diff = np.abs(np.rad2deg(wrapToPi(lon - lon_fit)))
            med_lat = np.median(lat_diff)
            med_lon = np.median(lon_diff)
            p95_lat = np.percentile(lat_diff, 95)
            p95_lon = np.percentile(lon_diff, 95)
            p99_lat = np.percentile(lat_diff, 99)
            p99_lon = np.percentile(lon_diff, 99)

            # ax1: latitude
            if mask_southpole:
                ax[0].plot(
                    t_data[idx_visible] / 24 / 60,
                    lat_diff[idx_visible],
                    "bo",
                    markersize=2,
                    label="Latitude Err (Elev over 5 deg)",
                )
                ax[0].plot(
                    t_data[~idx_visible] / 24 / 60,
                    lat_diff[~idx_visible],
                    "ro",
                    markersize=2,
                    label="Latitude Err (Elev below 5 deg)",
                )
            else:
                ax[0].plot(
                    t_data / 24 / 60, lat_diff, "bo", markersize=2, label="True - Fit (latitude)"
                )
            ax[0].set_xlabel("Time [days]")
            ax[0].set_ylabel("Error [deg]")
            ax[0].legend()
            ax[0].grid(True)
            # ax[0].set_yscale('log')

            # ax2: longitude
            if mask_southpole:
                ax[1].plot(
                    t_data[idx_visible] / 24 / 60,
                    lon_diff[idx_visible],
                    "bo",
                    markersize=2,
                    label="Longitude Err (Elev over 5 deg)",
                )
                ax[1].plot(
                    t_data[~idx_visible] / 24 / 60,
                    lon_diff[~idx_visible],
                    "ro",
                    markersize=2,
                    label="Longitude Err (Elev below 5 deg)",
                )
            else:
                ax[1].plot(
                    t_data / 24 / 60, lon_diff, "bo", markersize=2, label="True - Fit (longitude)"
                )
            ax[1].set_xlabel("Time [days]")
            ax[1].set_ylabel("Error [deg]")
            ax[1].legend()
            ax[1].grid(True)
            # ax[1].set_yscale('log')

            # ax3: histogram
            bins = np.linspace(0, 90, 360)
            ax[2].hist(
                lat_diff,
                bins=bins,
                alpha=0.7,
                label="Latitude Error (Median: {:.2f})".format(med_lat),
                density=True,
            )
            ax[2].hist(
                lon_diff,
                bins=bins,
                alpha=0.7,
                label="Longitude Error (Median: {:.2f})".format(med_lon),
                density=True,
            )
            if plot_prctile:
                ax[2].axvline(
                    p95_lat,
                    color="r",
                    linestyle="--",
                    label="Latitude 95%: {:.2f} deg".format(p95_lat),
                )
                ax[2].axvline(
                    p95_lon,
                    color="b",
                    linestyle="--",
                    label="Longitude 95%: {:.2f} deg".format(p95_lon),
                )
                ax[2].axvline(
                    p99_lat,
                    color="r",
                    linestyle=":",
                    label="Latitude 99%: {:.2f} deg".format(p99_lat),
                )
                ax[2].axvline(
                    p99_lon,
                    color="b",
                    linestyle=":",
                    label="Longitude 99%: {:.2f} deg".format(p99_lon),
                )
            ax[2].set_xlabel("Error [deg]")
            ax[2].set_ylabel("Density")
            ax[2].set_xlim([0, azel_xlim])
            if plot_azel_legend:
                ax[2].legend()
            # ax[1].set_yscale('log')
            ax[2].grid(True)
            ax[2].set_title("All data")

            # ax4: histogram (visible only)
            ax[3].hist(
                lat_diff[idx_visible],
                bins=bins,
                alpha=0.7,
                label="Latitude Error (Median: {:.2f})".format(np.median(lat_diff[idx_visible])),
                density=True,
            )
            ax[3].hist(
                lon_diff[idx_visible],
                bins=bins,
                alpha=0.7,
                label="Longitude Error (Median: {:.2f})".format(np.median(lon_diff[idx_visible])),
                density=True,
            )
            if plot_prctile:
                ax[3].axvline(
                    np.percentile(lat_diff[idx_visible], 95),
                    color="r",
                    linestyle="--",
                    label="Latitude 95%: {:.2f} deg".format(
                        np.percentile(lat_diff[idx_visible], 95)
                    ),
                )
                ax[3].axvline(
                    np.percentile(lon_diff[idx_visible], 95),
                    color="b",
                    linestyle="--",
                    label="Longitude 95%: {:.2f} deg".format(
                        np.percentile(lon_diff[idx_visible], 95)
                    ),
                )
                ax[3].axvline(
                    np.percentile(lat_diff[idx_visible], 99),
                    color="r",
                    linestyle=":",
                    label="Latitude 99%: {:.2f} deg".format(
                        np.percentile(lat_diff[idx_visible], 99)
                    ),
                )
                ax[3].axvline(
                    np.percentile(lon_diff[idx_visible], 99),
                    color="b",
                    linestyle=":",
                    label="Longitude 99%: {:.2f} deg".format(
                        np.percentile(lon_diff[idx_visible], 99)
                    ),
                )
            ax[3].set_xlabel("Error [deg]")
            ax[3].set_ylabel("Density")
            ax[3].set_xlim([0, azel_xlim])
            if plot_azel_legend:
                ax[3].legend()
            ax[3].grid(True)
            ax[3].set_title("South Pole Visible (elev over 5 deg)")

            plt.tight_layout()

            if figname_azel is not None:
                plt.savefig(figname_azel, dpi=300)

            plt.show()

        # -----------------------------------------------------------------------
        if plot_azel_range_doppler:
            fig, ax = plt.subplots(2, 2, figsize=(10, 6))
            ax = ax.flatten()

            # Plot1: elevation error
            r_u2m = -user_pos_tile
            los_u2s_true = rvf_data[:, :3] - user_pos_tile
            los_u2s_true = los_u2s_true / np.linalg.norm(
                los_u2s_true, axis=1, keepdims=True
            )  # N x 3
            los_u2s_fit = fit_y[:, :3] - user_pos_tile
            los_u2s_fit = los_u2s_fit / np.linalg.norm(los_u2s_fit, axis=1, keepdims=True)  # N x 3

            lon = 0.0
            lat = -np.pi / 2
            rot_mcmf_to_enu = np.array(
                [
                    [-np.sin(lon), np.cos(lon), 0],
                    [-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)],
                    [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)],
                ]
            )  # N x 3 x 3
            los_u2s_true_enu = (rot_mcmf_to_enu @ los_u2s_true.T).T  # N x 3
            los_u2s_fit_enu = (rot_mcmf_to_enu @ los_u2s_fit.T).T  # N x 3

            elev_true = np.arcsin(los_u2s_true_enu[:, 2])
            elev_fit = np.arcsin(los_u2s_fit_enu[:, 2])
            elev_error = np.rad2deg(elev_fit - elev_true)
            az_true = np.arctan2(los_u2s_true_enu[:, 1], los_u2s_true_enu[:, 0])
            az_fit = np.arctan2(los_u2s_fit_enu[:, 1], los_u2s_fit_enu[:, 0])
            az_error = np.rad2deg(wrapToPi(az_fit - az_true))

            idx_visible = elev_true > np.deg2rad(5)  # above 5 deg elevation
            idx_est_visible = elev_fit > np.deg2rad(5)  # above 5 deg elevation
            idx_non_visible = ~idx_visible
            idx_est_non_visible = ~idx_est_visible

            num_vis_est_vis = np.sum(idx_visible & idx_est_visible) / len(t_data)
            num_vis_est_nonvis = np.sum(idx_visible & idx_est_non_visible) / len(t_data)
            num_nonvis_est_vis = np.sum(idx_non_visible & idx_est_visible) / len(t_data)
            num_nonvis_est_nonvis = np.sum(idx_non_visible & idx_est_non_visible) / len(t_data)

            fontsize = 12
            fontsize_title = 14

            # Plot 1: elevation error
            elev_error = np.rad2deg(elev_fit - elev_true)
            ax[0].plot(
                t_data[idx_visible & idx_est_visible] / 24 / 60,
                elev_error[idx_visible & idx_est_visible],
                "bo",
                markersize=ms,
                label="Visible & Estimated Visible ({:.2f}%)".format(num_vis_est_vis * 100),
            )
            ax[0].plot(
                t_data[idx_visible & idx_est_non_visible] / 24 / 60,
                elev_error[idx_visible & idx_est_non_visible],
                "ro",
                markersize=ms,
                label="Visible & Estimated Non-Visible ({:.2f}%)".format(num_vis_est_nonvis * 100),
            )
            ax[0].plot(
                t_data[idx_non_visible & idx_est_visible] / 24 / 60,
                elev_error[idx_non_visible & idx_est_visible],
                "go",
                markersize=ms,
                label="Non-Visible & Estimated Visible ({:.2f}%)".format(num_nonvis_est_vis * 100),
            )
            ax[0].plot(
                t_data[idx_non_visible & idx_est_non_visible] / 24 / 60,
                elev_error[idx_non_visible & idx_est_non_visible],
                "ko",
                markersize=ms,
                label="Non-Visible & Estimated Non-Visible ({:.2f}%)".format(
                    num_nonvis_est_nonvis * 100
                ),
            )
            ax[0].set_xlabel("Time [days]", fontsize=fontsize)
            ax[0].set_ylabel("Elevation Error [deg]", fontsize=fontsize)
            ax[0].set_title("Error in Estimated Elevation [deg]", fontsize=fontsize_title)
            ax[0].legend(loc="upper right", fontsize=10)
            ax[0].grid(True)

            # Plot 2: Azimuth error
            az_error = np.rad2deg(wrapToPi(az_fit - az_true))
            ax[1].plot(
                t_data[idx_visible] / 24 / 60,
                az_error[idx_visible],
                "bo",
                markersize=ms,
                label="Azimuth Err",
            )
            ax[1].set_xlabel("Time [days]", fontsize=fontsize)
            ax[1].set_ylabel("Azimuth Error [deg]", fontsize=fontsize)
            ax[1].set_title("Error in Estimated Azimuth (Visible) [deg]", fontsize=fontsize_title)
            ax[1].legend(loc="upper right", fontsize=10)
            ax[1].grid(True)
            ax[1].set_ylim([-10, 10])

            # Plot 3: range error
            range_pole_to_sat_true = np.linalg.norm(rvfw_data[:, :3] - user_pos_tile, axis=1)
            range_pole_to_sat_fit = np.linalg.norm(fit_y[:, :3] - user_pos_tile, axis=1)
            range_error = (range_pole_to_sat_fit - range_pole_to_sat_true) / 1000  # convert to km
            min_true_range = np.min(range_pole_to_sat_true[idx_visible]) / 1000
            max_true_range = np.max(range_pole_to_sat_true[idx_visible]) / 1000
            ax[2].plot(
                t_data[idx_visible] / 24 / 60,
                range_error[idx_visible],
                "bo",
                markersize=ms,
                label="Range Err (True Range {:.2f} to {:.2f} km)".format(
                    min_true_range, max_true_range
                ),
            )
            # ax[0].plot(t_data[~idx_visible]/24/60, range_error[~idx_visible], "ro", markersize=ms, label="Range Err (Elev below 5 deg)")
            ax[2].set_xlabel("Time [days]", fontsize=fontsize)
            ax[2].set_ylabel("Range Error [km]", fontsize=fontsize)
            ax[2].set_title("Error in Estimated Range (Visible)", fontsize=fontsize_title)
            ax[2].legend(loc="upper right", fontsize=10)
            ax[2].grid(True)

            # Plot doppler error
            fc = 2492e6  # Hz (carrier frequency)
            vel_norm_true = np.linalg.norm(rvfw_data[:, 3:], axis=1)
            vel_norm_fit = np.linalg.norm(fit_y[:, 3:], axis=1)
            unit_pos_true = rvfw_data[:, :3] / np.linalg.norm(
                rvfw_data[:, :3], axis=1, keepdims=True
            )
            unit_pos_fit = fit_y[:, :3] / np.linalg.norm(fit_y[:, :3], axis=1, keepdims=True)
            vel_radial_true = np.sum(rvfw_data[:, 3:] * unit_pos_true, axis=1)
            vel_radial_fit = np.sum(fit_y[:, 3:] * unit_pos_fit, axis=1)
            true_doppler = -fc / pnt.C * vel_radial_true
            fit_doppler = -fc / pnt.C * vel_radial_fit
            doppler_error = fit_doppler - true_doppler
            min_true_doppler = np.min(true_doppler[idx_visible])
            max_true_doppler = np.max(true_doppler[idx_visible])

            ax[3].plot(
                t_data[idx_visible] / 24 / 60,
                doppler_error[idx_visible],
                "bo",
                markersize=ms,
                label="Doppler Err (True Doppler Shift: {:.2f} to {:.2f} Hz)".format(
                    min_true_doppler, max_true_doppler
                ),
            )
            # ax[3].plot(t_data[~idx_visible]/24/60, doppler_error[~idx_visible], "ro", markersize=ms, label="Doppler Err (Elev below 5 deg)")
            ax[3].set_xlabel("Time [days]", fontsize=fontsize)
            ax[3].set_ylabel("Doppler Shift Error [Hz]", fontsize=fontsize)
            ax[3].set_title(f"Error in Estimated Doppler Shift (Visible)", fontsize=fontsize_title)
            ax[3].legend(loc="upper right", fontsize=10)
            ax[3].grid(True)

            #  Median 95% and 99% error in Elev, Azimuth (visible), Range (visible), Doppler (visible)
            elev_error = np.abs(elev_error)
            az_error = np.abs(az_error)
            range_error = np.abs(range_error)
            doppler_error = np.abs(doppler_error)
            print(
                "Elevation error: median = {:.2f} deg, 95% = {:.2f} deg, 99% = {:.2f} deg".format(
                    np.median(elev_error),
                    np.percentile(elev_error, 95),
                    np.percentile(elev_error, 99),
                )
            )
            print(
                "Azimuth error (visible): median = {:.2f} deg, 95% = {:.2f} deg, 99% = {:.2f} deg".format(
                    np.median(az_error[idx_visible]),
                    np.percentile(az_error[idx_visible], 95),
                    np.percentile(az_error[idx_visible], 99),
                )
            )
            print(
                "Range error (visible): median = {:.2f} km, 95% = {:.2f} km, 99%  = {:.2f} km".format(
                    np.median(range_error[idx_visible]),
                    np.percentile(range_error[idx_visible], 95),
                    np.percentile(range_error[idx_visible], 99),
                )
            )
            print(
                "Doppler error (visible): median = {:.2f} Hz, 95% = {:.2f} Hz, 99% = {:.2f} Hz".format(
                    np.median(doppler_error[idx_visible]),
                    np.percentile(doppler_error[idx_visible], 95),
                    np.percentile(doppler_error[idx_visible], 99),
                )
            )

            plt.tight_layout()

            if figname_azel_range_doppler is not None:
                plt.savefig(figname_azel_range_doppler, dpi=300)

            plt.show()

        # print the error statics (rms, 95%, 99%) for R, T, N, pos, Vr, Vt, Vn, ve --------------------------------

        # Compute Elevation masks
        if mask_southpole:
            idx_rmspos = idx_visible
            idx_rmsvel = idx_rmspos
            print(
                f"Number of points used for RMS calculation: {np.sum(idx_rmspos)} / {diff_y.shape[0]}"
            )
        else:
            idx_rmspos = np.arange(diff_y.shape[0])
            idx_rmsvel = np.arange(diff_y.shape[0])

        n_idxpos = idx_rmspos.size
        n_idxvel = idx_rmsvel.size

        rms_pos = np.sqrt(np.sum(diff_y[idx_rmspos, :3] ** 2, axis=0) / n_idxpos)
        median_pos = np.median(diff_y[idx_rmspos, :3], axis=0)
        std_pos = np.std(diff_y[idx_rmspos, :3], axis=0)
        p95_pos = np.percentile(np.abs(diff_y[idx_rmspos, :3]), 95, axis=0)
        p99_pos = np.percentile(np.abs(diff_y[idx_rmspos, :3]), 99, axis=0)

        rms_vel = np.sqrt(np.sum(diff_y[idx_rmsvel, 3:] ** 2, axis=0) / n_idxvel)
        median_vel = np.median(np.abs(diff_y[idx_rmsvel, 3:]), axis=0)
        std_vel = np.std(diff_y[idx_rmsvel, 3:], axis=0)
        p95_vel = np.percentile(np.abs(diff_y[idx_rmsvel, 3:]), 95, axis=0)
        p99_vel = np.percentile(np.abs(diff_y[idx_rmsvel, 3:]), 99, axis=0)

        # construct pandas dataframe
        rms_pos = np.append(rms_pos, np.linalg.norm(rms_pos))
        median_pos = np.append(median_pos, np.linalg.norm(median_pos))
        std_pos = np.append(std_pos, np.linalg.norm(std_pos))
        p95_pos = np.append(p95_pos, np.linalg.norm(p95_pos))
        p99_pos = np.append(p99_pos, np.linalg.norm(p99_pos))

        rms_vel = np.append(rms_vel, np.linalg.norm(rms_vel))
        median_vel = np.append(median_vel, np.linalg.norm(median_vel))
        std_vel = np.append(std_vel, np.linalg.norm(std_vel))
        p95_vel = np.append(p95_vel, np.linalg.norm(p95_vel))
        p99_vel = np.append(p99_vel, np.linalg.norm(p99_vel))

        stats = {
            "RMS": rms_pos,
            "Median": median_pos,
            "Std": std_pos,
            "95%": p95_pos,
            "99%": p99_pos,
        }

        stats_vel = {
            "RMS": rms_vel,
            "Median": median_vel,
            "Std": std_vel,
            "95%": p95_vel,
            "99%": p99_vel,
        }

        # set precision
        if use_rtn:
            df_pos = pd.DataFrame(stats, index=["R", "T", "N", "Total"])
            df_vel = pd.DataFrame(stats_vel, index=["Vr", "Vt", "Vn", "Total"])
        else:
            df_pos = pd.DataFrame(stats, index=["x", "y", "z", "Total"])
            df_vel = pd.DataFrame(stats_vel, index=["Vx", "Vy", "Vz", "Total"])

        # print as table
        if print_stats:
            pd.set_option("display.precision", 3)
            print("------------------------------")
            print("Position errors [{}]".format("km" if in_kms else "m"))
            print("------------------------------")
            print(df_pos)
            print("\n------------------------------")
            print("Velocity errors [{}]".format("m/s" if in_kms else "mm/s"))
            print("------------------------------")
            print(df_vel)

        return axes, df_pos, df_vel

    def eval_fit_error(
        self,
        t_data,
        rvf_data,
        rvfw_data,
        ephem,
        use_grad_for_velfit=False,
        print_stats=False,
        use_rtn=True,
        scale=None,
    ):
        """
        Evaluate the fit of the ephemeris

        Args:
            t_data: array of times
            rvf_data: array of position and velocity vectors in fixed frame
            ephem: ephemeris parameters
            axes: axes to plot the data
        """
        # compute fitted trajectory and errors --------------------------------
        lent = rvf_data.shape[0]
        if scale is None:
            fit_y = self.ephem2cart(t_data, ephem, compute_velocity=True)
        else:
            fit_y = self.ephem2cart(t_data, ephem, compute_velocity=True, scale=scale)

        # use gradient to compute velocity
        if use_grad_for_velfit:
            for i in range(3):
                fit_y[:, 3 + i] = np.gradient(fit_y[:, i], t_data, axis=0)

        # plot errors ---------------------------------------------------------
        diff_y = fit_y - rvfw_data

        # convert to rtn
        if use_rtn:
            for ti in range(lent):
                M0 = get_rtn_matrix(rvf_data[ti])
                diff_y[ti, :3] = np.dot(M0, diff_y[ti, :3])
                diff_y[ti, 3:] = np.dot(M0, diff_y[ti, 3:])

        # compute test statistics ---------------------------------------------
        diff_y_posnorm = np.linalg.norm(diff_y[:, :3], axis=1)
        diff_y_velnorm = np.linalg.norm(diff_y[:, 3:], axis=1)

        # diff_y_posnorm_p95 = np.percentile(diff_y_posnorm, 99)
        # diff_y_velnorm_p95 = np.percentile(diff_y_velnorm, 99)
        # idx_rmspos = np.where(diff_y_posnorm < diff_y_posnorm_p95)[0]
        # idx_rmsvel = np.where(diff_y_velnorm < diff_y_velnorm_p95)[0]
        # n_idxpos = idx_rmspos.size
        # n_idxvel = idx_rmsvel.size

        idx_rmspos = np.arange(lent)
        idx_rmsvel = np.arange(lent)
        n_idxpos = idx_rmspos.size
        n_idxvel = idx_rmsvel.size

        rms_pos = np.sqrt(np.sum(diff_y[idx_rmspos, :3] ** 2, axis=0) / n_idxpos)
        mean_pos = np.mean(diff_y[idx_rmspos, :3], axis=0)
        std_pos = np.std(diff_y[idx_rmspos, :3], axis=0)
        p95_pos = np.percentile(np.abs(diff_y[:, :3]), 95, axis=0)
        p99_pos = np.percentile(np.abs(diff_y[:, :3]), 99.7, axis=0)

        rms_vel = np.sqrt(np.sum(diff_y[idx_rmsvel, 3:] ** 2, axis=0) / n_idxvel)
        mean_vel = np.mean(diff_y[idx_rmsvel, 3:], axis=0)
        std_vel = np.std(diff_y[idx_rmsvel, 3:], axis=0)
        p95_vel = np.percentile(np.abs(diff_y[:, 3:]), 95, axis=0)
        p99_vel = np.percentile(np.abs(diff_y[:, 3:]), 99.7, axis=0)

        # construct pandas dataframe
        rms_pos = np.append(rms_pos, np.linalg.norm(rms_pos))
        mean_pos = np.append(mean_pos, np.linalg.norm(mean_pos))
        std_pos = np.append(std_pos, np.linalg.norm(std_pos))
        p95_pos = np.append(p95_pos, np.linalg.norm(p95_pos))
        p99_pos = np.append(p99_pos, np.linalg.norm(p99_pos))

        rms_vel = np.append(rms_vel, np.linalg.norm(rms_vel))
        mean_vel = np.append(mean_vel, np.linalg.norm(mean_vel))
        std_vel = np.append(std_vel, np.linalg.norm(std_vel))
        p95_vel = np.append(p95_vel, np.linalg.norm(p95_vel))
        p99_vel = np.append(p99_vel, np.linalg.norm(p99_vel))

        stats = {"RMS": rms_pos, "Mean": mean_pos, "Std": std_pos, "95%": p95_pos, "99.7%": p99_pos}
        stats_vel = {
            "RMS": rms_vel,
            "Mean": mean_vel,
            "Std": std_vel,
            "95%": p95_vel,
            "99.7%": p99_vel,
        }

        # set precision
        if use_rtn:
            df_pos = pd.DataFrame(stats, index=["R", "T", "N", "Total"])
            df_vel = pd.DataFrame(stats_vel, index=["Vr", "Vt", "Vn", "Total"])
        else:
            df_pos = pd.DataFrame(stats, index=["x", "y", "z", "Total"])
            df_vel = pd.DataFrame(stats_vel, index=["Vx", "Vy", "Vz", "Total"])

        # print as table
        if print_stats:
            pd.set_option("display.precision", 3)
            print(" ")
            print("------------------------------")
            print("Position errors [m]")
            print("------------------------------")
            print(df_pos)
            print("\n------------------------------")
            print("Velocity errors [mm/s]")
            print("------------------------------")
            print(df_vel)
            print(" ")

        return diff_y, df_pos, df_vel

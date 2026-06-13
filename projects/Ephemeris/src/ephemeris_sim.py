import numpy as np
import os
from tqdm import tqdm
import pylupnt as pnt
from copy import deepcopy
import matplotlib.pyplot as plt

try:
    from .keplarian_ephemeris import KeplarianEphemeris
    from .cartesian_ephemeris import CartesianEphemeris
except ImportError:
    from keplarian_ephemeris import KeplarianEphemeris
    from cartesian_ephemeris import CartesianEphemeris


class EphemerisSimulation:
    def __init__(self, data_dir):
        self.datadir = data_dir

        # store dict
        self.orbdata = {}
        self.fit_results = {}
        self.ephem_class = {}

    def setup_orbit(
        self,
        orbit_manager,
        sample_M,
        fit_mins,
        use_cheby_sampling=True,
        dt_fit=None,
        dt_eval=None,
        overwrite=False,
    ):

        self.orbm = orbit_manager
        self.sample_M = sample_M
        self.fit_mins = fit_mins
        self.use_cheby_sampling = use_cheby_sampling
        self.orbit = self.orbm.orbit

        Marray = np.linspace(0, 360, self.sample_M + 1)
        Marray = Marray[:-1]
        self.M_array = Marray

        if dt_fit is None:
            dt_fit = self.orbm.dt

        if dt_eval is None:
            dt_eval = self.orbm.dt

        point_min_fit = int(60 / dt_fit)
        point_min_eval = int(60 / dt_eval)

        print("points per min (fit):", point_min_fit)
        print("points per min (eval):", point_min_eval)

        self.dt_fit = dt_fit
        self.dt_eval = dt_eval

        dM = Marray[1] - Marray[0]
        dM_min = dM / 360 * self.orbm.period / 60
        print("Mean anomaly step size in minutes:", dM_min)

        for M in Marray:
            self.orbdata[M] = {}
            for fmin in self.fit_mins:
                self.orbdata[M][fmin] = {
                    "t_data": None,
                    "rvi": None,
                    "rvbf": None,
                    "rvbf_w": None,
                    "t_data_eval": None,
                    "rvi_eval": None,
                    "rvbf_eval": None,
                    "rvbf_w_eval": None,
                }

        M_min_orbm = np.min(self.orbm.M_prop_stack)  # propagated mean anomaly
        M_max_orbm = np.max(self.orbm.M_prop_stack)  # propagated mean anomaly

        orbdir = self.datadir + "{}/orbit/".format(self.orbit)
        if not os.path.exists(orbdir):
            os.makedirs(orbdir, exist_ok=True)

        for M in tqdm(Marray):
            for fmin in self.fit_mins:

                filename_fit = orbdir + "/M_{}_fmin_{}_dt_eval_{}_cheby_{}.npz".format(
                    M, fmin, dt_fit, self.use_cheby_sampling
                )
                filename_eval = orbdir + "/M_{}_fmin_{}__dt_eval_{}.npz".format(M, fmin, dt_fit)

                prop_M_half = (fmin / 2 * 60) / self.orbm.period * 360
                M_min_prop = M - prop_M_half
                M_max_prop = M + prop_M_half

                if (M_min_prop < M_min_orbm) and (M_max_prop > M_max_orbm):
                    print("Fitting range {} min is too large, skipping M:{}".format(fmin, M))
                    continue
                elif M_max_prop > M_max_orbm:
                    center_M = M - 360
                    if M_min_prop < M_min_orbm:
                        print("Fitting range {} min is too large, skipping M:{}".format(fmin, M))
                        continue
                elif M_min_prop < M_min_orbm:
                    center_M = M + 360
                    if M_max_prop > M_max_orbm:
                        print("Fitting range {} min is too large, skipping M:{}".format(fmin, M))
                        continue
                else:
                    center_M = M

                # fitting points ----------------------------------------------------------------------------------------------------------------
                if os.path.exists(filename_fit) and not overwrite:
                    data = np.load(filename_fit, allow_pickle=True)
                    for key, value in data.items():
                        self.orbdata[M][fmin][key] = value
                else:
                    t_data, rvi, rvbf, rvbf_w = self.orbm.extract_orbit(
                        center_M,
                        fmin,
                        method="interp",
                        n_points=fmin * point_min_fit,
                        plot_points=False,
                        use_chebyshev_sample=self.use_cheby_sampling,
                    )
                    self.orbdata[M][fmin]["t_data"] = t_data
                    self.orbdata[M][fmin]["rvi"] = rvi
                    self.orbdata[M][fmin]["rvbf"] = rvbf
                    self.orbdata[M][fmin]["rvbf_w"] = rvbf_w

                    np.savez(filename_fit, t_data=t_data, rvi=rvi, rvbf=rvbf, rvbf_w=rvbf_w)

                # evaluation points ----------------------------------------------------------------------------------------------------------------
                if os.path.exists(filename_eval) and not overwrite:
                    data = np.load(filename_eval, allow_pickle=True)
                    for key, value in data.items():
                        self.orbdata[M][fmin][key] = value
                else:
                    t_data_eval, rvi_eval, rvbf_eval, rvbf_w_eval = self.orbm.extract_orbit(
                        center_M,
                        fmin,
                        method="interp",
                        n_points=fmin * point_min_eval,
                        plot_points=False,
                        use_chebyshev_sample=False,
                    )

                    self.orbdata[M][fmin]["t_data_eval"] = t_data_eval
                    self.orbdata[M][fmin]["rvi_eval"] = rvi_eval
                    self.orbdata[M][fmin]["rvbf_eval"] = rvbf_eval
                    self.orbdata[M][fmin]["rvbf_w_eval"] = rvbf_w_eval

                    # Save the data to a file
                    np.savez(
                        filename_eval,
                        t_data_eval=t_data_eval,
                        rvi_eval=rvi_eval,
                        rvbf_eval=rvbf_eval,
                        rvbf_w_eval=rvbf_w_eval,
                    )

    def setup_kepdict(self, config):

        # Keplarian ephemeris -------------------------------------------------------
        if config["use_cheby"]:
            ptype = "cheby"
        else:
            ptype = "monomial"

        l = config["order_linear"]
        f = config["order_fourier"]
        self.kep_dict = {
            # 0 (constant), 1 (t), 2 (t^2), 3(t^3), 4(t^4), sin-cos
            "a": {"linear": 0, "fourier": 0, "fourier_sidereal": 0, "polytype": ptype},
            "e": {"linear": 0, "fourier": 0, "fourier_sidereal": 0, "polytype": ptype},
            "w": {"linear": 0, "fourier": 0, "fourier_sidereal": 0, "polytype": ptype},
            "M": {"linear": 0, "fourier": 0, "fourier_sidereal": 0, "polytype": ptype},
            "r": {"linear": l, "fourier": f, "fourier_sidereal": 0, "polytype": ptype},
            "u": {"linear": l, "fourier": f, "fourier_sidereal": 0, "polytype": ptype},
            "i": {"linear": l, "fourier": f, "fourier_sidereal": 0, "polytype": ptype},
            "l": {"linear": l, "fourier": f, "fourier_sidereal": 0, "polytype": ptype},
        }

    def setup_orders(self, config):
        orbit = self.orbit
        if orbit == "ELFO":
            orders = [2, 3, 4]
            fit_mins = [30, 60, 120, 240]
        elif orbit == "NRHO":
            orders = [2, 3, 4]
            fit_mins = [60, 120, 240, 480]
        else:
            raise ValueError("Invalid orbit. Choose 'ELFO' or 'NRHO'.")

        return orders, fit_mins

    def setup_ephem(self, ephem_type, config, print_info):
        if ephem_type == "keplarian":
            self.setup_kepdict(config)
            eph = KeplarianEphemeris(self.kep_dict, body=pnt.MOON, print_info=print_info)

            ephem_config = "keplarian_l{}_f{}_cheby_{}".format(
                config["order_linear"], config["order_fourier"], config["use_cheby"]
            )

        elif ephem_type == "cartesian":
            eph = CartesianEphemeris(
                order=config["order"],
                use_kep=config["use_kep"],
                use_rsw=config["use_rsw"],
                use_fourier=config["use_fourier"],
                use_meq=config["use_meq"],
                poly_type=config["poly_type"],
                convert_to_coe=True,
                body=pnt.MOON,
                print_info=print_info,
            )

            if config["use_meq"] and config["use_kep"]:
                use_meq = True
            else:
                use_meq = False

            ephem_config = "cartesian_{}_kep_{}_rsw_{}_fourier_{}_meq_{}_{}_sampling_{}".format(
                config["order"],
                config["use_kep"],
                config["use_rsw"],
                config["use_fourier"],
                use_meq,
                config["poly_type"],
                config["sampling_type"],
            )

        else:
            raise ValueError("Invalid ephemeris type. Choose 'cartesian' or 'polynomial'.")

        return eph, ephem_config

    def fit_ephemeris(
        self,
        ephem_type,
        fit_obj="lsq",
        config=None,
        print_errors=False,
        print_opt_results=False,
        overwrite=False,
    ):
        """
        Fit the ephemeris using the specified type.

        Parameters
        ----------
        ephem_type : str
            Type of ephemeris to fit. Options are 'chebyshev' or 'polynomial'.
        """
        eph, ephem_config = self.setup_ephem(ephem_type, config, print_info=print_errors)

        # if self.orbit == 'NRHO':
        #     use_grad_vel = True
        # else:
        #     use_grad_vel = False
        use_grad_vel = True

        savedir = self.datadir + "{}/{}".format(self.orbit, ephem_config)
        if not os.path.exists(savedir):
            # print("Creating directory: {}".format(savedir))
            os.makedirs(savedir, exist_ok=True)

        self.fit_results[ephem_config] = {}
        self.ephem_class[ephem_config] = eph

        if print_errors:
            print(" ")
            print(" ")

        for fmin in self.fit_mins:
            self.fit_results[ephem_config][fmin] = {}
            diff_pos_norm = np.zeros(0)
            diff_vel_norm = np.zeros(0)

            if print_errors:
                print(" ")
                print(" ----------------------------------------------------")
                print("  Fitting ephemeris for fmin: {}".format(fmin))
                print(" ----------------------------------------------------")

            for M in self.M_array:

                self.fit_results[ephem_config][fmin][M] = {}

                filename = savedir + "/fit_M_{}_fmin_{}.npz".format(M, fmin)

                if os.path.exists(filename) and not overwrite:
                    # if data exists, load it
                    data = np.load(filename, allow_pickle=True)
                    for key, value in data.items():
                        self.fit_results[ephem_config][fmin][M][key] = value
                    diff_xyz = self.fit_results[ephem_config][fmin][M]["diff_xyz"]
                    ephem_x = self.fit_results[ephem_config][fmin][M]["ephem"]
                else:
                    # fit ephemeris and evaluate error
                    t_data = self.orbdata[M][fmin]["t_data"]
                    rvbf = self.orbdata[M][fmin]["rvbf"]

                    t_data_eval = self.orbdata[M][fmin]["t_data_eval"]
                    rvbf_eval = self.orbdata[M][fmin]["rvbf_eval"]
                    rvbf_w_eval = self.orbdata[M][fmin]["rvbf_w_eval"]

                    # fit the ephemeris
                    ephem_x = eph.fit(t_data, rvbf, fit_obj=fit_obj, print_result=print_opt_results)

                    # evaluate error
                    diff_xyz, _, _ = eph.eval_fit_error(
                        t_data_eval,
                        rvbf_eval,
                        rvbf_w_eval,
                        ephem_x,
                        use_grad_for_velfit=use_grad_vel,
                        print_stats=False,
                        use_rtn=False,
                    )

                    # scale to meters and miillimeters
                    # diff_xyz[:, 0:3] *= 1000
                    # diff_xyz[:, 3:6] *= 1e6
                    diff_xyz[:, 3:6] *= 1e3

                    # store to dict
                    self.fit_results[ephem_config][fmin][M] = {
                        "diff_xyz": diff_xyz,
                        "ephem": ephem_x,
                    }

                    # store results
                    np.savez(filename, diff_xyz=diff_xyz, ephem=ephem_x)

                # end of if os.path.exists(filename)

                # compute error norm and display results -----------------------------------
                diff_pos_norm_tmp = np.linalg.norm(diff_xyz[:, 0:3], axis=1)
                diff_vel_norm_tmp = np.linalg.norm(diff_xyz[:, 3:6], axis=1)

                diff_pos_norm = np.append(diff_pos_norm, diff_pos_norm_tmp)
                diff_vel_norm = np.append(diff_vel_norm, diff_vel_norm_tmp)

                if print_errors:
                    print(
                        "  M: {0:.1f} | pos_norm: {2:.4f}, vel_norm: {3:.4f}".format(
                            M, fmin, np.mean(diff_pos_norm_tmp), np.mean(diff_vel_norm_tmp)
                        )
                    )

            # compute rms value
            if diff_pos_norm.size == 0:
                print("No data for fmin: {}".format(fmin))
                diff_pos_rms = None
                diff_vel_rms = None
                diff_pos_p95 = None
                diff_vel_p95 = None
            else:
                diff_pos_rms = np.sqrt(np.mean(diff_pos_norm**2))
                diff_vel_rms = np.sqrt(np.mean(diff_vel_norm**2))
                diff_pos_p95 = np.percentile(diff_pos_norm, 95)
                diff_vel_p95 = np.percentile(diff_vel_norm, 95)
                diff_pos_p997 = np.percentile(diff_pos_norm, 99.7)
                diff_vel_p997 = np.percentile(diff_vel_norm, 99.7)

            self.fit_results[ephem_config][fmin]["pos_rms"] = diff_pos_rms
            self.fit_results[ephem_config][fmin]["vel_rms"] = diff_vel_rms
            self.fit_results[ephem_config][fmin]["pos_p95"] = diff_pos_p95
            self.fit_results[ephem_config][fmin]["vel_p95"] = diff_vel_p95
            self.fit_results[ephem_config][fmin]["pos_p997"] = diff_pos_p997
            self.fit_results[ephem_config][fmin]["vel_p997"] = diff_vel_p997

            if print_errors:
                print(" ")
                print("  fit min: {}".format(fmin))
                print(
                    "  RMS pos error: {0:.4f} m, RMS vel error: {1:.4f} mm/s".format(
                        diff_pos_rms, diff_vel_rms
                    )
                )
                print(
                    "  95% pos error: {0:.4f} m, 95% vel error: {1:.4f} mm/s".format(
                        diff_pos_p95, diff_vel_p95
                    )
                )
                print(
                    "  99.7% pos error: {0:.4f} m, 99.7% vel error: {1:.4f} mm/s".format(
                        diff_pos_p997, diff_vel_p997
                    )
                )
                print(" ")

        return ephem_config

    def generate_label(self, ephem_config):
        """
        Generate a label for the ephemeris configuration.

        Parameters
        ----------
        ephem_config : str
            Ephemeris configuration string.

        Returns
        -------
        str
            Label for the ephemeris configuration.
        """
        if "cartesian" in ephem_config:
            is_kep = ephem_config.split("_")[3]
            is_rsw = ephem_config.split("_")[5]
            is_fourier = ephem_config.split("_")[7]

            if is_rsw == "True":
                label = "RSW"
            else:
                label = "XYZ"

            if is_fourier == "True":
                label = label + " (cheby: {} + fourier + kep) ".format(ephem_config.split("_")[1])
            elif is_kep == "True":
                label = label + " (cheby: {} + kep) ".format(ephem_config.split("_")[1])
            else:
                label = label + " (cheby: {}) ".format(ephem_config.split("_")[1])

        elif "keplarian" in ephem_config:
            if "True" in ephem_config:
                label = "RUIL (cheby:{} + fourier:{})".format(
                    ephem_config.split("_")[1][1:], ephem_config.split("_")[2][1:]
                )
            else:
                label = "RUIL (poly:{} + fourier:{})".format(
                    ephem_config.split("_")[1][1:], ephem_config.split("_")[2][1:]
                )
        else:
            label = ephem_config

        return label

    def plot_fit_error(
        self,
        configs=None,
        fit_mins=None,
        axes=None,
        ylim_fitmins=None,
        plot_tinv=60,
        plot_Minv=2,
        use_ma=False,
        use_bar=False,
    ):

        n_fmin = len(fit_mins)
        colors = ["b", "r", "g", "y", "m", "c", "k"]

        if axes is None:
            n_rows = 2
            n_cols = n_fmin
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_fmin, 10))

        # change text size
        plt.rcParams.update({"font.size": 14})

        if configs is None:
            configs = self.fit_results.keys()

        if fit_mins is None:
            fit_mins = self.fit_mins

        if use_bar:
            plot_Minv = 1

        n_configs = len(configs)

        for i, fmin in enumerate(fit_mins):
            axes[0][i].set_title("Fitting range: {} min".format(fmin))
            if use_ma:
                axes[0][i].set_xlabel("Mean Anomaly (deg)")
                axes[1][i].set_xlabel("Mean Anomaly (deg)")
            else:
                axes[0][i].set_xlabel("Time (min)")
                axes[1][i].set_xlabel("Time (min)")

            if use_bar:
                axes[0][i].set_ylabel("95% 3D Position Error (m)")
                axes[1][i].set_ylabel("95% 3D Velocity Error (mm/s)")
            else:
                axes[0][i].set_ylabel("3D Position Error (m)")
                axes[1][i].set_ylabel("3D Velocity Error (mm/s)")

            axes[0][i].grid(True)
            axes[1][i].grid(True)

            pos_95s = np.zeros(n_configs)
            vel_95s = np.zeros(n_configs)

            for k, ephem_config in enumerate(configs):

                # first run the fit_ephemeris function to get the data
                if ephem_config not in self.fit_results.keys():
                    print("Ephemeris config {} not found in fit results.".format(ephem_config))
                    continue

                if fmin not in self.fit_results[ephem_config].keys():
                    print("Fitting range {} not found in fit results.".format(fmin))
                    continue

                label = self.generate_label(ephem_config)

                diff_pos_p95 = self.fit_results[ephem_config][fmin]["pos_p95"]
                diff_vel_p95 = self.fit_results[ephem_config][fmin]["vel_p95"]

                pos_95s[k] = diff_pos_p95
                vel_95s[k] = diff_vel_p95

                numM = self.M_array[::plot_Minv].shape[0]
                pos_95_M = np.zeros(numM)
                vel_95_M = np.zeros(numM)
                bin_x = np.zeros(numM)

                dM = self.M_array[1] - self.M_array[0]
                dM_min = dM / 360 * self.orbm.period / 60

                # plot the data
                for j, M in enumerate(self.M_array[::plot_Minv]):
                    t_center = self.fit_results[ephem_config][fmin][M]["ephem"][0] / 60
                    t_min = t_center - fmin / 2
                    t_max = t_center + fmin / 2

                    diff_xyz = self.fit_results[ephem_config][fmin][M]["diff_xyz"]
                    n_points = diff_xyz.shape[0]
                    diff_pos_norm = np.linalg.norm(diff_xyz[:, :3], axis=1)
                    diff_vel_norm = np.linalg.norm(diff_xyz[:, 3:6], axis=1)

                    tspan = np.linspace(t_min, t_max, n_points)
                    if use_ma:
                        M_span = tspan / (self.orbm.period / 60) * 360
                        xspan = M_span
                    else:
                        xspan = tspan

                    if use_bar:
                        bin_x[j] = (xspan[0] + xspan[-1]) / 2
                        pos_95_M[j] = np.percentile(diff_pos_norm, 95)
                        vel_95_M[j] = np.percentile(diff_vel_norm, 95)
                    else:
                        if j == 0:
                            axes[0][i].plot(
                                xspan[::plot_tinv],
                                diff_pos_norm[::plot_tinv],
                                "-",
                                label=label,
                                color=colors[k],
                                alpha=0.1,
                            )
                            axes[1][i].plot(
                                xspan[1:-1:plot_tinv],
                                diff_vel_norm[1:-1:plot_tinv],
                                "-",
                                label=label,
                                color=colors[k],
                                alpha=0.1,
                            )
                        else:
                            axes[0][i].plot(
                                xspan[::plot_tinv],
                                diff_pos_norm[::plot_tinv],
                                "-",
                                color=colors[k],
                                alpha=0.1,
                            )
                            axes[1][i].plot(
                                xspan[1:-1:plot_tinv],
                                diff_vel_norm[1:-1:plot_tinv],
                                "-",
                                color=colors[k],
                                alpha=0.1,
                            )

                # plot the histogram
                if use_bar:
                    bar_width = dM if use_ma else dM_min
                    axes[0][i].bar(
                        bin_x, pos_95_M, color=colors[k], alpha=0.4, label=label, width=bar_width
                    )
                    axes[1][i].bar(
                        bin_x, vel_95_M, color=colors[k], alpha=0.4, label=label, width=bar_width
                    )

            for k, ephem_config in enumerate(configs):
                axes[0][i].axhline(y=pos_95s[k], color=colors[k], linestyle="--", alpha=1)
                axes[1][i].axhline(y=vel_95s[k], color=colors[k], linestyle="--", alpha=1)
                if use_ma:
                    xtext = 100
                else:
                    xtext = self.orbm.period / 60 * (1 / 3)

                axes[0][i].text(
                    xtext, pos_95s[k], "95%", color=colors[k], fontsize=10, ha="left", va="bottom"
                )
                axes[1][i].text(
                    xtext, vel_95s[k], "95%", color=colors[k], fontsize=10, ha="left", va="bottom"
                )

            axes[0][i].legend(loc="upper left")
            axes[1][i].legend(loc="upper left")

            axes[0][i].set_yscale("log")
            axes[1][i].set_yscale("log")

            if ylim_fitmins is not None:
                axes[0][i].set_ylim(ylim_fitmins["pos"])
                axes[1][i].set_ylim(ylim_fitmins["vel"])

        plt.tight_layout()
        plt.show()

        return axes

    def compute_resolution(self, eph, config, fmin, precision=1e-2, debug=False):

        nM = self.M_array.shape[0]
        ephem_size = self.fit_results[config][fmin][self.M_array[0]]["ephem"].shape[0]
        ephems = np.zeros((nM, ephem_size))
        k_store = np.zeros((nM, ephem_size))

        for i, M in enumerate(self.M_array):
            t_data = self.orbdata[M][fmin]["t_data"]
            n_data = t_data.shape[0]
            ephem = self.fit_results[config][fmin][M]["ephem"]
            ephems[i, :] = ephem

            # reference trajectory
            fit_y = eph.ephem2cart(t_data, ephem, compute_velocity=False)

            # identify the required fidelity
            pert_ephem = deepcopy(ephem)

            # ignore the first two parameters ["t_ref", "t_fit"]
            for j in range(2, ephem_size):
                eps = 1e-3
                pert_ephem[j] += eps
                pert_y = eph.ephem2cart(t_data, pert_ephem, compute_velocity=False)
                df_dparam = np.mean(np.linalg.norm(fit_y - pert_y, axis=1)) / eps  # derivative
                pert_ephem[j] -= eps  # reset parameter

                if df_dparam == 0:
                    df_dparam = eps
                    k = 10.0
                else:
                    k = -np.ceil(np.log2(precision / df_dparam))  # rough estimate of the resolution
                # print("  param: {}, df_dparam:{} init k: {}".format(j, df_dparam, k))
                errs = {}

                # decrease the bias until the error is small enough
                max_iter = 100
                for iter in range(max_iter):
                    eps = np.power(2, -k)
                    pert_ephem[j] += eps  # perturb the parameter to plus

                    pert_y_plus = eph.ephem2cart(
                        t_data, pert_ephem, compute_velocity=False
                    )  # t x 3
                    err_plus = np.mean(np.linalg.norm(fit_y - pert_y_plus, axis=1))

                    pert_ephem[j] -= 2 * eps  # bring to minus

                    pert_y_minus = eph.ephem2cart(t_data, pert_ephem, compute_velocity=False)

                    pert_ephem[j] += eps  # reset parameter

                    err_minus = np.mean(np.linalg.norm(fit_y - pert_y_minus, axis=1))
                    # print("    k: {}, err_plus: {}, err_minus: {}".format(k, err_plus, err_minus))

                    # check the error on both sides
                    intk = int(k)
                    errs[intk] = max(err_plus, err_minus)

                    if errs[intk] <= precision:  # within the tolerance
                        if (intk - 1) in errs.keys():
                            if errs[intk - 1] > precision:  # k-1 is not in the tolerance, but k is
                                min_k = k
                                break
                        k = k - 1  # increase the bias
                    else:  # errs[k] > precision  (not within the tolerance)
                        if (intk + 1) in errs.keys():
                            if errs[intk + 1] <= precision:  # k+1 is in the tolerance, but k is not
                                min_k = k + 1
                                break
                        k = k + 1  # decrease the bias

                    if iter == (max_iter - 1):
                        print(
                            "     Warning: max iterations reached for M: {}, fmin: {}, param: {}, k: {}".format(
                                M, fmin, j, k
                            )
                        )
                        min_k = k

                # store the k value
                k_store[i, j] = min_k

            # end of for j in range(ephem_size)
            if debug:
                print("  M: {}, k: {}".format(M, k_store[i, :]))

        # end of for i in range(nM)
        max_k = np.max(k_store, axis=0)  # maximum resolution value for each parameter

        return max_k, ephems, k_store

    def compute_bits(self, eph, max_k, ephems, config, fmin, filename, sma_keys=[], debug=False):

        ephem_size = int(self.fit_results[config][fmin][self.M_array[0]]["ephem"].shape[0])
        keys_list = eph.keys_list

        angle_keys = ["i_0", "u_0", "l_0", "w_0", "M_0", "L_0"]  # the ranges are fixed to 0-2pi
        meq_keys = ["f_0", "g_0", "h_0", "k_0"]
        ecc_keys = ["e_0"]
        bits = np.zeros(ephem_size)
        range_vals = np.zeros(ephem_size)
        max_vals = np.zeros(ephem_size)
        min_vals = np.zeros(ephem_size)

        for j in range(2, ephem_size):
            key = keys_list[j]
            if key in sma_keys:
                # use square root for the semi-major axis
                min_val = 0
                max_val = np.max(np.abs(ephems[:, j]))
                # max_k[j] = int(max_k[j]/2)
                sign_bit = 0
                margin_bit = 1
            elif key in angle_keys:  # [0, 2pi]
                min_val = 0
                max_val = 2 * np.pi
                sign_bit = 1
                margin_bit = 0
            elif key in ecc_keys:  # [0, 1]
                min_val = 0
                max_val = 1
                sign_bit = 0
                margin_bit = 0
            elif key in meq_keys:  # [-1, 1]
                min_val = 0
                max_val = 1
                sign_bit = 1
                margin_bit = 0
            else:  # coefficients
                min_val = 0
                max_val = np.max(np.abs(ephems[:, j]))
                sign_bit = 1
                margin_bit = 1

            range_val = max_val - min_val

            range_vals[j] = range_val
            max_vals[j] = max_val
            min_vals[j] = min_val

            bits[j] = (
                max(int(np.ceil(np.log2(range_val)) + max_k[j]), 1) + sign_bit + margin_bit
            )  # range of the parameter

        if debug:
            for j in range(2, ephem_size):
                print(
                    "    key: {0}  max: {1:.3e} min:{2:.3e} range:{3:.3e}  k:{4}, bits:{5}".format(
                        keys_list[j], max_vals[j], min_vals[j], range_vals[j], max_k[j], bits[j]
                    )
                )
            print("    -----------------------------------------")
            print("    total bits: {}".format(np.sum(bits)))

        self.datasizes[config][fmin]["keys"] = keys_list
        self.datasizes[config][fmin]["range"] = range_vals
        self.datasizes[config][fmin]["max"] = max_vals
        self.datasizes[config][fmin]["min"] = min_vals
        self.datasizes[config][fmin]["bits"] = bits
        self.datasizes[config][fmin]["scale"] = max_k
        self.datasizes[config][fmin]["total_bits"] = np.sum(bits)

        # save to datafile
        if filename is not None:
            np.savez(
                filename,
                keys=keys_list,
                max=max_vals,
                min=min_vals,
                range=range_vals,
                scale=max_k,
                bits=bits,
            )

    def compute_datasize(
        self, configs=None, fit_mins=None, precision=1e-2, debug=False, overwrite=False
    ):
        """
        Compute the size of the data for each ephemeris configuration and fitting range.

        Parameters
        ----------
        configs : list
            List of ephemeris configurations to compute the data size for. If None, all configurations are used.
        fit_mins : list
            List of fitting ranges to compute the data size for. If None, all fitting ranges are used.
        precision : float
            Required precision for the 3D position error. Default is 1e-5 [km] = 1 [cm]
        debug : bool
            If True, print debug information. Default is False.
        """
        if configs is None:
            configs = self.fit_results.keys()

        if fit_mins is None:
            fit_mins = self.fit_mins

        self.datasizes = {}
        for config in configs:
            self.datasizes[config] = {}

        for config in configs:
            if config not in self.fit_results.keys():
                print("Ephemeris config {} not found in fit results.".format(config))
                continue

            eph = self.ephem_class[config]

            savedir = self.datadir + "{}/{}".format(self.orbit, config)
            if not os.path.exists(savedir):
                # print("Creating directory: {}".format(savedir))
                os.makedirs(savedir, exist_ok=True)

            if debug:
                print("-----------------------------------")
                print("config: {}".format(config))
                print("-----------------------------------")

            for fmin in fit_mins:
                self.datasizes[config][fmin] = {}

                ephem_size = self.fit_results[config][fmin][self.M_array[0]]["ephem"].shape[0]

                if debug:
                    print("  [fit min: {}]".format(fmin))

                filename = savedir + "/datasize_fmin_{}.npz".format(fmin)

                if os.path.exists(filename) and not overwrite:
                    data = np.load(filename, allow_pickle=True)
                    for key, value in data.items():
                        self.datasizes[config][fmin][key] = value

                    if debug:
                        keys = data["keys"]
                        max_vals = data["max"]
                        min_vals = data["min"]
                        range_vals = data["range"]
                        max_k = data["scale"]
                        bits = data["bits"]
                        for j in range(2, ephem_size):
                            print(
                                "    key: {0}  max: {1:.3e} min:{2:.3e} range:{3:.3e}  k:{4}, bits:{5}".format(
                                    keys[j],
                                    max_vals[j],
                                    min_vals[j],
                                    range_vals[j],
                                    max_k[j],
                                    bits[j],
                                )
                            )
                        print("    -----------------------------------------")
                        print("    total bits: {}".format(np.sum(bits)))
                        print(" ")

                    bits = data["bits"]
                    self.datasizes[config][fmin]["total_bits"] = np.sum(bits)

                    continue
                else:
                    # semi-major axis keys
                    if "keplarian" in config:
                        sma_keys = ["a_0"]
                    elif "cartesian" in config:
                        if "meq" in config:
                            sma_keys = ["p_0"]
                        else:
                            sma_keys = ["A_0"]
                    else:
                        sma_keys = []

                    max_k, ephems, k_store = self.compute_resolution(
                        eph, config, fmin, precision, debug
                    )
                    self.compute_bits(
                        eph, max_k, ephems, config, fmin, filename, sma_keys=sma_keys, debug=debug
                    )

                    if debug:
                        print(" ")

    def plot_datasize_M(self, M_array, k_store, ephems, ephem_keys, config):
        """
        Plot the data size for each ephemeris configuration and fitting range.

        Parameters
        ----------
        M_array : np.ndarray
            Array of mean anomalies.
        k_store : np.ndarray  (M, ephemeris_size)
            Array of resolution values for each mean anomaly and ephemeris parameter.
        ephems : np.ndarray
            Array of ephemeris parameters.
        """
        nM = M_array.shape[0]
        n_params = k_store.shape[1]

        k_plot = k_store[:, 2:]  # skip the first two parameters (t_ref, t_fit)
        ephems_plot = ephems[:, 2:]  # skip the first two parameters (t_ref, t_fit)
        ephem_keys_plot = ephem_keys[2:]  # skip the first two parameters (t_ref, t_fit)

        # create labels for the plot ----------------------------------------------------------
        labels_sorted = []
        labels_x = []
        x_idx = []
        y_idx = []
        z_idx = []
        xf_idx = []
        yf_idx = []
        zf_idx = []
        kep_idx = []

        use_kep = config["use_kep"]
        use_rsw = config["use_rsw"]
        use_fourier = config["use_fourier"]
        use_meq = config["use_meq"]
        order = config["order"]

        # append the polynomial elements
        if use_rsw:
            poly_labels = ["a", "c", "r"]
        else:
            poly_labels = ["x", "y", "z"]

        pidx = len(kep_idx)  # start index for polynomial elements
        for pl in poly_labels:
            for j in range(order + 1):
                labels_sorted.append("{}_{}".format(pl, j))
                labels_x.append(r"${}_{{{}}}$".format(pl, j))
                if pl == "x":
                    x_idx.append(pidx)
                elif pl == "y":
                    y_idx.append(pidx)
                elif pl == "z":
                    z_idx.append(pidx)
                pidx += 1

        # then append the Fourier elements
        fidx = len(labels_sorted)  # start index for Fourier elements
        if use_fourier:
            for pl in poly_labels:
                labels_sorted.append("C_{}_{}".format(pl, order))
                labels_sorted.append("S_{}_{}".format(pl, order))
                labels_x.append(r"$c_{}$".format(pl))
                labels_x.append(r"$s_{}$".format(pl))
                if pl == "x":
                    xf_idx.append(fidx)
                    xf_idx.append(fidx + 1)
                elif pl == "y":
                    yf_idx.append(fidx)
                    yf_idx.append(fidx + 1)
                elif pl == "z":
                    zf_idx.append(fidx)
                    zf_idx.append(fidx + 1)
                fidx += 2

        index_label_sorted = []
        for label_s in labels_sorted:
            # find the index of the label in the ephemeris keys
            if label_s in ephem_keys_plot:
                index = ephem_keys_plot.index(label_s)
                index_label_sorted.append(index)
            else:
                print(f"Label {label_s} not found in ephemeris keys.")

        n_keys = len(ephem_keys_plot)

        # sort the parameters
        k_plot = k_plot[:, index_label_sorted]
        ephems_plot = ephems_plot[:, index_label_sorted]

        # plot the data size for each ephemeris parameter -----------------------
        if use_fourier:
            n_rows = 2  # first_x, mid_x, last_x, fourier_C, fourier_S
            plot_idxs = x_idx + [xf_idx[0], xf_idx[-1]]
            k_plot = k_plot[:, plot_idxs]
            ephems_plot = ephems_plot[:, plot_idxs]
            labels_x = [labels_x[i] for i in plot_idxs]
        else:
            n_rows = 1
            plot_idxs = x_idx
            k_plot = k_plot[:, plot_idxs]
            ephems_plot = ephems_plot[:, plot_idxs]
            labels_x = [labels_x[i] for i in plot_idxs]

        fig, axes = plt.subplots(n_rows, 2, figsize=(10, 4 * n_rows))

        if n_rows == 1:
            # convert axest to 1x2 array
            axes = np.array([axes])

        for i in range(n_rows):
            # In the first row plot the polynomial
            if i == 0:
                axes[i, 0].set_title("Polynomial Coefficients")

                for j, label in enumerate(labels_x):
                    if j < len(x_idx):
                        axes[i, 0].plot(M_array, k_plot[:, j], "o-", label=label, alpha=0.5)
                        axes[i, 1].plot(
                            M_array,
                            np.log2(np.abs(ephems_plot[:, j])),
                            "o-",
                            label=label,
                            alpha=0.5,
                        )

            elif n_rows == 2:
                axes[i, 0].set_title("Fourier Coefficients")

                for j, label in enumerate(labels_x):
                    if j >= len(x_idx):
                        axes[i, 0].plot(M_array, k_plot[:, j], "o-", label=label)
                        axes[i, 1].plot(
                            M_array, np.log2(np.abs(ephems_plot[:, j])), "o-", label=label
                        )

            axes[i, 0].set_ylabel("Resolution (log2)")
            axes[i, 1].set_ylabel("Absolute Value (log2)")

            axes[i, 0].set_xlabel("Mean Anomaly (deg)")
            axes[i, 1].set_xlabel("Mean Anomaly (deg)")

            axes[i, 0].grid(True)
            axes[i, 1].grid(True)

            axes[i, 0].legend(loc="upper left")
            axes[i, 1].legend(loc="upper left")

            axes[i, 0].set_ylim(10, 20)
            axes[i, 1].set_ylim(-40, 15)

        plt.tight_layout()

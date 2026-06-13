import numpy as np
import os
import matplotlib.colors as mcolors
from matplotlib import cm
import matplotlib as mpl
import matplotlib.pyplot as plt
from tqdm import tqdm
import pylupnt as pnt

try:
    from .orbit_manager import OrbitManager
    from .ephemeris_sim import EphemerisSimulation
except ImportError:
    from orbit_manager import OrbitManager
    from ephemeris_sim import EphemerisSimulation


def get_cart_ephem_config(config):
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

    return ephem_config


def get_method_label(poly_type, use_cheby_sampling, use_kep, use_fourier):

    if use_cheby_sampling:
        sampling_label = "chebyshev"
    else:
        sampling_label = "uniform"

    method_labels = [
        "coe + {} + fourier (sampling:{})".format(poly_type, sampling_label),
        "coe + {} (sampling:{})".format(poly_type, sampling_label),
        "{} (sampling:{})".format(poly_type, sampling_label),
    ]

    if use_kep and use_fourier:
        label = method_labels[0]
    elif use_kep and not use_fourier:
        label = method_labels[1]
    elif not use_kep and not use_fourier:
        label = method_labels[2]

    return label


def load_datasize_and_fit_data(
    orbit,
    use_meq,
    fit_mins,
    poly_types,
    use_kep_fouriers,
    use_cheby_samplings,
    esim_dir,
    sample_M=30,
    orders_max=20,
):

    M_array = np.linspace(0, 360, sample_M + 1)

    data = {}
    for fit_min in fit_mins:
        data[fit_min] = {}

    for fit_min in fit_mins:

        for cti, use_cheby_sampling in enumerate(use_cheby_samplings):

            for pti, poly_type in enumerate(poly_types):

                # for label in use_kep_fourier_labels:
                #     data[fit_min][label] = {}

                for j, use_kep_fourier in enumerate(use_kep_fouriers):
                    use_kep = use_kep_fourier[0]
                    use_fourier = use_kep_fourier[1]

                    label = get_method_label(poly_type, use_cheby_sampling, use_kep, use_fourier)

                    orders = np.zeros(0)
                    datasizes = np.zeros(0)
                    diff_pos_p95 = np.zeros(0)
                    diff_vel_p95 = np.zeros(0)
                    diff_pos_p997 = np.zeros(0)
                    diff_vel_p997 = np.zeros(0)

                    ephem_keys_list = []
                    range_bits_list = []
                    bits_list = []
                    scales_list = []
                    diff_pos_p95_M_list = []
                    diff_vel_p95_M_list = []

                    valid = []

                    for order in range(2, orders_max):
                        config = {
                            "order": order,
                            "use_kep": use_kep,
                            "use_rsw": False,
                            "use_fourier": use_fourier,
                            "use_meq": use_meq,
                            "poly_type": poly_type,
                            "sampling_type": "cheby" if use_cheby_sampling else "uniform",
                        }

                        ephem_config = get_cart_ephem_config(config)
                        datasize_fname = (
                            esim_dir
                            + orbit
                            + "/"
                            + ephem_config
                            + "/datasize_fmin_{}.npz".format(fit_min)
                        )

                        # load datasize data
                        if os.path.exists(datasize_fname):
                            data_ds = np.load(datasize_fname, allow_pickle=True)
                            bits = data_ds["bits"]
                            ephem_keys = data_ds["keys"]
                            range_bits = np.ceil(np.log2(data_ds["range"]))
                            scales = data_ds["scale"]
                            total_bits = np.sum(bits)
                            datasizes = np.append(datasizes, total_bits)

                            # Next load ephemeris accuracy data
                            diff_pos_norm = np.zeros(0)
                            diff_vel_norm = np.zeros(0)
                            diff_pos_p95_M = np.zeros(len(M_array))
                            diff_vel_p95_M = np.zeros(len(M_array))

                            for Mi, M in enumerate(M_array):
                                filename = (
                                    esim_dir
                                    + orbit
                                    + "/"
                                    + ephem_config
                                    + "/fit_M_{}_fmin_{}.npz".format(M, fit_min)
                                )

                                if os.path.exists(filename):
                                    # if data exists, load it
                                    data_fit = np.load(filename, allow_pickle=True)
                                    diff_xyz = data_fit["diff_xyz"]
                                    diff_pos_norm_tmp = np.linalg.norm(diff_xyz[:, 0:3], axis=1)
                                    diff_vel_norm_tmp = np.linalg.norm(diff_xyz[:, 3:6], axis=1)

                                    diff_pos_norm = np.append(diff_pos_norm, diff_pos_norm_tmp)
                                    diff_vel_norm = np.append(diff_vel_norm, diff_vel_norm_tmp)

                                    diff_pos_p95_tmp = np.percentile(diff_pos_norm_tmp, 95)
                                    diff_vel_p95_tmp = np.percentile(diff_vel_norm_tmp, 95)
                                    diff_pos_p95_M[Mi] = diff_pos_p95_tmp
                                    diff_vel_p95_M[Mi] = diff_vel_p95_tmp

                            diff_pos_p95_tmp = np.percentile(diff_pos_norm, 95)
                            diff_vel_p95_tmp = np.percentile(diff_vel_norm, 95)
                            diff_pos_p997_tmp = np.percentile(diff_pos_norm, 99.7)
                            diff_vel_p997_tmp = np.percentile(diff_vel_norm, 99.7)

                            diff_pos_p95 = np.append(diff_pos_p95, diff_pos_p95_tmp)
                            diff_vel_p95 = np.append(diff_vel_p95, diff_vel_p95_tmp)
                            diff_pos_p997 = np.append(diff_pos_p997, diff_pos_p997_tmp)
                            diff_vel_p997 = np.append(diff_vel_p997, diff_vel_p997_tmp)
                            orders = np.append(orders, order)

                            range_bits_list.append(range_bits)
                            scales_list.append(scales)
                            bits_list.append(bits)
                            ephem_keys_list.append(ephem_keys)
                            diff_pos_p95_M_list.append(diff_pos_p95_M)
                            diff_vel_p95_M_list.append(diff_vel_p95_M)

                            if (
                                total_bits <= 900
                                and diff_pos_p95_tmp <= 10
                                and diff_vel_p95_tmp <= 2.5
                            ):
                                valid.append(True)
                            else:
                                valid.append(False)

                    print("fit_min:{} {} | orders:{}".format(fit_min, label, orders))
                    data[fit_min][label] = {
                        "orbit": orbit,
                        "orders": orders,
                        "M_array": M_array,
                        "use_kep": use_kep,
                        "use_fourier": use_fourier,
                        "use_meq": use_meq,
                        "use_rsw": False,
                        "datasizes": datasizes,
                        "ephem_keys": ephem_keys_list,
                        "range_bits": range_bits_list,
                        "bits": bits_list,
                        "scales": scales_list,
                        "valid": np.array(valid),
                        "diff_pos_p95": diff_pos_p95,
                        "diff_vel_p95": diff_vel_p95,
                        "diff_pos_p99.7": diff_pos_p997,
                        "diff_vel_p99.7": diff_vel_p997,
                        "diff_pos_p95_M": diff_pos_p95_M_list,
                        "diff_vel_p95_M": diff_vel_p95_M_list,
                    }
                # end of for loop of use_kep_fourier
            # end of for loop of poly_types
        # end of for loop of use_cheby_sampling
    # end of for loop of fit_mins
    return data


def plot_fit_results_M(data, fit_min, order, config_plot, orbit, figname=None):

    poly_type = config_plot["poly_type"]
    use_cheby_sampling = config_plot["use_cheby_sampling"]
    use_kep = config_plot["use_kep"]
    use_fourier = config_plot["use_fourier"]

    # Generate the label for the method based on the config ------------------
    if use_cheby_sampling:
        sampling_label = "chebyshev"
    else:
        sampling_label = "uniform"

    method_labels = [
        "coe + {} + fourier (sampling:{})".format(poly_type, sampling_label),
        "coe + {} (sampling:{})".format(poly_type, sampling_label),
        "{} (sampling:{})".format(poly_type, sampling_label),
    ]

    if use_kep and use_fourier:
        label = method_labels[0]
    elif use_kep and not use_fourier:
        label = method_labels[1]
    elif not use_kep and not use_fourier:
        label = method_labels[2]

    # Find the order with the smallest data size that meets the accuracy requirements
    orders = data[fit_min][label]["orders"]
    order_idx = int(np.where(orders == order)[0])

    bits = data[fit_min][label]["bits"][order_idx]
    p95_pos = data[fit_min][label]["diff_pos_p95"][order_idx]
    p95_vel = data[fit_min][label]["diff_vel_p95"][order_idx]
    print(
        "Selected order: {}, data size: {} bits, 95p position error: {:.2f} m, 95p velocity error: {:.2f} mm/s".format(
            order, bits, p95_pos, p95_vel
        )
    )

    # Plot the position and velocity error at different mean anomalies for the selected order
    M_array = data[fit_min][label]["M_array"]
    diff_pos_p95_M_list = data[fit_min][label]["diff_pos_p95_M"][order_idx]
    diff_vel_p95_M_list = data[fit_min][label]["diff_vel_p95_M"][order_idx]

    # plot the arrow centered with the mean anomaly, and the length with the fit min

    plt.figure(figsize=(6, 3))
    plt.plot(M_array, diff_pos_p95_M_list, "o", label="95th Position error (m)")
    plt.plot(M_array, diff_vel_p95_M_list, "s", label="95th Velocity error (mm/s)")
    plt.xlabel("Central Mean Anomaly (degrees)")
    plt.ylabel("Fitting Error")
    plt.xlim(M_array[0], M_array[-1])
    # plt.xticks([0, 30 , 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360])
    plt.xticks([0, 36, 72, 108, 144, 180, 216, 252, 288, 324, 360])
    plt.yscale("log")
    plt.yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10])
    plt.title(
        "{} / {} mins / {} bits (order: {})".format(
            orbit,
            int(fit_min),
            int(data[fit_min][label]["datasizes"][order_idx]),
            int(orders[order_idx]),
        )
    )
    plt.grid()
    plt.legend()
    plt.tight_layout()

    if figname is not None:
        plt.savefig(figname, dpi=300)
        print("Figure saved to {}".format(figname))

    plt.show()


def plot_datasize_and_fit_data(
    data,
    pos_limit=10,
    vel_limit=2.5,
    max_bits=900,
    percentile="95",
    use_log=True,
    ylim=None,
    figname=None,
    plot_legends=True,
    labels=None,
):

    fit_mins = data.keys()
    n_fmin = len(fit_mins)

    if plot_legends:
        fig, axes = plt.subplots(2, n_fmin + 1, figsize=(4.2 * (n_fmin + 1), 4.2 * 2))
    else:
        fig, axes = plt.subplots(2, n_fmin, figsize=(4.2 * n_fmin, 4.2 * 2))

    # change text size
    plt.rcParams.update({"font.size": 14})

    max_orders = 2
    for i, fmin in enumerate(fit_mins):
        for k, ephem_config in enumerate(data[fmin].keys()):
            max_orders = max(max_orders, np.max(data[fmin][ephem_config]["orders"]))
    max_orders = int(max_orders)

    obj_poss = []
    obj_vels = []

    config_logs = {}

    for i, fmin in enumerate(fit_mins):
        # title with bold
        axes[0][i].set_title("Fitting range: {} min".format(fmin), fontweight="bold")

        configs = data[fmin].keys()
        n_configs = len(configs)
        cmaps = ["Reds", "Greens", "Blues", "Purples"]

        print(" ")
        print("Fitting range: {} min".format(fmin))

        for k, ephem_config in enumerate(configs):

            config_logs[ephem_config] = {
                "success_orders": [],
                "success_p95": [],
                "success_v95": [],
                "success_bits": [],
            }

            min_pos_err = np.inf
            min_vel_err = None
            min_posvel_order = None
            min_posvel_bits = None

            orders = data[fmin][ephem_config]["orders"]
            total_bits = data[fmin][ephem_config]["datasizes"]
            pos_p95 = data[fmin][ephem_config]["diff_pos_p{}".format(percentile)]
            vel_p95 = data[fmin][ephem_config]["diff_vel_p{}".format(percentile)]

            color = plt.get_cmap(cmaps[k])(np.linspace(0.4, 1, max_orders + 1))

            print("  {}".format(ephem_config))

            # row0: poisition vs data size
            for j in range(orders.size):
                if total_bits[j] <= 900 and pos_p95[j] <= pos_limit and vel_p95[j] <= vel_limit:
                    success = True
                    markeredgecolor = "black"
                    markeredgewidth = 1.5
                    config_logs[ephem_config]["success_orders"].append(orders[j])
                    config_logs[ephem_config]["success_p95"].append(pos_p95[j])
                    config_logs[ephem_config]["success_v95"].append(vel_p95[j])
                    config_logs[ephem_config]["success_bits"].append(total_bits[j])
                else:
                    success = False
                    markeredgecolor = "none"
                    markeredgewidth = 0

                if total_bits[j] <= 900:
                    if pos_p95[j] < min_pos_err:
                        min_pos_err = pos_p95[j]
                        min_vel_err = vel_p95[j]
                        min_posvel_order = orders[j]
                        min_posvel_bits = total_bits[j]

                markersize = 10

                if i == (n_fmin - 1) and j == 0:
                    (obj_pos,) = axes[0][i].plot(
                        total_bits[j],
                        pos_p95[j],
                        "o",
                        markeredgecolor=markeredgecolor,
                        markeredgewidth=markeredgewidth,
                        markersize=markersize,
                        label=labels[k] if labels is not None else ephem_config,
                        color=color[int(orders[j])],
                    )
                    # row1: velocity vs data size
                    (obj_vel,) = axes[1][i].plot(
                        total_bits[j],
                        vel_p95[j],
                        "o",
                        markeredgecolor=markeredgecolor,
                        markeredgewidth=markeredgewidth,
                        markersize=markersize,
                        label=labels[k] if labels is not None else ephem_config,
                        color=color[int(orders[j])],
                    )

                    obj_poss.append(obj_pos)
                    obj_vels.append(obj_vel)

                else:
                    axes[0][i].plot(
                        total_bits[j],
                        pos_p95[j],
                        "o",
                        markeredgecolor=markeredgecolor,
                        markeredgewidth=markeredgewidth,
                        markersize=markersize,
                        color=color[j],
                    )
                    # row1: velocity vs data size
                    axes[1][i].plot(
                        total_bits[j],
                        vel_p95[j],
                        "o",
                        markeredgecolor=markeredgecolor,
                        markeredgewidth=markeredgewidth,
                        markersize=markersize,
                        color=color[j],
                    )

            # end of for loop of orders

            # print the config logs
            suc_orders = config_logs[ephem_config]["success_orders"]
            suc_orders_str = ", ".join([str(int(o)) for o in suc_orders])
            print(f"    Sucess Orders: {suc_orders_str}")
            # position percentile with .2f
            suc_postion_str = ", ".join(
                [f"{p:.4f}" for p in config_logs[ephem_config]["success_p95"]]
            )
            print(f"    Sucess Position {percentile}p: {suc_postion_str}")
            # velocity percentile with .2f
            suc_velocity_str = ", ".join(
                [f"{v:.4f}" for v in config_logs[ephem_config]["success_v95"]]
            )
            print(f"    Sucess Velocity {percentile}p: {suc_velocity_str}")
            print(
                f"    Success Bits: {[int(b) for b in config_logs[ephem_config]['success_bits']]}"
            )

            print(" ")
            print("   Minimum position error configuration ========================")
            print(
                f"     Order: {min_posvel_order}, Bits: {int(min_posvel_bits)}"
                + f"     Pos {percentile}p: {min_pos_err:.2e} m"
                f"     Vel {percentile}p: {min_vel_err:.2e} mm/s"
            )
            print("   ==============================================================")
            print(" ")

        if i == (n_fmin - 1) and plot_legends:
            # make axes [0][i+1] and [1][i+1] invisible
            # do not show the boxes, except the legend
            for rows in range(2):
                axes[rows][n_fmin].spines["top"].set_visible(False)
                axes[rows][n_fmin].spines["bottom"].set_visible(False)
                axes[rows][n_fmin].spines["left"].set_visible(False)
                axes[rows][n_fmin].spines["right"].set_visible(False)
                # no xticks and yticks
                axes[rows][n_fmin].set_xticks([])
                axes[rows][n_fmin].set_yticks([])

            axes[0][n_fmin].legend(
                handles=obj_poss, loc="upper left", fontsize=16, bbox_to_anchor=(0, 1)
            )
            axes[1][n_fmin].legend(
                handles=obj_vels, loc="upper left", fontsize=16, bbox_to_anchor=(0, 1)
            )

            # colorbar
            colorbar1 = {}
            colorbar2 = {}
            for ci in range(n_configs):
                truncated_cmap = mpl.colors.ListedColormap(
                    plt.get_cmap(cmaps[ci])(np.linspace(0.4, 1, max_orders + 1))
                )

                colorbar1[ci] = fig.colorbar(
                    cm.ScalarMappable(
                        norm=mpl.colors.Normalize(vmin=0, vmax=max_orders + 1), cmap=truncated_cmap
                    ),
                    ax=axes[0][4],
                    orientation="horizontal",
                    location="bottom",
                    pad=0.12 * (ci + 1),
                )
                colorbar2[ci] = fig.colorbar(
                    cm.ScalarMappable(
                        norm=mpl.colors.Normalize(vmin=0, vmax=max_orders + 1), cmap=truncated_cmap
                    ),
                    ax=axes[1][4],
                    orientation="horizontal",
                    location="bottom",
                    pad=0.12 * (ci + 1),
                )
                if ci == 0:
                    colorbar1[ci].set_label("Order")
                    colorbar2[ci].set_label("Order")

        axes[0][i].set_xlabel("Ephemeris size (bits)", fontsize=14)
        axes[0][i].set_ylabel("95% Position Error (m)", fontsize=14)
        axes[0][i].set_xlim(0, 1400)
        axes[0][i].set_xticks(np.arange(0, 1500, 200))
        axes[0][i].axhline(y=pos_limit, color="k", linestyle="--")
        axes[0][i].axvline(x=max_bits, color="k", linestyle="--")
        if use_log:
            axes[0][i].set_yscale("log")
            if ylim is not None:
                axes[0][i].set_ylim(ylim)
            else:
                axes[0][i].set_ylim(1e-4, 1e4)
            axes[0][i].set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 10000])
        else:
            axes[0][i].set_ylim(ylim)
        axes[0][i].grid(True)
        if i == (n_fmin - 1):
            axes[0][i].legend()

        axes[1][i].set_xlabel("Data size (bits)")
        axes[1][i].set_ylabel("95p Velocity Error (mm/s)")
        axes[1][i].set_xlim(0, 1400)
        axes[1][i].set_xticks(np.arange(0, 1500, 200))
        if use_log:
            axes[1][i].set_yscale("log")
            if ylim is not None:
                axes[1][i].set_ylim(ylim)
            else:
                axes[1][i].set_ylim(1e-3, 1e6)
            axes[0][i].set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 10000])
        else:
            axes[1][i].set_ylim(ylim)
        axes[1][i].axhline(y=vel_limit, color="k", linestyle="--")
        axes[1][i].axvline(x=max_bits, color="k", linestyle="--")
        axes[1][i].grid(True)
        if i == (n_fmin - 1):
            axes[1][i].legend()
        # colorbar
        # colorbar = plt.colorbar(color, ax=axes[0][n_fmin-1], orientation='horizontal')
        # colorbar.set_label("Order")

    plt.tight_layout()
    if figname is not None:
        plt.savefig(figname, dpi=300)
        print("Figure saved to {}".format(figname))
    plt.show()


def plot_bits(
    data,
    fit_min,
    bit_only=False,
    use_max_order=False,
    use_order=None,
    show_legend=True,
    figdir=None,
    add_title=True,
):
    """
    Plot the data bits for each ephemeris message type.
    """

    data_fitmin = data[fit_min]
    n_data = len(data_fitmin)

    n_fig = len(data_fitmin.items())

    for i, (label, values) in enumerate(data_fitmin.items()):
        valid = values["valid"]
        bits = values["bits"]
        total_bits = values["datasizes"]
        range_bits = values["range_bits"]
        ephem_keys = values["ephem_keys"]
        scales = values["scales"]
        orders = values["orders"]
        use_kep = values["use_kep"]
        use_fourier = values["use_fourier"]
        use_meq = values["use_meq"]
        use_rsw = values["use_rsw"]
        p95_pos = values["diff_pos_p95"]
        p95_vel = values["diff_vel_p95"]

        # find the least order that satisfies the constraints
        if use_max_order:
            max_valid_index = len(bits) - 1
        elif use_order is not None:
            max_valid_index = np.where(orders <= use_order)[0][-1]
        else:
            valid_indices = np.where(valid)[0]
            if valid_indices.size > 0:
                max_valid_index = valid_indices[-1]
            else:
                max_valid_index = np.where(total_bits <= 900)[0][-1]

        max_order = orders[max_valid_index]
        max_bits = bits[max_valid_index][2:]
        max_range_bits = range_bits[max_valid_index][2:]
        max_scales = scales[max_valid_index][2:]
        max_ephem_keys = ephem_keys[max_valid_index][2:]
        max_total_bits = total_bits[max_valid_index]
        max_p95_pos = p95_pos[max_valid_index]
        max_p95_vel = p95_vel[max_valid_index]

        max_order = int(max_order)
        max_ephem_keys = max_ephem_keys.tolist()

        # create labels for the plot
        labels_sorted = []
        labels_x = []
        x_idx = []
        y_idx = []
        z_idx = []
        xf_idx = []
        yf_idx = []
        zf_idx = []
        kep_idx = []

        # first append keplarian elements if they are used
        if use_kep:
            if use_meq:
                kep_labels = ["p_0", "f_0", "g_0", "h_0", "k_0", "L_0"]
                labels_x_kep = ["p", "f", "g", "h", "k", "L"]
            else:
                kep_labels = ["A_0", "e_0", "i_0", "l_0", "w_0", "M_0"]
                labels_x_kep = ["a", "e", "i", "l", "w", "M"]
            labels_sorted.extend(kep_labels)
            labels_x.extend(labels_x_kep)
            kep_idx = [i for i in range(6)]
        else:
            kep_idx = []

        # then append the polynomial elements
        if use_rsw:
            poly_labels = ["a", "c", "r"]
        else:
            poly_labels = ["x", "y", "z"]

        pidx = len(kep_idx)  # start index for polynomial elements
        for pl in poly_labels:
            for j in range(max_order + 1):
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
                labels_sorted.append("C_{}_{}".format(pl, max_order))
                labels_sorted.append("S_{}_{}".format(pl, max_order))
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

        # print(f"Labels sorted ({len(labels_sorted)}): {labels_sorted}")
        # print(f"Labels x ({len(labels_x)}): {labels_x}")
        # print(f"max_ephem_keys ({len(max_ephem_keys)}): {max_ephem_keys}")
        # print(f"bits ({len(max_bits)}): {max_bits}")
        # print(f"range_bits ({len(max_range_bits)}): {max_range_bits}")
        # print(f"scales ({len(max_scales)}): {max_scales}")

        index_label_sorted = []
        for label_s in labels_sorted:
            # find the index of the label in the ephemeris keys
            if label_s in max_ephem_keys:
                index = max_ephem_keys.index(label_s)
                index_label_sorted.append(index)
            else:
                print(f"Label {label_s} not found in ephemeris keys.")

        n_keys = len(max_ephem_keys)

        bits_sorted = max_bits[index_label_sorted]
        range_bits_sorted = max_range_bits[index_label_sorted]
        scales_sorted = max_scales[index_label_sorted]

        # create a new figure for each ephemeris type
        if bit_only:
            fig, ax = plt.subplots(figsize=(8, 4))
        else:
            fig, axes = plt.subplots(1, 3, figsize=(24, 5))

        plt.rcParams["text.usetex"] = True

        if add_title:
            if not bit_only:
                fig.suptitle(
                    f"Ephemeris Type: {label}  \\ Fit Min: {fit_min} \\ Order: {max_order} \n Total Bits: {int(max_total_bits)} \\ 95\% Position SISE: {max_p95_pos:.2f} m \\ 95\% Velocity SISE: {max_p95_vel:.2f} mm/s",
                    fontsize=16,
                    fontweight="bold",
                )
            else:
                fig.suptitle(
                    f"Ephemeris Type: {label}  \n Fit Min: {fit_min}  \\ Order: {max_order} \\ Total Bits: {int(max_total_bits)} \n 95\% Position SISE: {max_p95_pos:.2f} m \\ 95\% Velocity SISE: {max_p95_vel:.2f} mm/s",
                    fontsize=16,
                    fontweight="bold",
                )

        colors = ["purple", "red", "green", "blue", "darkred", "darkgreen", "darkblue"]

        if max_order >= 15:
            fs = 12
        elif max_order >= 10:
            fs = 14
        else:
            fs = 16

        # in first subplot, plot the bits
        if not bit_only:
            ax = axes[0]
        ax.bar(kep_idx, bits_sorted[kep_idx], color=colors[0], alpha=0.6, label="Keplarian")
        ax.bar(x_idx, bits_sorted[x_idx], color=colors[1], alpha=0.6, label="X Polynomial")
        ax.bar(y_idx, bits_sorted[y_idx], color=colors[2], alpha=0.6, label="Y Polynomial")
        ax.bar(z_idx, bits_sorted[z_idx], color=colors[3], alpha=0.6, label="Z Polynomial")
        if use_fourier:
            ax.bar(xf_idx, bits_sorted[xf_idx], color=colors[4], alpha=0.6, label="X Fourier")
            ax.bar(yf_idx, bits_sorted[yf_idx], color=colors[5], alpha=0.6, label="Y Fourier")
            ax.bar(zf_idx, bits_sorted[zf_idx], color=colors[6], alpha=0.6, label="Z Fourier")
        ax.set_xticks(range(n_keys))
        ax.set_xticklabels(labels_x, rotation=90, fontsize=fs)
        ax.set_ylabel("Bits")
        if not bit_only:
            ax.set_title("Bits")
        ax.grid(True)

        if show_legend:
            ax.legend(loc="upper right", fontsize=8)

        if not bit_only:
            # in second subplot, plot the range bits
            ax = axes[1]
            ax.bar(
                kep_idx, range_bits_sorted[kep_idx], color=colors[0], alpha=0.6, label="Range Bits"
            )
            ax.bar(
                x_idx,
                range_bits_sorted[x_idx],
                color=colors[1],
                alpha=0.6,
                label="X Polynomial Range Bits",
            )
            ax.bar(
                y_idx,
                range_bits_sorted[y_idx],
                color=colors[2],
                alpha=0.6,
                label="Y Polynomial Range Bits",
            )
            ax.bar(
                z_idx,
                range_bits_sorted[z_idx],
                color=colors[3],
                alpha=0.6,
                label="Z Polynomial Range Bits",
            )
            if use_fourier:
                ax.bar(
                    xf_idx,
                    range_bits_sorted[xf_idx],
                    color=colors[4],
                    alpha=0.6,
                    label="X Fourier Range Bits",
                )
                ax.bar(
                    yf_idx,
                    range_bits_sorted[yf_idx],
                    color=colors[5],
                    alpha=0.6,
                    label="Y Fourier Range Bits",
                )
                ax.bar(
                    zf_idx,
                    range_bits_sorted[zf_idx],
                    color=colors[6],
                    alpha=0.6,
                    label="Z Fourier Range Bits",
                )
            ax.set_xticks(range(n_keys))
            ax.set_xticklabels(labels_x, rotation=90, fontsize=fs)
            ax.set_ylabel("Range (log2)")
            ax.set_title("Range (log2)")
            ax.grid(True)

            # in third subplot, plot the scales
            ax = axes[2]
            ax.bar(kep_idx, scales_sorted[kep_idx], color=colors[0], alpha=0.6, label="Scales")
            ax.bar(
                x_idx, scales_sorted[x_idx], color=colors[1], alpha=0.6, label="X Polynomial Scales"
            )
            ax.bar(
                y_idx, scales_sorted[y_idx], color=colors[2], alpha=0.6, label="Y Polynomial Scales"
            )
            ax.bar(
                z_idx, scales_sorted[z_idx], color=colors[3], alpha=0.6, label="Z Polynomial Scales"
            )
            if use_fourier:
                ax.bar(
                    xf_idx,
                    scales_sorted[xf_idx],
                    color=colors[4],
                    alpha=0.6,
                    label="X Fourier Scales",
                )
                ax.bar(
                    yf_idx,
                    scales_sorted[yf_idx],
                    color=colors[5],
                    alpha=0.6,
                    label="Y Fourier Scales",
                )
                ax.bar(
                    zf_idx,
                    scales_sorted[zf_idx],
                    color=colors[6],
                    alpha=0.6,
                    label="Z Fourier Scales",
                )
            ax.set_xticks(range(n_keys))
            ax.set_xticklabels(labels_x, rotation=90, fontsize=fs)
            ax.set_ylabel("Scales")
            ax.set_title("Scales (log2)")
            ax.grid(True)

        plt.tight_layout()
        if figdir is not None:
            fig.savefig(
                figdir
                + "/bits_fitmin_{}_type_{}.pdf".format(
                    fit_min,
                    label.replace(" ", "_")
                    .replace("+", "")
                    .replace("(", "")
                    .replace(")", "")
                    .replace(":", "_")
                    .replace("__", "_"),
                ),
                dpi=300,
            )


def create_esim(orbit, config, fit_min, basedir):

    # basedir = "/Users/keidaiiiyama/Documents/sw_navlab/LuPNT/projects/LunarEphem/"
    orbm_save_dir = basedir + "data/orbits"
    esim_dir = basedir + "data/ephemeris/"

    if orbit == "ELFO" or orbit == "Polar" or orbit == "CLFO":
        dt = 0.1
        n_period = 3
    else:
        dt = 10.0
        n_period = 2

    dt_fit = 60.0
    dt_eval = 1.0
    use_cheby_sampling = config["sampling_type"] == "cheby"

    orbm = OrbitManager(
        orbit, dyn=None, n_period=n_period, dt=dt, data_dir=orbm_save_dir, overwrite=False
    )
    esim = EphemerisSimulation(data_dir=esim_dir)
    esim.setup_orbit(
        orbm,
        sample_M=30,
        fit_mins=[fit_min],
        use_cheby_sampling=use_cheby_sampling,
        dt_fit=dt_fit,
        dt_eval=dt_eval,
        overwrite=False,
    )

    return esim


def plot_ephemsize_M(orbit, config, fit_min, basedir):

    esim = create_esim(orbit, config, fit_min, basedir)
    ephem_config = esim.fit_ephemeris(
        ephem_type="cartesian",
        config=config,
        print_errors=False,
        print_opt_results=False,
        overwrite=False,
    )
    eph = esim.ephem_class[ephem_config]
    max_k, ephems, k_store = esim.compute_resolution(
        eph, ephem_config, fit_min, precision=1e-5, debug=False
    )
    esim.plot_datasize_M(esim.M_array, k_store, ephems, eph.keys_list, config)


def test_fitting_L2(orbit, config, fit_min, alphas, basedir, fit_obj="lsq"):
    esim = create_esim(orbit, config, fit_min, basedir)
    eph, ephem_config = esim.setup_ephem(ephem_type="cartesian", config=config, print_info=False)

    fit_results = {}
    for alpha in alphas:
        fit_results[alpha] = {}
        fit_results[alpha][fit_min] = {}

        for i, M in enumerate(esim.M_array):
            fit_results[alpha][fit_min][M] = {"diff_pos_norm": [], "diff_vel_norm": [], "ephem": []}

    # Fit ephemeris -------------------------------------------
    print("Fitting ephemeris for each M...")
    nM = esim.M_array.size
    for i in tqdm(range(nM), desc="Fitting ephemeris for each M"):
        M = esim.M_array[i]
        # fit ephemeris and evaluate error
        t_data = esim.orbdata[M][fit_min]["t_data"]
        rvbf = esim.orbdata[M][fit_min]["rvbf"]
        t_data_eval = esim.orbdata[M][fit_min]["t_data_eval"]
        rvbf_eval = esim.orbdata[M][fit_min]["rvbf_eval"]
        rvbf_w_eval = esim.orbdata[M][fit_min]["rvbf_w_eval"]

        # fit the ephemeris
        ephem_x_vec = eph.fit(t_data, rvbf, fit_obj=fit_obj, print_result=False, L2_alpha=alphas)

        # evaluate error
        for j, ephem_x in enumerate(ephem_x_vec):
            if np.isnan(ephem_x).any():
                fit_results[alphas[j]][fit_min][M]["diff_pos_norm"] = np.array([np.nan])
                fit_results[alphas[j]][fit_min][M]["diff_vel_norm"] = np.array([np.nan])
                fit_results[alphas[j]][fit_min][M]["ephem"] = np.array([np.nan])
                continue

            diff_xyz, df_pos, df_vel = eph.eval_fit_error(
                t_data_eval,
                rvbf_eval,
                rvbf_w_eval,
                ephem_x,
                use_grad_for_velfit=True,
                print_stats=False,
                use_rtn=False,
            )

            # scale to meters and miillimeters
            # diff_xyz[:, 0:3] *= 1000
            # diff_xyz[:, 3:6] *= 1e6
            diff_xyz[:, 3:6] *= 1e3

            diff_pos_norm_tmp = np.linalg.norm(diff_xyz[:, 0:3], axis=1)
            diff_vel_norm_tmp = np.linalg.norm(diff_xyz[:, 3:6], axis=1)

            fit_results[alphas[j]][fit_min][M]["diff_pos_norm"] = diff_pos_norm_tmp
            fit_results[alphas[j]][fit_min][M]["diff_vel_norm"] = diff_vel_norm_tmp
            fit_results[alphas[j]][fit_min][M]["ephem"] = ephem_x

            # print(
            #     "alpha: {0:.3f}  M: {1:.1f}  | pos: {2:.3f} m  vel: {3:.3f} mm/s".format(
            #         alphas[j],
            #         M,
            #         np.percentile(diff_pos_norm_tmp, 95),
            #         np.percentile(diff_vel_norm_tmp, 95),
            #     )
            # )

    # compute 95% error for each alpha ------------------------------
    print("Computing 95 percentile error for each alpha...")
    M_all_valid = {}

    for alpha in alphas:
        diff_pos_norm = np.zeros(0)
        diff_vel_norm = np.zeros(0)
        M_all_valid[alpha] = True

        for M in esim.M_array:
            if np.isnan(fit_results[alpha][fit_min][M]["diff_pos_norm"]).all():
                M_all_valid[alpha] = False
                print(f"alpha: {alpha:.3e}  M: {M:.1f} has NaN values. Skipping this alpha.")
                break

            diff_pos_norm = np.append(
                diff_pos_norm, fit_results[alpha][fit_min][M]["diff_pos_norm"]
            )
            diff_vel_norm = np.append(
                diff_vel_norm, fit_results[alpha][fit_min][M]["diff_vel_norm"]
            )

        if not M_all_valid[alpha]:
            fit_results[alpha][fit_min]["pos_p95"] = np.nan
            fit_results[alpha][fit_min]["vel_p95"] = np.nan
            fit_results[alpha][fit_min]["pos_rms"] = np.nan
            fit_results[alpha][fit_min]["vel_rms"] = np.nan
        else:
            fit_results[alpha][fit_min]["pos_p95"] = np.percentile(diff_pos_norm, 95)
            fit_results[alpha][fit_min]["vel_p95"] = np.percentile(diff_vel_norm, 95)
            fit_results[alpha][fit_min]["pos_rms"] = np.sqrt(np.mean(diff_pos_norm**2))
            fit_results[alpha][fit_min]["vel_rms"] = np.sqrt(np.mean(diff_vel_norm**2))

    esim.fit_results = fit_results
    esim.datasizes = {}

    # compute datasize for each alpha ------------------------------
    print("Computing datasize for each alpha...")
    for alpha in alphas:
        if not M_all_valid[alpha]:
            esim.datasizes[alpha] = {}
            esim.datasizes[alpha][fit_min] = {
                "keys": np.array([np.nan]),
                "max": np.nan,
                "min": np.nan,
                "scale": np.array([np.nan]),
                "range ": np.array([np.nan]),
                "bits": np.array([np.nan]),
                "total_bits": np.nan,
            }
        else:
            esim.datasizes[alpha] = {}
            esim.datasizes[alpha][fit_min] = {}
            max_k, ephems, k_store = esim.compute_resolution(
                eph, alpha, fit_min, precision=1e-2, debug=False
            )
            esim.compute_bits(
                eph, max_k, ephems, alpha, fit_min, filename=None, sma_keys=["A_0"], debug=False
            )

    return esim, eph


def plot_fitting_l2_results(fig, axes, esim_data, fit_min, alphas, plotlabel, plot_thresholds=True):
    """
    Plot the results of the fitting L2 method.
    """
    fit_results = esim_data["fit_results"]
    datasizes = esim_data["datasizes"]
    M_array = esim_data["M_array"]
    orbdata = esim_data["orbdata"]
    config = esim_data["config"]
    orbit = esim_data["orbit"]
    basedir = esim_data["basedir"]

    esim = create_esim(orbit, config, fit_min, basedir)
    eph, ephem_config = esim.setup_ephem(ephem_type="cartesian", config=config, print_info=False)

    nM = M_array.size

    # fig, axes = plt.subplots(2, 3, figsize=(18, 8))

    # recompute fit error using descretized ephemeris
    is_valid = {}

    for alpha in alphas:
        if np.isnan(datasizes[alpha][fit_min]["scale"]).any():
            is_valid[alpha] = False
        else:
            is_valid[alpha] = True

    recompute = False
    if recompute:
        for i in range(nM):
            M = M_array[i]
            t_data_eval = orbdata[M][fit_min]["t_data_eval"]
            rvbf_eval = orbdata[M][fit_min]["rvbf_eval"]
            rvbf_w_eval = orbdata[M][fit_min]["rvbf_w_eval"]

            for alpha in alphas:
                if is_valid[alpha] == True:
                    ephem_x = fit_results[alpha][fit_min][M]["ephem"]
                    scale = datasizes[alpha][fit_min]["scale"]

                    diff_xyz, df_pos, df_vel = eph.eval_fit_error(
                        t_data_eval,
                        rvbf_eval,
                        rvbf_w_eval,
                        ephem_x,
                        use_grad_for_velfit=True,
                        print_stats=False,
                        use_rtn=False,
                        scale=scale,
                    )
                    diff_xyz[:, 3:6] *= 1e3

                    diff_pos_norm_tmp = np.linalg.norm(diff_xyz[:, 0:3], axis=1)
                    diff_vel_norm_tmp = np.linalg.norm(diff_xyz[:, 3:6], axis=1)

                    fit_results[alpha][fit_min][M]["diff_pos_norm_s"] = diff_pos_norm_tmp
                    fit_results[alpha][fit_min][M]["diff_vel_norm_s"] = diff_vel_norm_tmp
                    fit_results[alpha][fit_min][M]["ephem_s"] = ephem_x
    else:
        for i in range(nM):
            M = M_array[i]
            for alpha in alphas:
                if is_valid[alpha] == True:
                    fit_results[alpha][fit_min][M]["diff_pos_norm_s"] = fit_results[alpha][fit_min][
                        M
                    ]["diff_pos_norm"]
                    fit_results[alpha][fit_min][M]["diff_vel_norm_s"] = fit_results[alpha][fit_min][
                        M
                    ]["diff_vel_norm"]
                    fit_results[alpha][fit_min][M]["ephem_s"] = fit_results[alpha][fit_min][M][
                        "ephem"
                    ]

    # compute 95% error for each alpha ------------------------------
    # print("Computing 95 percentile error for each alpha...")
    for alpha in alphas:
        if is_valid[alpha] == True:
            diff_pos_norm = np.zeros(0)
            diff_vel_norm = np.zeros(0)
            for M in M_array:
                diff_pos_norm = np.append(
                    diff_pos_norm, fit_results[alpha][fit_min][M]["diff_pos_norm_s"]
                )
                diff_vel_norm = np.append(
                    diff_vel_norm, fit_results[alpha][fit_min][M]["diff_vel_norm_s"]
                )
            fit_results[alpha][fit_min]["pos_p95s"] = np.percentile(diff_pos_norm, 95)
            fit_results[alpha][fit_min]["vel_p95s"] = np.percentile(diff_vel_norm, 95)
            fit_results[alpha][fit_min]["pos_rmss"] = np.sqrt(np.mean(diff_pos_norm**2))
            fit_results[alpha][fit_min]["vel_rmss"] = np.sqrt(np.mean(diff_vel_norm**2))
        else:
            fit_results[alpha][fit_min]["pos_p95s"] = np.nan
            fit_results[alpha][fit_min]["vel_p95s"] = np.nan
            fit_results[alpha][fit_min]["pos_rmss"] = np.nan
            fit_results[alpha][fit_min]["vel_rmss"] = np.nan

    # for the first column, plot the alpha vs pos_p95 (1st row) and vel_p95 (2nd row)
    axes[0].set_title(r"Position Error (95 perc) vs Optimization weight $\gamma$")
    axes[0].set_xlabel(r"Optimization weight $\gamma$")
    axes[0].set_ylabel("Position Error (m)")
    axes[0].set_xscale("log")
    axes[0].grid(True)
    pos_p95_alpha = np.zeros(alphas.size)
    vel_p95_alpha = np.zeros(alphas.size)
    pos_p95_alpha_s = np.zeros(alphas.size)
    vel_p95_alpha_s = np.zeros(alphas.size)
    for i, alpha in enumerate(alphas):
        pos_p95_alpha[i] = fit_results[alpha][fit_min]["pos_p95"]
        vel_p95_alpha[i] = fit_results[alpha][fit_min]["vel_p95"]
        pos_p95_alpha_s[i] = fit_results[alpha][fit_min]["pos_p95s"]
        vel_p95_alpha_s[i] = fit_results[alpha][fit_min]["vel_p95s"]
    # axes[0][0].plot(alphas, pos_p95_alpha, 'ro--', label='Fitted')
    axes[0].plot(alphas, pos_p95_alpha_s, "o-", label=plotlabel)

    # for the second column, plot the alpha vs total bits (1st row) and range bits (2nd row)
    axes[1].set_title(r"Total Bits vs Optimization weight $\gamma$")
    axes[1].set_xlabel(r"Optimization weight $\gamma$")
    axes[1].set_ylabel("Total Bits")
    axes[1].grid(True)
    axes[1].set_xscale("log")
    total_bits_alpha = np.zeros(alphas.size)
    for i, alpha in enumerate(alphas):
        total_bits_alpha[i] = datasizes[alpha][fit_min]["total_bits"]

    nonnanidx = ~np.isnan(total_bits_alpha)
    axes[1].plot(alphas[nonnanidx], total_bits_alpha[nonnanidx], "o-", label=plotlabel)

    # for the third column, plot the pos95 vs total_bits_alpha (1st row) and vel95 vs total_bits_alpha (2nd row)
    axes[2].set_title("Position Error (95 perc) vs Total Bits")
    axes[2].set_xlabel("Total Bits")
    axes[2].set_ylabel("Position Error (m)")
    axes[2].grid(True)
    if plot_thresholds:
        axes[2].axhline(y=10, color="k", linestyle="--")
        axes[2].axvline(x=900, color="k", linestyle="--")

    nonnanidx = ~np.isnan(pos_p95_alpha_s) & ~np.isnan(total_bits_alpha)
    axes[2].plot(total_bits_alpha[nonnanidx], pos_p95_alpha[nonnanidx], "o-", label=plotlabel)
    # # plot gamma as text on the plot
    # for i, alpha in enumerate(alphas):
    #     axes[2].text(
    #         total_bits_alpha[i] + 12,
    #         pos_p95_alpha[i],
    #         rf"$\gamma$ = {alpha:.3f}",
    #         fontsize=10,
    #         ha="center",
    #         va="bottom",
    #     )

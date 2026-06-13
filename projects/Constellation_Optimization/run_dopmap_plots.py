import numpy as np
import matplotlib.pyplot as plt
import pylupnt as pnt
from src.gridsearch import compute_dop_map
import pickle
from src.postprocess import plot_user_positions
import os
from src.constellation_design import compute_enu_matrices, fibonacci_sphere


def simulate_dop_map(
    sma,
    inc,
    walker_pattern,
    x_user_grid,
    lat_users_vec,
    lon_users_vec,
    enu_mats_grid,
    dt_sim,
    use_hybrid,
    plot_ratio_under=None,
    recompute=False,
    debug=False,
):

    label = "ma{0:.0f}_inc{1:.1f}_walker{2:d}x{3:d}_hybrid{4:d}".format(
        sma, np.rad2deg(inc), walker_pattern[0], walker_pattern[1], int(use_hybrid)
    )
    pickle_file = "data/gridsearch_coverage/dopsim_{0}.pkl".format(label)
    figdop_file = "figs/coverage/dop_map_{0}.pdf".format(label)
    figorb_file = "figs/coverage/orbits_{0}.pdf".format(label)
    figsky_file = "figs/coverage/skyplot_{0}.pdf".format(label)

    if os.path.exists(pickle_file) and not recompute:
        res = pickle.load(open(pickle_file, "rb"))
        print("Loaded existing results from", pickle_file)
    else:
        res = None

    figdop, figorb, figsky, res = compute_dop_map(
        sma,
        inc,
        walker_pattern,
        et0,
        dt_sim,
        x_user_grid,
        lat_users_vec,
        lon_users_vec,
        enu_mats_grid,
        min_elev_deg=5.0,
        parallel=False,
        use_hybrid=use_hybrid,
        plot_fig=True,
        sat_fault=False,
        res=res,
        plot_ratio_under=plot_ratio_under,
        debug=debug,
    )
    figdop.savefig(figdop_file, dpi=300)
    if figorb is not None:
        figorb.write_image(figorb_file)
    if figsky is not None:
        figsky.savefig(figsky_file, dpi=300)

    with open(pickle_file, "wb") as f:
        pickle.dump(res, f)
    print("Saved results to", pickle_file)


# Main script
if __name__ == "__main__":
    # Simulation parameters
    et0 = pnt.convert_time(pnt.gregorian_to_time(2030, 1, 1, 12, 0, 0), pnt.UTC, pnt.TAI)
    dt_sim = 300.0
    sma_24hr = (pnt.GM_MOON * (24 * 3600) ** 2 / (4 * np.pi**2)) ** (1 / 3)

    # Generate user points
    # precompute user positions and ENU matrices for DOP calculation
    # Users & ENU matrices (shared, read-only)
    lat_users_vec = np.linspace(-np.pi / 2, np.pi / 2, 30 + 1)
    lon_users_vec = np.linspace(-np.pi, np.pi, 60 + 1)
    lat_grid, lon_grid = np.meshgrid(lat_users_vec, lon_users_vec)
    x_user_grid = np.zeros((lat_grid.size, 3))
    x_user_grid[:, 0] = pnt.R_MOON * np.cos(lat_grid.ravel()) * np.cos(lon_grid.ravel())
    x_user_grid[:, 1] = pnt.R_MOON * np.cos(lat_grid.ravel()) * np.sin(lon_grid.ravel())
    x_user_grid[:, 2] = pnt.R_MOON * np.sin(lat_grid.ravel())
    #
    enu_mats_grid = compute_enu_matrices(x_user_grid)  # (n_users, 3, 3)
    # plot_user_positions(x_user_grid, lat_range=(-90, 90), marker_size=2)

    # Case 1: DOP Map for South ELFO families
    sma = sma_24hr  # km
    inc = np.deg2rad(55)
    walker_pattern = (3, 3)  # total 9
    use_hybrid = False

    print("--------------------------------------")
    print(
        "Case 1: South ELFO families, sma={:.1f}km, inc={:.1f}deg, walker {}x{}, hybrid={}".format(
            sma, np.rad2deg(inc), walker_pattern[0], walker_pattern[1], int(use_hybrid)
        )
    )
    print("---------------------------------------")
    simulate_dop_map(
        sma,
        inc,
        walker_pattern,
        x_user_grid,
        lat_users_vec,
        lon_users_vec,
        enu_mats_grid,
        dt_sim,
        use_hybrid,
        plot_ratio_under=6,
        recompute=False,
        debug=False,
    )

    # Case 2: DOP Map for South ELFO families with hybrid method
    walker_pattern = (2, 4)  # total 16
    use_hybrid = True
    print("--------------------------------------")
    print(
        "Case 2: South ELFO families with hybrid method, sma={:.1f}km, inc={:.1f}deg, walker {}x{}, hybrid={}".format(
            sma, np.rad2deg(inc), walker_pattern[0], walker_pattern[1], int(use_hybrid)
        )
    )
    print("---------------------------------------")
    simulate_dop_map(
        sma,
        inc,
        walker_pattern,
        x_user_grid,
        lat_users_vec,
        lon_users_vec,
        enu_mats_grid,
        dt_sim,
        use_hybrid,
        plot_ratio_under=6,
        recompute=False,
        debug=True,
    )

    # Case 3: Circular Walker pattern with 16 satellites, inc=90 deg
    eps = 1e-6
    inc = np.arccos(np.sqrt(3 / 5)) + eps  # critical inclination
    ecc = np.sqrt(1 - 5 / 3 * np.cos(inc) ** 2)
    print("Eccentricity for inc={0:.1f} deg: {1:.4f}".format(np.rad2deg(inc), ecc))
    walker_pattern = (4, 4)
    use_hybrid = False

    print("--------------------------------------")
    print(
        "Case 3: Circular Walker pattern with 16 satellites, sma={:.1f}km, inc={:.1f}deg, walker {}x{}, hybrid={}".format(
            sma, np.rad2deg(inc), walker_pattern[0], walker_pattern[1], int(use_hybrid)
        )
    )
    print("---------------------------------------")
    simulate_dop_map(
        sma,
        inc,
        walker_pattern,
        x_user_grid,
        lat_users_vec,
        lon_users_vec,
        enu_mats_grid,
        dt_sim,
        use_hybrid,
        plot_ratio_under=6,
        recompute=False,
        debug=False,
    )

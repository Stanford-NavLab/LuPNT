from src.odsim import fit_oderr_params_parallel
import os
import numpy as np
import pylupnt as pnt

if __name__ == "__main__":
    et0 = pnt.convert_time(pnt.gregorian_to_time(2030, 1, 1, 12, 0, 0), pnt.UTC, pnt.TAI)

    # constants
    n_sma = 14
    n_ecc = 15
    n_omega = 37
    n_orbit = 4  # number of orbit to simulate
    use_int_sma = True  # whether to use integer sma values

    # results
    results = np.ones((n_sma, n_ecc, n_omega, 8)) * np.nan

    stats_ratio = (n_orbit - 1) / n_orbit  # ratio of the stats period

    # initial conditions
    if use_int_sma:
        T_ratios = np.array(
            [1 / 1, 2 / 3, 1 / 2, 2 / 5, 1 / 3, 2 / 7, 1 / 4]
        )  # integer multiples of 48 hr
        T_orbit = 48 * 3600 * T_ratios
        smas = ((T_orbit / (2 * np.pi)) ** 2 * pnt.GM_MOON) ** (1 / 3)
    else:
        smas = np.linspace(3000, 16000, n_sma)
    eccs = np.linspace(0.0, 0.70, n_ecc)
    omegas = np.linspace(0, 2 * np.pi, n_omega)

    # make data directory
    os.makedirs("data/gridsearch_od", exist_ok=True)

    if use_int_sma:
        filename = "data/gridsearch_od/fit_oderr_intsma.npy"
    else:
        filename = "data/gridsearch_od/fit_oderr.npy"

    if os.path.exists(filename):
        results = np.load(filename)
        print("Loaded existing results from", filename)
    else:
        results = fit_oderr_params_parallel(et0, smas, eccs, omegas, n_orbit=4, max_workers=10)

        # save results
        if use_int_sma:
            np.save("data/gridsearch_od/fit_oderr_intsma.npy", results)
            np.save("data/gridsearch_od/smas_ints.npy", np.arange(len(smas)))
            np.save("data/gridsearch_od/eccs.npy", eccs)
            np.save("data/gridsearch_od/omegas.npy", omegas)
        else:
            np.save("data/gridsearch_od/fit_oderr.npy", results)
            np.save("data/gridsearch_od/smas.npy", smas)
            np.save("data/gridsearch_od/eccs.npy", eccs)
            np.save("data/gridsearch_od/omegas.npy", omegas)

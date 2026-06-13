import numpy as np
from joblib import Parallel, delayed
import pylupnt as pnt
from src.orbit_manager import OrbitManager


def mc_worker(orbm, base_rv, prop_sec, t_epoch, pos_sigma, vel_sigma, seed=None):
    """
    One Monte Carlo draw: set epoch, add RV noise, propagate, return (nprop, 6).
    """
    # Reproducible randomness per worker if seed provided
    dyn = orbm.generate_dynamics()

    rng = np.random.default_rng(seed)
    d_rv = np.zeros(6)
    d_rv[:3] = rng.normal(0.0, pos_sigma, 3)
    d_rv[3:] = rng.normal(0.0, vel_sigma, 3)

    # IMPORTANT: reset epoch for each propagation
    pnt.set_lupnt_epoch(t_epoch)

    return dyn.propagate(base_rv + d_rv, prop_sec)


# Optional: independent seeds for each MC run for reproducibility
def run_pertubation_analysis(
    start_M, prop_mins, init_pos_errs, init_vel_errs, n_MC, orbm: OrbitManager, n_jobs=-1
):
    """
    Run Monte Carlo error growth analysis for the given OrbitManager instance.
    The instance should have its orbit already propagated and available as orbm.rv_prop_mci.
    """
    base_seed = 12345
    seeds = np.random.SeedSequence(base_seed).spawn(n_MC)

    # --------------------------
    # Pre-allocations (unchanged)
    # --------------------------
    nM = len(start_M)
    nprop = len(prop_mins)
    npose = len(init_pos_errs)

    prop_rv_true = np.zeros((nM, nprop, 6))
    prop_rv_pert = np.zeros((npose, n_MC, nM, nprop, 6))
    error_rnorm = np.zeros((npose, n_MC, nM, nprop))
    error_vnorm = np.zeros((npose, n_MC, nM, nprop))

    start_rv = np.zeros((nM, 6))
    t_tai = np.zeros(nM)

    # If thread-safe, keep one dynamics object; if not, create inside worker
    dyn = orbm.generate_dynamics()
    prop_sec_all = prop_mins * 60  # seconds

    for i, M in enumerate(start_M):
        print(f"Propagating M = {M:.1f} deg ...")
        idx = np.argmin(np.abs(orbm.M_prop_stack - M))
        start_rv[i, :] = orbm.rv_prop_mci[idx, :]
        t_tai[i] = orbm.t_tai[idx]

        # First propagate without error
        pnt.set_lupnt_epoch(t_tai[i])
        prop_rv_true[i, :, :] = dyn.propagate(start_rv[i, :], prop_sec_all)  # (nprop, 6)

        for j, init_pos_err in enumerate(init_pos_errs):
            print(
                f"  Init pos err = {init_pos_err:.1f} m, vel err = {init_vel_errs[j]:.4f} m/s ..."
            )
            # Parallelize over Monte Carlo runs
            results = Parallel(n_jobs=n_jobs)(
                delayed(mc_worker)(
                    orbm,
                    base_rv=start_rv[i, :],
                    prop_sec=prop_sec_all,
                    t_epoch=t_tai[i],
                    pos_sigma=init_pos_err,
                    vel_sigma=init_vel_errs[j],
                    seed=int(seeds[k].generate_state(1)[0]),  # per-run seed
                )
                for k in range(n_MC)
            )
            # Stack to (n_MC, nprop, 6)
            results = np.asarray(results)

            # Store into arrays
            prop_rv_pert[j, :, i, :, :] = results  # (n_MC, nprop, 6)

            # Vectorized error norms against truth (broadcast truth over MC axis)
            dr = results[:, :, :3] - prop_rv_true[i, np.newaxis, :, :3]  # (n_MC, nprop, 3)
            dv = results[:, :, 3:] - prop_rv_true[i, np.newaxis, :, 3:]  # (n_MC, nprop, 3)

            error_rnorm[j, :, i, :] = np.linalg.norm(dr, axis=2)  # (n_MC, nprop)
            error_vnorm[j, :, i, :] = np.linalg.norm(dv, axis=2)  # (n_MC, nprop)

    results_all = {
        "start_M": start_M,
        "prop_mins": prop_mins,
        "init_pos_errs": init_pos_errs,
        "init_vel_errs": init_vel_errs,
        "n_MC": n_MC,
        "error_rnorm": error_rnorm,
        "error_vnorm": error_vnorm,
        "prop_rv_true": prop_rv_true,
        "prop_rv_pert": prop_rv_pert,
    }

    return results_all

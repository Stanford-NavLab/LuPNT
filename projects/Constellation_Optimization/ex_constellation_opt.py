import os
import time
import pickle
import numpy as np
import multiprocessing

import pylupnt as pnt
from src.failure_model import NavSatFailureModel
from src.constellation_design import setup_problem_config, ConstellationOptimization
from src.init_population import create_initial_population

from pymoo.termination import get_termination
from pymoo.optimize import minimize
from pymoo.core.mixed import MixedVariableGA, MixedVariableSampling
from pymoo.algorithms.moo.nsga2 import RankAndCrowdingSurvival
from pymoo.core.problem import StarmapParallelization
from pymoo.visualization.scatter import Scatter
from src.postprocess import plot_history_hv
from pymoo.indicators.hv import Hypervolume
from pymoo.core.population import Population
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting


class SaveCallback:
    def __init__(self, save_dir, n_obj):
        self.save_dir = save_dir
        self.n_obj = n_obj
        self.prev_time = time.time()
        self.F_all = []
        self.fronts = None

    def __call__(self, algorithm):
        # Save entire population each generation
        pop = algorithm.pop
        gen = algorithm.n_gen

        X = pop.get("X")
        F = pop.get("F")
        n_obj = self.n_obj

        metric = Hypervolume(
            ref_point=np.ones(n_obj),
            norm_ref_point=False,
            zero_to_one=True,
            ideal=np.zeros(n_obj),
            nadir=np.ones(n_obj),
        )

        # accumulate all F values
        self.F_all.extend(F.tolist())

        # compute non-dominated front
        nds = NonDominatedSorting()

        if self.fronts is not None:
            # combine previous fronts with current F
            combined_F = np.vstack((np.array(self.F_all)[self.fronts], F))
            fronts = nds.do(combined_F, only_non_dominated_front=True)
        else:
            fronts = nds.do(np.array(self.F_all), only_non_dominated_front=True)

        # store fronts for hv calculation
        self.fronts = fronts

        # compute hv
        hv = metric.do(np.array(self.F_all)[fronts])

        data = {
            "X": pop.get("X"),
            "F": pop.get("F"),
            "F_all": np.array(self.F_all),
            "fronts": fronts,
            "hv": hv,
        }

        new_time = time.time()
        elapsed = new_time - self.prev_time

        elapsed_min = int(elapsed // 60)
        elapsed_sec = elapsed % 60

        print(
            f"Generation {gen}: {X.shape[0]} individuals, HV = {data['hv']}, Time since last gen: {elapsed_min} min {elapsed_sec:.2f} sec"
        )
        self.prev_time = new_time

        filename = os.path.join(self.save_dir, f"gen_{gen:d}.npz")
        np.savez(filename, **data)
        print(f"Saved generation {gen} data to {filename}")


def main(
    n_processes=6,
    parallel=True,
    n_gen=1,
    phases=[1, 2, 3],
    pop_size=100,
    use_float_sma=False,
    fmodel_type="long",
):
    # ------------------- config -------------------
    verbose = True
    objs = ["dop", "nsat"]
    et0 = pnt.convert_time(pnt.gregorian_to_time(2030, 1, 1, 12, 0, 0), pnt.UTC, pnt.TAI)
    # tspan = np.linspace(0, sim_days * pnt.SECS_DAY, sim_days * 24 * steps_per_hr + 1)
    dt = 15 * 60  # [s]

    if fmodel_type == "short":
        failmodel = NavSatFailureModel(
            tau1=0.25,
            S_tau1=0.95,
            beta0=0.6,
            lam1=0.008,
            tau2=5.0,
            targets_years=(6.0, 8.0),
            targets_survival=(0.8, 0.5),
        )
        max_fail = 2
    elif fmodel_type == "long":
        failmodel = NavSatFailureModel(
            tau1=0.25,
            S_tau1=0.99,
            beta0=0.6,
            lam1=0.008,
            tau2=10.0,
            targets_years=(12.0, 15.0),
            targets_survival=(0.9, 0.75),
        )
        max_fail = 1
    else:
        raise ValueError("invalid failure model type")

    phases_str = "_".join(map(str, phases))

    basedir = "/Users/keidaiiiyama/Dropbox/Research/ResearchPapers/2025/25_09_IONGNSS/ConstellationDesign/Data/"
    datadir = basedir + f"phase_{phases_str}_gen_{n_gen}_pop_{pop_size}_fmodel_{fmodel_type}"
    os.makedirs(datadir, exist_ok=True)

    if phases == [1, 2, 3]:
        sat_range_phases = [[3, 6], [6, 20], [20, 35]]
        config = setup_problem_config(
            objs,
            n_walker=3,
            n_phase=3,
            sat_range_phases=sat_range_phases,
            float_sma=use_float_sma,
            n_users=200,
            elev_mask_deg=5.0,
            cn0_thresh_user=30.0,
            dop_type_phase=["hdop", "pdop", "pdop"],
            use_wdop=True,
            dop_target=[12.0, 10.0, 8.0],
            compute_link_budget=True,
            lat_masks=[[-90, -70], [-90, -30], [-90, 90]],
            launch_years=[0, 5, 10],
            eval_years=[1, 6, 11],
            fail_model=failmodel,
            max_fail_sat=max_fail,
            epoch=et0,
            dt=dt,
        )
    else:
        raise ValueError("invalid phase configuration")

    # savefile ---------------------------------------------------------------------------
    # save the config
    with open(os.path.join(datadir, "config.pkl"), "wb") as f:
        pickle.dump(config, f)

    # callback to save each generation data
    callback = SaveCallback(save_dir=datadir, n_obj=int(len(objs) * len(phases)))

    # ------------------- problem / parallel -------------------
    ctx = multiprocessing.get_context("spawn")  # explicit for macOS / conda

    if parallel:
        with ctx.Pool(processes=n_processes) as pool:
            runner = StarmapParallelization(pool.starmap)
            problem = ConstellationOptimization(
                config=config, verbose=verbose, elementwise_runner=runner
            )

            if len(phases) == 3:
                X0 = create_initial_population(problem, sat_range_phases, pop_size, use_float_sma)
                pop = Population.new("X", X0)
                algorithm = MixedVariableGA(
                    pop_size=pop_size, survival=RankAndCrowdingSurvival(), sampling=pop
                )
            else:
                algorithm = MixedVariableGA(pop_size=pop_size, survival=RankAndCrowdingSurvival())

            termination = get_termination("n_gen", n_gen)
            start_time = time.time()
            print("start optimization...")
            res = minimize(
                problem,
                algorithm,
                termination,
                seed=1,
                save_history=False,
                callback=callback,
                verbose=True,
            )
            print(f"Optimization completed in {time.time() - start_time:.2f} seconds")
            # leaving the with-block cleanly closes & joins the pool
    else:
        problem = ConstellationOptimization(config=config, verbose=verbose)
        if len(phases) == 3:
            X0 = create_initial_population(problem, sat_range_phases, pop_size, use_float_sma)
            pop = Population.new("X", X0)
            algorithm = MixedVariableGA(
                pop_size=pop_size, survival=RankAndCrowdingSurvival(), sampling=pop
            )
        else:
            algorithm = MixedVariableGA(pop_size=pop_size, survival=RankAndCrowdingSurvival())
        termination = get_termination("n_gen", n_gen)
        start_time = time.time()
        print("start optimization...")
        res = minimize(
            problem,
            algorithm,
            termination,
            seed=1,
            save_history=False,
            callback=callback,
            verbose=True,
        )
        print(f"Optimization completed in {time.time() - start_time:.2f} seconds")


if __name__ == "__main__":
    # Important on macOS / conda:
    multiprocessing.set_start_method("spawn", force=True)

    # main(n_processes=6, parallel=True, n_gen=30, phases=[1, 2, 3], pop_size=100, fmodel_type="long")

    main(
        n_processes=10,
        parallel=True,
        n_gen=101,
        phases=[1, 2, 3],
        pop_size=100,
        use_float_sma=False,
        fmodel_type="short",
    )

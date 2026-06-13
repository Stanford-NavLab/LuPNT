import numpy as np
from pymoo.core.problem import StarmapParallelization
from pymoo.core.problem import ElementwiseProblem
import multiprocessing
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.sampling.rnd import FloatRandomSampling
from pymoo.termination import get_termination
from pymoo.optimize import minimize
import time


class MyProblem(ElementwiseProblem):

    def __init__(self, **kwargs):
        super().__init__(
            n_var=2, n_obj=2, n_ieq_constr=2, xl=np.array([-2, -2]), xu=np.array([2, 2]), **kwargs
        )

    def _evaluate(self, x, out, *args, **kwargs):
        f1 = 100 * (x[0] ** 2 + x[1] ** 2)
        f2 = (x[0] - 1) ** 2 + x[1] ** 2

        g1 = 2 * (x[0] - 0.1) * (x[0] - 0.9) / 0.18
        g2 = -20 * (x[0] - 0.4) * (x[0] - 0.6) / 4.8

        # pause for 0.01 s
        time.sleep(0.01)  # to measure the effect of parallelization

        out["F"] = [f1, f2]
        out["G"] = [g1, g2]


if __name__ == "__main__":
    # Define the number of processes and create a multiprocessing pool
    n_proccess = 10
    parallel = True

    if parallel:
        pool = multiprocessing.Pool(n_proccess)
        runner = StarmapParallelization(pool.starmap)
        problem = MyProblem(elementwise_runner=runner)
    else:
        problem = MyProblem()

    algorithm = NSGA2(
        pop_size=100,
        n_offsprings=10,
        sampling=FloatRandomSampling(),
        crossover=SBX(prob=0.9, eta=15),
        mutation=PM(eta=20),
        eliminate_duplicates=True,
    )

    termination = get_termination("n_gen", 40)

    start_time = time.time()
    res = minimize(problem, algorithm, termination, seed=1, save_history=True, verbose=True)

    end_time = time.time()
    print(f"Optimization completed in {end_time - start_time:.2f} seconds")

    X = res.X
    F = res.F

import numpy as np
from tqdm.notebook import tqdm
from multiprocessing import Pool
import hashlib
from problem import PntSchedulingProblem, State, Action
from solvers import Solver
from itertools import product
import os
import pickle


def load_or_recompute(filepath: os.PathLike, func, *args, **kwargs):
    # Check if directory exists
    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    # Check if file exists
    if os.path.exists(filepath):
        with open(filepath, "rb") as f:
            return pickle.load(f)
    else:
        result = func(*args, **kwargs)
        with open(filepath, "wb") as f:
            pickle.dump(result, f)
        return result


def iterate_params(params: dict) -> list[dict]:
    return [dict(zip(params.keys(), x)) for x in product(*params.values())]


def get_metrics(
    problem: PntSchedulingProblem,
    policy: list[tuple[State, Action]],
    time: float,
    gamma: float,
) -> dict:
    # Metrics
    percentages = [round(x, 4) for x in problem.percentage_completed(policy)]
    percentage = round(problem.duration_fullfilled(policy), 4)
    reward = round(problem.total_reward(policy, gamma=gamma), 4)
    time = round(time, 5)
    return {"total": percentage, "user": percentages, "reward": reward, "time": time}


def timed(func, *args, **kwargs):
    import time

    start = time.time()
    result = func(*args, **kwargs)
    end = time.time()
    return result, end - start


def solve(
    params: dict, problem: PntSchedulingProblem, solver: Solver, seed: int
) -> dict:
    np.random.seed(seed)
    s = problem.initial_state()
    solver_ = solver(problem)
    try:
        policy, time = timed(solver_.solve, s, **params)
    except Exception as e:
        policy, time = [(s, None)], 0

    metric = get_metrics(problem, policy, time, gamma=1)
    metric["params"] = params
    return metric


def solve_online(
    params: dict,
    problem: PntSchedulingProblem,
    solver: Solver,
    seed: int,
    return_all=False,
) -> dict:
    np.random.seed(seed)
    solver = solver(problem)

    old_policy = None
    problem.set_current_policy(None)

    if return_all:
        policies, metrics = [], []
    for ta in problem.get_arrival_times():

        # Set current time and policy
        problem.set_current_time(ta)
        if old_policy is not None:
            problem.set_current_policy(old_policy)

        # Solve
        s = problem.initial_state()
        policy, time = timed(solver.solve, s, **params)
        policy = policy if old_policy is None else problem.merge_policies(policy)
        metric = get_metrics(problem, policy, time, gamma=1)
        metric["params"] = params
        old_policy = policy

        # Save results
        if return_all:
            policies.append(policy)
            metrics.append(metric)

    if return_all:
        return metrics, policies
    return metric, policy


def run_parallel(func, it, n_jobs=4, desc="Multiprocessing", progress=True):
    if n_jobs > 1:
        with Pool(n_jobs) as p:
            if progress:
                results = list(tqdm(p.imap(func, it), total=len(it), desc=desc))
            else:
                results = p.map(func, it)
    else:
        results = [func(x) for x in tqdm(it, desc=desc)]

    return results


def cross_norm(a, b):
    return np.cross(a, b) / np.linalg.norm(np.cross(a, b), axis=2)[:, :, None]


def create_hash(*args):
    return hashlib.sha1(str(args).encode()).hexdigest()


def normalize(v: np.ndarray) -> np.ndarray:
    if v.ndim == 1:
        return v / np.linalg.norm(v)
    else:
        return v / np.linalg.norm(v, axis=1)[:, None]


def get_start_end_indexes(sequence: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    tmp = np.concatenate(([False], sequence.astype(bool), [False]))
    starts = np.where(np.logical_and(~tmp[:-1], tmp[1:]))[0]
    ends = np.where(np.logical_and(tmp[:-1], ~tmp[1:]))[0]
    return starts, ends


def my_draw_networkx_edge_labels(
    G,
    pos,
    edge_labels=None,
    label_pos=0.5,
    font_size=10,
    font_color="k",
    font_family="sans-serif",
    font_weight="normal",
    alpha=None,
    bbox=None,
    horizontalalignment="center",
    verticalalignment="center",
    ax=None,
    rotate=True,
    clip_on=True,
    rad=0,
):
    """Draw edge labels.

    Parameters
    ----------
    G : graph
        A networkx graph

    pos : dictionary
        A dictionary with nodes as keys and positions as values.
        Positions should be sequences of length 2.

    edge_labels : dictionary (default={})
        Edge labels in a dictionary of labels keyed by edge two-tuple.
        Only labels for the keys in the dictionary are drawn.

    label_pos : float (default=0.5)
        Position of edge label along edge (0=head, 0.5=center, 1=tail)

    font_size : int (default=10)
        Font size for text labels

    font_color : string (default='k' black)
        Font color string

    font_weight : string (default='normal')
        Font weight

    font_family : string (default='sans-serif')
        Font family

    alpha : float or None (default=None)
        The text transparency

    bbox : Matplotlib bbox, optional
        Specify text box properties (e.g. shape, color etc.) for edge labels.
        Default is {boxstyle='round', ec=(1.0, 1.0, 1.0), fc=(1.0, 1.0, 1.0)}.

    horizontalalignment : string (default='center')
        Horizontal alignment {'center', 'right', 'left'}

    verticalalignment : string (default='center')
        Vertical alignment {'center', 'top', 'bottom', 'baseline', 'center_baseline'}

    ax : Matplotlib Axes object, optional
        Draw the graph in the specified Matplotlib axes.

    rotate : bool (deafult=True)
        Rotate edge labels to lie parallel to edges

    clip_on : bool (default=True)
        Turn on clipping of edge labels at axis boundaries

    Returns
    -------
    dict
        `dict` of labels keyed by edge

    Examples
    --------
    >>> G = nx.dodecahedral_graph()
    >>> edge_labels = nx.draw_networkx_edge_labels(G, pos=nx.spring_layout(G))

    Also see the NetworkX drawing examples at
    https://networkx.org/documentation/latest/auto_examples/index.html

    See Also
    --------
    draw
    draw_networkx
    draw_networkx_nodes
    draw_networkx_edges
    draw_networkx_labels
    """
    import matplotlib.pyplot as plt
    import numpy as np

    if ax is None:
        ax = plt.gca()
    if edge_labels is None:
        labels = {(u, v): d for u, v, d in G.edges(data=True)}
    else:
        labels = edge_labels
    text_items = {}
    for (n1, n2), label in labels.items():
        (x1, y1) = pos[n1]
        (x2, y2) = pos[n2]
        (x, y) = (
            x1 * label_pos + x2 * (1.0 - label_pos),
            y1 * label_pos + y2 * (1.0 - label_pos),
        )
        pos_1 = ax.transData.transform(np.array(pos[n1]))
        pos_2 = ax.transData.transform(np.array(pos[n2]))
        linear_mid = 0.5 * pos_1 + 0.5 * pos_2
        d_pos = pos_2 - pos_1
        rotation_matrix = np.array([(0, 1), (-1, 0)])
        ctrl_1 = linear_mid + rad * rotation_matrix @ d_pos
        ctrl_mid_1 = 0.5 * pos_1 + 0.5 * ctrl_1
        ctrl_mid_2 = 0.5 * pos_2 + 0.5 * ctrl_1
        bezier_mid = 0.5 * ctrl_mid_1 + 0.5 * ctrl_mid_2
        (x, y) = ax.transData.inverted().transform(bezier_mid)

        if rotate:
            # in degrees
            angle = np.arctan2(y2 - y1, x2 - x1) / (2.0 * np.pi) * 360
            # make label orientation "right-side-up"
            if angle > 90:
                angle -= 180
            if angle < -90:
                angle += 180
            # transform data coordinate angle to screen coordinate angle
            xy = np.array((x, y))
            trans_angle = ax.transData.transform_angles(
                np.array((angle,)), xy.reshape((1, 2))
            )[0]
        else:
            trans_angle = 0.0
        # use default box of white with white border
        if bbox is None:
            bbox = dict(boxstyle="round", ec=(1.0, 1.0, 1.0), fc=(1.0, 1.0, 1.0))
        if not isinstance(label, str):
            label = str(label)  # this makes "1" and 1 labeled the same

        t = ax.text(
            x,
            y,
            label,
            size=font_size,
            color=font_color,
            family=font_family,
            weight=font_weight,
            alpha=alpha,
            horizontalalignment=horizontalalignment,
            verticalalignment=verticalalignment,
            rotation=trans_angle,
            transform=ax.transData,
            bbox=bbox,
            zorder=1,
            clip_on=clip_on,
        )
        text_items[(n1, n2)] = t

    ax.tick_params(
        axis="both",
        which="both",
        bottom=False,
        left=False,
        labelbottom=False,
        labelleft=False,
    )

    return text_items

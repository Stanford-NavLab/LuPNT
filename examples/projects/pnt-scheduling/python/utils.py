import matplotlib.pyplot as plt
from typing import List, Dict, Tuple
import numpy as np
import matplotlib.colors as mcolors
from dataclasses import dataclass, field
from enum import Enum
import hashlib

TABLEAU_COLORS = list(mcolors.TABLEAU_COLORS.values())


def create_hash(obj):
    return hashlib.sha1(str(obj).encode()).hexdigest()


def normalize(v):
    if v.ndim == 1:
        return v / np.linalg.norm(v)
    else:
        return v / np.linalg.norm(v, axis=1)[:, None]


def plot_windows(
    targets: np.array,  # Target ids
    vtws: np.array,  # Visibility time windows (start, end)
    durations: np.array,  # Task durations
    #
    task_otws: np.array,  # Opportunity time windows (start, end)
    task_idxs: np.array,  # Target ids
    #
    rewards: np.array = None,  # Rewards
    selected_tasks: np.array = None,  # Selected tasks
    plot_labels: bool = True,
):
    color = TABLEAU_COLORS[0]
    if selected_tasks is None:
        selected_targets = None
    else:
        assert selected_tasks.dtype == bool
        assert task_idxs.dtype == int
        selected_targets = task_idxs[selected_tasks]

    # Sort tasks by start time
    indices = sorted(range(len(vtws)), key=lambda i: vtws[i][0])

    # Add tasks to the minimum number of rows for plotting
    rows = []
    for idx in indices:
        for row in rows:
            if vtws[idx][0] > vtws[row[-1]][1]:
                row.append(idx)
                break
        else:
            rows.append([idx])

    # Plot tasks
    for r, row in enumerate(rows):

        def add_text(ts, te, txt):
            if plot_labels:
                plt.text((ts + te) / 2, r + 0.15, txt, ha="center", va="center")

        def add_rect(ts, te):
            eps = 0
            plt.fill_between(
                [ts + eps, te - eps], r - 0.3, r + 0.3, alpha=0.5, color=color
            )

        def add_empty_rect(ts, te):
            plt.fill_between(
                [ts, te],
                r - 0.3,
                r + 0.3,
                alpha=0.15,
                hatch="/",
                edgecolor="black",
                color=color,
            )

        for idx in row:
            # Opportunity window
            vtw_s = vtws[idx][0]
            vtw_e = vtws[idx][1]
            dur = durations[idx]

            if True:
                # Plot lines for visibility time windows
                plt.hlines(r, vtw_s, vtw_e, colors="k")
                plt.vlines(vtw_s, r - 0.15, r + 0.15, colors="k")
                plt.vlines(vtw_e, r - 0.15, r + 0.15, colors="k")

                # Plot lines for opportunity time windows
                otws = task_otws[task_idxs == idx]
                if len(otws) > 0:
                    for otw in otws:
                        plt.vlines(otw[0], r - 0.05, r + 0.05, colors="k")

                # Plot rectangles for tasks
                if selected_targets is None:
                    # Start times not resolved
                    ts = (vtw_s + vtw_e - dur) / 2
                    te = (vtw_s + vtw_e + dur) / 2
                    add_rect(ts, te)
                    add_text(vtw_s, vtw_e, str(targets[idx]))
                else:

                    if idx not in selected_targets:
                        # Start times resolved and task is not selected
                        ts = (vtw_s + vtw_e - dur) / 2
                        te = (vtw_s + vtw_e + dur) / 2
                        add_empty_rect(ts, te)
                        add_text(vtw_s, vtw_e, str(targets[idx]))
                    else:
                        # Start times resolved and task is selected
                        for ts, te in task_otws[selected_tasks & (task_idxs == idx)]:
                            add_rect(ts, te)
                            add_text(ts, te, str(targets[idx]))

            else:
                # Sun pointing or downlink opportunity
                if start_end_times is None:
                    add_rect(vtw_start, t_end)
                    add_text(vtw_start, t_end, str(opp.id))
                elif opp.id in start_end_times:
                    for ts, te in zip(*start_end_times[opp.id]):
                        add_rect(ts, te)
                        add_text(ts, te, str(opp.id))
                        add_text(ts, te, str(opp.id))

    plt.xlabel("Time")
    plt.xlim(0, np.max(vtws[:, 1]))
    # plt.gca().spines[["left", "top", "right"]].set_visible(False)
    # if start_end_times is not None:
    #     legend = [
    #         plt.Line2D(
    #             [0], [0], color=COLORS[tt], lw=5, label=MODE_NAMES[tt], alpha=0.5
    #         )
    #         for tt in TaskType
    #         if tt not in [TaskType.START]
    #     ]
    #     # Legend with 3 columns at the bottom outside the plot
    #     plt.legend(
    #         handles=legend,
    #         loc="upper center",
    #         bbox_to_anchor=(0.5, 1.3),
    #         ncol=3,
    #     )
    plt.yticks([])
    plt.ylabel("Tasks")


def get_start_end_indexes(sequence: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
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

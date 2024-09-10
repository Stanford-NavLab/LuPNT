import os
import pickle
import numpy as np
import time
import matplotlib.pyplot as plt
from h5py import File

# import array_to_latex as a2l
from matplotlib.gridspec import GridSpec

LUPNT_DATA_PATH = os.getenv("LUPNT_DATA_PATH")
assert LUPNT_DATA_PATH is not None, "Environment variable LUPNT_DATA_PATH not set"

LUPNT_OUTPUT_PATH = os.getenv("LUPNT_OUTPUT_PATH")


def get_output_path(output_dirs=None):
    if LUPNT_OUTPUT_PATH is None:
        output_path = os.path.join(LUPNT_DATA_PATH, "output")
    else:
        output_path = LUPNT_OUTPUT_PATH
    if output_dirs is not None:
        output_path = os.path.join(output_path, output_dirs)
    os.makedirs(output_path, exist_ok=True)
    return output_path


def normalize(vec):
    return vec / np.linalg.norm(vec)


def find_file(filename, path=LUPNT_DATA_PATH):
    for root, dirs, files in os.walk(path):
        if filename in files:
            return os.path.join(root, filename)
    return None


def timed(func, *args, **kwargs):
    start = time.time()
    result = func(*args, **kwargs)
    end = time.time()
    return result, end - start


# Load data
def load_data(directory):
    data = {}
    output_path = os.path.join(LUPNT_DATA_PATH, "output", directory)
    for filename in os.listdir(output_path):
        name, ext = os.path.splitext(filename)
        if ext == ".csv":
            try:
                data[name] = np.loadtxt(
                    os.path.join(output_path, filename), delimiter=","
                )
            except ValueError:
                # Each row has a different number of elements
                tmp = []
                with open(os.path.join(output_path, filename), "r") as f:
                    for line in f:
                        tmp.append([float(x) for x in line.strip().split(",")])
                max_len = max([len(x) for x in tmp])
                data[name] = np.array([x + [np.nan] * (max_len - len(x)) for x in tmp])
    return data


def timer_func(func):
    # This function shows the execution time of
    # the function object passed
    def wrap_func(*args, **kwargs):
        t1 = time()
        result = func(*args, **kwargs)
        t2 = time()
        print(f"Function {func.__name__!r} executed in {(t2-t1):.4f}s")
        return result

    return wrap_func


def get_basepath():
    return LUPNT_DATA_PATH


def load_pickle(path):
    with open(path, "rb") as f:
        return pickle.load(f)


def dump_pickle(obj, path):
    with open(path, "wb") as f:
        pickle.dump(obj, f)


def set_axes_equal(ax):
    """Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    """

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


# def printState(df, idx):
#     row = df.iloc[idx]
#     print(row.iloc[0])
#     vec = row.iloc[1:].to_numpy().reshape(-1, 1)
#     a2l.to_ltx(vec, frmt="{:g}", arraytype="bmatrix", print_out=True)
#     return vec


# def printPosVel(df, idx):
#     row = df.iloc[idx]
#     print(row.iloc[0])
#     a2l.to_ltx(
#         row.iloc[1:4].to_numpy().reshape(-1, 1),
#         frmt="{:g}",
#         arraytype="bmatrix",
#         print_out=True,
#     )
#     a2l.to_ltx(
#         row.iloc[4:].to_numpy().reshape(-1, 1),
#         frmt="{:g}",
#         arraytype="bmatrix",
#         print_out=True,
#     )


def plot_RTN(
    rv_RTN, labels=None, legend_text=None, init=False, final=True, center=True
):
    delta = 5
    # Find limits
    if labels is None:
        labels = [""]
    if legend_text is None:
        legend_text = labels

    limits = {
        "rR": [
            min(np.min(rv_RTN[f"rR{s}"]) for s in labels) - delta,
            max(np.max(rv_RTN[f"rR{s}"]) for s in labels) + delta,
        ],
        "rT": [
            min(np.min(rv_RTN[f"rT{s}"]) for s in labels) - delta,
            max(np.max(rv_RTN[f"rT{s}"]) for s in labels) + delta,
        ],
        "rN": [
            min(np.min(rv_RTN[f"rN{s}"]) for s in labels) - delta,
            max(np.max(rv_RTN[f"rN{s}"]) for s in labels) + delta,
        ],
    }

    # Create figure
    wr = (limits["rT"][1] - limits["rT"][0]) / (limits["rR"][1] - limits["rR"][0])
    gs = GridSpec(
        2, 2, height_ratios=[1, 1], width_ratios=[wr, 1], hspace=0.25, wspace=0.25
    )
    fig = plt.figure(figsize=(4.5 * (wr + 1), 9))

    # Plotting function
    def plot_planes(ax, p1, p2):
        for s, lt in zip(labels, legend_text):
            ax.plot(rv_RTN[f"r{p1}{s}"], rv_RTN[f"r{p2}{s}"], lw=0.5, label=lt)
            if init:
                ax.plot(
                    rv_RTN[f"r{p1}{s}"].iloc[0],
                    rv_RTN[f"r{p2}{s}"].iloc[0],
                    "o",
                    label=f"Deputy {lt}",
                )
            if final:
                ax.plot(
                    rv_RTN[f"r{p1}{s}"].iloc[-1],
                    rv_RTN[f"r{p2}{s}"].iloc[-1],
                    "o",
                    label=f"Final {lt}",
                )
            if center:
                ax.plot(0, 0, "o", label=f"Chief {lt}")

        ax.set_xlabel(p1 + " [m]")
        ax.set_ylabel(p2 + " [m]")
        ax.set_xlim(limits["r" + p1])
        ax.set_ylim(limits["r" + p2])
        ax.set_aspect("equal")
        ax.grid()

    # Plot 3D trajectory
    ax = fig.add_subplot(gs[0, 1], projection="3d")
    for s, lt in zip(labels, legend_text):
        ax.plot(
            rv_RTN[f"rR{s}"],
            rv_RTN[f"rT{s}"],
            rv_RTN[f"rN{s}"],
            lw=0.5,
            label=lt,
        )

    for s, lt in zip(labels, legend_text):
        if init:
            ax.plot(
                rv_RTN[f"rR{s}"].iloc[0],
                rv_RTN[f"rT{s}"].iloc[0],
                rv_RTN[f"rN{s}"].iloc[0],
                "o",
                label=f"Initial {lt}",
            )
        if final:
            ax.plot(
                rv_RTN[f"rR{s}"].iloc[-1],
                rv_RTN[f"rT{s}"].iloc[-1],
                rv_RTN[f"rN{s}"].iloc[-1],
                "o",
                label=f"Deputy {lt}",
            )
        if center:
            ax.plot(0, 0, 0, "o", label=f"Chief {lt}")

    ax.legend(
        bbox_to_anchor=(0.5, -0.25),
        loc="lower center",
        ncol=2,
        frameon=False,
    )

    ax.set_xlabel("R [m]")
    ax.set_ylabel("T [m]")
    ax.set_zlabel("N [m]")
    ax.set_xlim3d(limits["rR"])
    ax.set_ylim3d(limits["rT"])
    ax.set_zlim3d(limits["rN"])
    # ax.view_init(elev=30, azim=-20)

    ax = fig.add_subplot(gs[0, 0])
    plot_planes(ax, "T", "N")
    ax = fig.add_subplot(gs[1, 0])
    plot_planes(ax, "T", "R")
    ax = fig.add_subplot(gs[1, 1])
    plot_planes(ax, "N", "R")


def format_element(x, fmt="{}"):
    if x == 0.0:
        return "0.0"
    else:
        return fmt.format(x)


def print_aligned(matrix):
    if len(matrix.shape) == 1:
        num_columns = 1
        matrix = matrix.reshape((len(matrix), 1))
    else:
        num_columns = len(matrix[0])
    max_before_decimal = [0] * num_columns
    max_after_decimal = [0] * num_columns

    # Step 1: Find the maximum length before and after decimal for each column
    for row in matrix:
        for i, value in enumerate(row):
            parts = str(value).split(".")
            max_before_decimal[i] = max(max_before_decimal[i], len(parts[0]))
            if len(parts) > 1:
                max_after_decimal[i] = max(max_after_decimal[i], len(parts[1]))

    # Step 2: Format each number with padding and accumulate in a string
    txt = "\n["
    for r, row in enumerate(matrix):
        line = "[" if r == 0 else " ["
        for i, value in enumerate(row):
            value = float(format_element(value))
            before, after = (
                str(value).split(".") if "." in str(value) else (str(value), "")
            )
            padding_left = " " * (max_before_decimal[i] - len(before))
            padding_right = " " * (max_after_decimal[i] - len(after))
            formatted_value = f"{padding_left}{before}.{after}{padding_right}"
            line += formatted_value + (" " if i < num_columns - 1 else "")
        txt += line + "]\n"

    print(txt)

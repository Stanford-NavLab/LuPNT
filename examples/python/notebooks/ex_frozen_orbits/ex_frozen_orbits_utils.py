# Other imports
import os
import numpy as np
import pylupnt as pnt
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import plotly.io as pio

if "VSCODE_PID" in os.environ:
    pio.renderers.default = "notebook"

cov_cmap_tmp = plt.get_cmap("Blues", 3 + 1)
cov_cmap = ListedColormap(cov_cmap_tmp([1, 2, 3]))
np.set_printoptions(precision=4, suppress=True)

# Data path
# output_path = pnt.get_output_path("ex_frozen_orbits")
output_path = os.path.join(os.path.dirname(__file__))
data_part1_path = os.path.join(output_path, "data_part1.h5")
data_part2_path = os.path.join(output_path, "data_part2.h5")
data_part1_exists = os.path.exists(data_part1_path)
data_part2_exists = os.path.exists(data_part2_path)

case0_0_text = """
<b>Satellite orbit</b><br>
Frame: Moon OP<br>
Dynamics: Two body (numerical)<br>
Propagation time: {:.2f} hours
"""

case0_0_text_args = dict(align="left", xanchor="left", x=0.01, y=0.98, showarrow=False)

case1_0_text = """
<b>Moon orbit</b><br>
Frame: Moon OP<br>
Dynamics: Two body (analytical)<br>
Propagation time: {:.2f} days
"""

case1_0_text_args = dict(align="left", xanchor="left", x=0.01, y=0.98, showarrow=False)

case1_1_text = """
<b>Satellite orbits</b><br>
Frame: Moon OP<br>
Dynamics: three body (Earth circular orbit)<br>
Propagation time: {:.2f} days
"""

case1_1_text_args = dict(align="left", xanchor="left", x=0.01, y=0.98, showarrow=False)


def plot_mean_anomaly_differences(tfs: np.ndarray, rvs: np.ndarray, set_lims=True):
    coes = np.zeros_like(rvs)
    for i in range(3):
        coes[i] = pnt.cart2classical(rvs[i], pnt.GM_MOON)

    plt.figure(figsize=(8, 4))
    plt.suptitle(
        """Mean anomaly differences between the second and first spacecraft (left)
        and the third and first spacecraft (right), with no phasing adjustments"""
    )
    # lims = [[115, 145], [-400, -100]]
    for i in range(2):
        plt.subplot(1, 2, i + 1)
        delta_M_tmp = pnt.wrap2Pi(coes[i + 1, :, 5] - coes[0, :, 5]) * pnt.DEG
        delta_M = np.zeros_like(delta_M_tmp)
        delta_M[0] = delta_M_tmp[0]
        for j in range(1, len(delta_M)):
            delta_M[j] = delta_M[j - 1] + pnt.wrap2Pi(delta_M_tmp[j] - delta_M[j - 1])
        plt.plot((tfs - tfs[0]) / pnt.SECS_DAY, delta_M)
        plt.xlabel("Days past " + pnt.time2gregorian_string(tfs[0]) + " TAI")
        plt.ylabel(f"$M_{i + 2} - M_1$ [deg]")
        plt.xlim(0, (tfs[-1] - tfs[0]) / pnt.SECS_DAY)
        if set_lims:
            if i == 0:
                plt.ylim(110, 130)
            else:
                plt.ylim(-130, -110)
        plt.grid()
    plt.tight_layout()
    plt.show()

"""Plot built-in antenna gain patterns.

This exploratory plotting example opens a Matplotlib window when run with an
interactive backend.  Use ``MPLBACKEND=Agg`` to smoke-test it without a GUI.
"""

import pylupnt as pnt
import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    antenna_names = [
        "Parabora_S_d10",
        "Parabora_S_d100",
        "Block-IIA_ACE",
        "Block-IIR-M_ACE",
        "BEIDOU_IGSO",
        "BEIDOU_MEO",
        "GALLILEO",
        "DSN-S",
        "DSN-X",
        "Patch_22_RHCP_8025MHz",
    ]

    fig = plt.figure(figsize=(15, 15))
    axs = fig.subplots(4, 4).flatten()
    for i, name in enumerate(antenna_names):
        print(f"Computing gain pattern for {name}")
        antenna = pnt.Antenna(name)
        phi = np.linspace(-400, 400, 100)  # [deg]
        theta = np.linspace(-400, 400, 5)  # [deg]
        plt.sca(axs[i])
        for az in theta:
            gain = antenna.compute_gain(az * pnt.RAD, phi * pnt.RAD)
            wrapped_phi = np.vectorize(pnt.wrap2pi)(phi * pnt.RAD)
            plt.plot(pnt.DEG * wrapped_phi, gain, label=f"theta = {az:.0f} deg")
        plt.title(name)
        plt.xlabel("Phi [deg]")
        plt.ylabel("Gain [dB]")
        plt.text(
            0.98,
            0.95,
            f"Max {antenna.get_gain_matrix().max():.2f} dB",
            transform=plt.gca().transAxes,
            ha="right",
            va="top",
        )
        if i == 0:
            plt.legend()
        plt.grid()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()

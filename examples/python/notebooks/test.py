import pylupnt as pnt
import plotly.graph_objects as go
import plotly.io as pio
import matplotlib.pyplot as plt
import numpy as np


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
    "Patch_22_RHCP_8025MHz",  # "LGPS", "moongpsr",
]

# Subplots
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
        plt.plot(pnt.DEG * pnt.wrap2pi(phi * pnt.RAD), gain, label=f"theta = {az:.0f}°")
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

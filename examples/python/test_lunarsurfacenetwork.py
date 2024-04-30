import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylupnt as pnt

import pylupnt.plots as plots
import pylupnt.utils as u

# mpl.rcParams.update(mpl.rcParamsDefault)

colors = plots.COLORS
fig = plots.Plot3D(elev=25, azim=-50)
fig.plot_surface("MOON")
fig.label_axis()
plt.title("Satellite orbit (MI)")
plt.show()

print(colors)
# print(fig)

def degToRad(deg):
    return deg * np.pi / 180


# Define satellite parameters (ELFO, Lunar Pathfinder)
a = 5760 # km
e = 0.58
i = degToRad(54.856)
Omega = 0
w = degToRad(86.322)
M = degToRad(180)
# nu = degToRad(92.335)
# M = pnt.true_to_mean_anomaly(nu, e)

oe = np.array([a, e, i, Omega, w, M])

# State
x_oe = pnt.ClassicalOE(oe, coord_sys=pnt.CoordSystem.MI)
print(" ")
print("Classical orbital elements:")
print(x_oe.vector)




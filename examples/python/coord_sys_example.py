import pylupnt as pnt
import numpy as np

# Create a coordinate converter object
rv_in = np.array([5102.5096, 6123.01152, 6378.1368, -4.7432196, 0.7905366, 5.55337561])
epoch = 2000.0

print(
    "GCRF:",
)
print(rv_in)
print(" ")

rv_out = pnt.CoordConverter.convert(epoch, rv_in, coord_sys_in=pnt.CoordSystem.GCRF, coord_sys_out=pnt.CoordSystem.ITRF)

print("ITRF:")
print(rv_out)
print(" ")

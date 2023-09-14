import pylupnt as pnt
import numpy as np

# Create a coordinate converter object
rv_in_ad = pnt.VectorXreal(6)
rv_in_vec = [5102.5096, 6123.01152, 6378.1368, -4.7432196, 0.7905366, 5.55337561]
for i in range(6):
    rv_in_ad[i] = rv_in_vec[i]

coord_in = "GCRF"
coord_out = "ITRF"
epoch = 2000.0

print("GCRF:",)
print(rv_in_ad)
print(" ")

rv_out_ad = pnt.CoordConverter.convert(rv_in_ad, epoch, coord_in, coord_out)

print("ITRF:")
print(rv_out_ad)
print(" ")
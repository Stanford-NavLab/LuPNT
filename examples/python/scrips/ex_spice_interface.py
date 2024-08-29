import pylupnt as pnt
from pylupnt import SpiceInterface as sp
import numpy as np

sp.load_spice_kernel()
tdb = sp.string2tdb("2000/01/01 12:00:00")
tai = sp.string2tai("2000/01/01 12:00:00")
utc = sp.tdb_to_string_utc(tdb, 3)
et = sp.convert_time(tdb, "TDB", "ET")
earth_et0 = sp.get_body_pos_vel(tai, 10, 399)
rot_j2i = sp.get_frame_conversion_matrix(et, pnt.GCRF, pnt.ITRF)

print("tdb at 2000/01/01 12:00:00  :", tdb)
print("tai at 2000/01/01 12:00:00  :", tai)
print("UTC at above TDB            :", utc)
print("et at above TDB             :", et)
print("Earth at et0                :")
print(earth_et0)
print("Rotation from J2000 to ITRF :")
print(rot_j2i)

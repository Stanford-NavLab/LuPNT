"""Convert a single Cartesian state between LuPNT reference frames.

Frame conversion epochs are expressed as TDB seconds from J2000.  The sample
state below uses SI units: meters and meters per second.
"""

import numpy as np
import pylupnt as pnt


def main() -> None:
    # Six-element states are ordered as [x, y, z, vx, vy, vz].
    rv_in = np.array(
        [5102.5096e3, 6123.01152e3, 6378.1368e3, -4.7432196e3, 0.7905366e3, 5.55337561e3]
    )
    t_tdb = 2000.0

    print("GCRF:")
    print(rv_in)
    print()

    rv_out = pnt.convert_frame(t_tdb, rv_in, frame_in=pnt.Frame.GCRF, frame_out=pnt.Frame.ITRF)

    print("ITRF:")
    print(rv_out)
    print()


if __name__ == "__main__":
    main()

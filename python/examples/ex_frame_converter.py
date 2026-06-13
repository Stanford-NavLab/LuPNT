"""Batch frame conversion for multiple Cartesian states.

This example is similar to ``ex_coord_sys.py`` but passes a matrix of states to
show the vectorized Python API.
"""

import numpy as np
import pylupnt as pnt


def main() -> None:
    # Each row is a six-element Cartesian state [r, v].
    rv_in = np.array(
        [
            [5102.5096e3, 6123.01152e3, 6378.1368e3, -4.7432196e3, 0.7905366e3, 5.55337561e3],
            [5102.5096e3, 6123.01152e3, 6378.1368e3, -4.7432196e3, 0.7905366e3, 5.55337561e3],
        ]
    )
    t_tdb = 2000.0

    print("GCRF:")
    print(rv_in)
    print()

    rv_out = pnt.convert_frame(t_tdb, rv_in, pnt.Frame.GCRF, pnt.Frame.ITRF)

    print("ITRF:")
    print(rv_out)
    print()


if __name__ == "__main__":
    main()

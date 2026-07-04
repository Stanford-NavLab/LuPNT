"""Convert classical orbital elements to Cartesian form.

The previous version of this example used the `ClassicalOE` state wrapper, which
is no longer exported by the current Python package.  The supported conversion
path is the NumPy-vector API shown here.
"""

import numpy as np
import pylupnt as pnt


def deg2rad(deg: float) -> float:
    return deg * np.pi / 180


def main() -> None:
    x_cart = np.array([6524.834e3, 6862.875e3, 6448.296e3, 4.901327e3, 5.533756e3, -1.976341e3])
    print()
    print("Cartesian state:")
    print(x_cart)

    p = 11067.790e3
    e = 0.83285
    a = p / (1 - pow(e, 2))
    i = deg2rad(87.87)
    Omega = deg2rad(227.89)
    w = deg2rad(53.38)
    nu = deg2rad(92.335)
    M = pnt.true2mean(nu, e)

    # Elements are [a, e, i, RAAN, argument of periapsis, mean anomaly].
    x_oe = np.array([a, e, i, Omega, w, M])
    print()
    print("Classical orbital elements:")
    print(x_oe)

    print()
    print("a = ", x_oe[0])

    # Conversion functions operate on array-like vectors.
    x_cart_from_oe = pnt.classical_to_cart(x_oe, pnt.GM_EARTH)
    print()
    print("Converted Cartesian State:")
    print(x_cart_from_oe)


if __name__ == "__main__":
    main()

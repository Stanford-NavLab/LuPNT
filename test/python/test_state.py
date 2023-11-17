import pylupnt as pnt
import numpy as np
import sys
import pytest

sys.path.append("pylupntutil")
from setup_gmat import gmat


def test_StateUtils():
    # Classical Orbital elements
    a = 6541.4
    e = 0.6
    i = np.deg2rad(65.5)
    Omega = np.deg2rad(90.0)
    w = np.deg2rad(0.0)
    M = np.deg2rad(0.0)

    coe_state = pnt.ClassicalOE(a, e, i, Omega, w, M, cs=pnt.CoordSystem.MI)
    coe_state.get_vector()
    tmp = pnt.Vector3real([1, 2, 3])

    mu = pnt.MU_MOON
    coe_vec = coe_state.get_vector()

    np.testing.assert_array_almost_equal(coe_vec, [a, e, i, Omega, w, M])

    cart_vec_gmat = gmat.StateConversionUtil.KeplerianToCartesian(
        mu, gmat.Rvector6(*coe_vec), "MA"
    )
    cart_vec_gmat = np.array([cart_vec_gmat.Get(i) for i in range(6)])
    cart_vec_gmat[2:6] = np.deg2rad(cart_vec_gmat[2:6])

    cart_vec = pnt.coe_to_cart(coe_vec, mu)
    # np.testing.assert_array_almost_equal(cart_vec_gmat, cart_vec)


if __name__ == "__main__":
    test_StateUtils()

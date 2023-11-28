import pylupnt as pnt
import numpy as np
import pytest

try:
    from . import utils
    from .gmat import gmat
except ImportError:
    import utils
    from gmat import gmat


class TestNumericalDynamics:
    def test_CartesianTwoBodyDynamics(self):
        # Constructor
        cart_array = pnt.classical_to_cartesian(utils.coe_array_elfo, pnt.MU_MOON)
        dyn = pnt.CartesianTwoBodyDynamics(pnt.MU_MOON)

        # Propagate
        t0 = 0.0
        tf = 1000.0
        dt = 1.0

        cart_array2 = dyn.propagate(cart_array, t0, tf, dt)

    def test_FullDynamics(self):
        pass


if __name__ == "__main__":
    pytest.main([__file__])

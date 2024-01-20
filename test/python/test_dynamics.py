import pylupnt as pnt
import numpy as np
import pytest

try:
    from .utils import data, gmat_helpers
    from .utils.gmat import gmat
except ImportError:
    from utils import data, gmat_helpers
    from utils.gmat import gmat


class TestNumericalOrbitDynamics:
    def test_CartesianTwoBodyDynamics(self):
        # Constructor
        cart_array = pnt.classical_to_cartesian(data.coe_array_elfo, pnt.MU_MOON)
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

import pylupnt as pnt
import numpy as np
import pytest
from ..gmat.utils import data
from test.gmat.utils import data

try:
    from .utils import gmat_helpers
    from .utils.gmat import gmat
except ImportError:
    from utils import gmat_helpers
    from utils.gmat import gmat


class TestStateUtils:
    def test_classical(self):
        # Cartesian to Classical

        cart_array = pnt.classical_to_cartesian(data.coe_array_elfo, pnt.MU_MOON)
        cart_state = pnt.CartesianOrbitState(cart_array, pnt.CoordSystem.MI)

        coe_array = pnt.cartesian_to_classical(cart_array, pnt.MU_MOON)
        coe_state = pnt.cartesian_to_classical(cart_state, pnt.MU_MOON)

        cart_gmat = cart_array
        coe_gmat = gmat_helpers.unpack_rvector(
            gmat.StateConversionUtil.CartesianToKeplerian(
                pnt.MU_MOON, gmat.Rvector6(*cart_gmat), "MA"
            )
        )
        coe_gmat[2:6] = pnt.wrapToPi(np.deg2rad(coe_gmat[2:6]))

        assert np.allclose(coe_array, coe_gmat)
        assert np.allclose(coe_state.vector, coe_gmat)

        # Classical to Cartesian

        coe_array = data.coe_array_elfo
        coe_state = pnt.ClassicalOE(coe_array, pnt.CoordSystem.MI)

        cart_array = pnt.classical_to_cartesian(coe_array, pnt.MU_MOON)
        cart_state = pnt.classical_to_cartesian(coe_state, pnt.MU_MOON)

        coe_gmat = coe_array
        coe_gmat[2:6] = np.rad2deg(coe_gmat[2:6])
        cart_gmat = gmat_helpers.unpack_rvector(
            gmat.StateConversionUtil.KeplerianToCartesian(
                pnt.MU_MOON, gmat.Rvector6(*coe_gmat), "MA"
            )
        )
        assert np.allclose(cart_array, cart_gmat)
        assert np.allclose(cart_state.vector, cart_gmat)

    def test_equinoctial(self):
        pass


if __name__ == "__main__":
    pytest.main([__file__])

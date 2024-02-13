import pylupnt as pnt
import numpy as np
import pytest

try:
    from .utils import gmat_helpers
    from .utils.gmat import gmat
except ImportError:
    from utils import gmat_helpers
    from utils.gmat import gmat


class TestCoordConverter:
    def test_conversions(self):
        cart_MI = pnt.classical_to_cartesian(data.coe_array_elfo, pnt.MU_MOON)
        epoch = pnt.SpiceInterface.string_to_tai("2020/07/20 12:00:00.000")

        # coord_systems = pnt.CoordSystem.__members__.values()
        coord_systems = [pnt.ITRF, pnt.GCRF, pnt.ICRF, pnt.MI, pnt.PA]
        coord_systems = [pnt.ITRF]

        for coord_sys_from in coord_systems:
            print("From", pnt.MI, "to", coord_sys_from)
            cart_from = pnt.CoordConverter.convert(
                epoch, cart_MI, pnt.MI, coord_sys_from
            )
            cart_from_gmat = gmat_helpers.convert_coord(
                epoch, cart_MI, pnt.MI, coord_sys_from
            )
            assert np.allclose(cart_from, cart_from_gmat)

            for coord_sys_to in coord_systems:
                print("From", coord_sys_from, "to", coord_sys_to)
                cart_to = pnt.CoordConverter.convert(
                    epoch, cart_from, coord_sys_from, coord_sys_to
                )
                cart_to_gmat = gmat_helpers.convert_coord(
                    epoch, cart_from_gmat, coord_sys_from, coord_sys_to
                )
                assert np.allclose(cart_to, cart_to_gmat)


if __name__ == "__main__":
    # pytest.main([__file__])
    TestCoordConverter().test_conversions()

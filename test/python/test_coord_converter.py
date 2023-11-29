import pylupnt as pnt
import numpy as np
import pytest

try:
    from .utils import data, gmat_helpers
    from .utils.gmat import gmat
except ImportError:
    from utils import data, gmat_helpers
    from utils.gmat import gmat


class TestCoordConverter:
    def test_conversions(self):
        cart_MI = pnt.classical_to_cartesian(data.coe_array_elfo, pnt.MU_MOON)

        # for coord_sys_from in pnt.CoordSystem.__members__.values():


if __name__ == "__main__":
    pytest.main([__file__])

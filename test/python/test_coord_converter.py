import pylupnt as pnt
import numpy as np
import pytest
import pickle
import os


def map_string_pnt_coordinate(str):
    if str == "ITRF":
        return pnt.Frame.ITRF
    elif str == "GCRF":
        return pnt.Frame.GCRF
    elif str == "ICRF":
        return pnt.Frame.ICRF
    elif str == "MI":
        return pnt.Frame.MI
    elif str == "PA":
        return pnt.Frame.PA
    else:
        raise ValueError("Unknown coordinate system: " + str)


class TestCoordConverter:
    def test_conversions(self):
        # load data
        basepath = pnt.utils.get_basepath()
        gmat_data_path = os.path.join(basepath, "..", "test", "gmat", "data")
        filename = os.path.join(gmat_data_path, "coord_conversions.pkl")
        with open(filename, "rb") as f:
            data = pickle.load(f)

        cart_from = data["cart_from"]
        coord_froms = data["coord_froms"]
        coord_tos = data["coord_tos"]
        cart_to = data["cart_to"]
        epoch = data["epoch"]

        for i, coord_sys_from in enumerate(coord_froms):
            coord_sys_from = map_string_pnt_coordinate(coord_sys_from)
            for j, coord_sys_to in enumerate(coord_tos):
                coord_sys_to = map_string_pnt_coordinate(coord_sys_to)
                print("From", coord_sys_from, "to", coord_sys_to)
                cart_to_pnt = pnt.CoordConverter.convert(
                    epoch, cart_from, coord_sys_from, coord_sys_to
                )
                cart_to_gmat = cart_to[i][j]

                assert np.allclose(cart_to_pnt, cart_to_gmat)


if __name__ == "__main__":
    # pytest.main([__file__])
    TestCoordConverter().test_conversions()

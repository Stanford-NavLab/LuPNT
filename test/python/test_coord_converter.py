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
        return pnt.Frame.MOON_CI
    elif str == "PA":
        return pnt.Frame.MOON_PA
    else:
        raise ValueError("Unknown coordinate system: " + str)


class TestFrameConverter:
    def test_conversions(self):
        # load data
        basepath = pnt.utils.get_basepath()
        gmat_data_path = os.path.join(basepath, "..", "test", "gmat", "data")
        filename = os.path.join(gmat_data_path, "frame_conversions.pkl")
        with open(filename, "rb") as f:
            data = pickle.load(f)

        cart_from = data["cart_from"]
        frame_froms = data["frame_froms"]
        frame_tos = data["frame_tos"]
        cart_to = data["cart_to"]
        epoch = data["epoch"]

        for i, frame_from in enumerate(frame_froms):
            frame_from = map_string_pnt_coordinate(frame_from)
            for j, frame_to in enumerate(frame_tos):
                frame_to = map_string_pnt_coordinate(frame_to)
                print("From", frame_from, "to", frame_to)
                cart_to_pnt = pnt.FrameConverter.convert(
                    epoch, cart_from, frame_from, frame_to
                )
                cart_to_gmat = cart_to[i][j]

                assert np.allclose(cart_to_pnt, cart_to_gmat)


if __name__ == "__main__":
    # pytest.main([__file__])
    TestFrameConverter().test_conversions()

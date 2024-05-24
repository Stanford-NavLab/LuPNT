# generate test data using gmat
import numpy as np
import pickle

try:
    from .utils import gmat_helpers
    from .utils.gmat import gmat
except ImportError:
    from utils import gmat_helpers
    from utils.gmat import gmat


def generate_test_data():
    generate_frame_conversions()


###############################################################
# Coordinate Conversion
###############################################################
def generate_frame_conversions():

    # parameter to be used
    epoch_str = "20 Jul 2020 12:00:00.000"
    cart_from = np.array(
        [5102.5096, 6123.01152, 6378.1368, -4.7432196, 0.7905366, 5.55337561]
    )

    # genererate gmat epoch
    time_sys_converter = gmat.TimeSystemConverter.Instance()
    epoch_gmat = time_sys_converter.Convert(
        time_sys_converter.ConvertGregorianToMjd(epoch_str),
        gmat.TimeSystemConverter.UTC,
        gmat.TimeSystemConverter.TAI,
    )
    print(epoch_gmat)

    # perform coordinate conversion
    frame_froms = ["GCRF"]
    frame_tos = ["ITRF", "ICRF", "MI", "PA"]

    # filename
    filename = "data/frame_conversions.pkl"
    cart_to_store = np.zeros((len(frame_froms), len(frame_tos), 6))

    for i, frame_from in enumerate(frame_froms):
        for j, frame_to in enumerate(frame_tos):
            cart_to = gmat_helpers.convert_coord(
                epoch_gmat, cart_from, frame_from, frame_to
            )
            # store in file
            cart_to_store[i, j, :] = cart_to

    # convert epoch to pylupnt epoch
    epoch_tai = gmat_helpers.convert_gmat_to_pylupnt_epoch(epoch_gmat)

    # store in file
    data = {}
    data["epoch"] = epoch_tai
    data["cart_from"] = cart_from
    data["frame_froms"] = frame_froms
    data["frame_tos"] = frame_tos
    data["cart_to"] = cart_to_store

    with open(filename, "wb") as f:
        pickle.dump(data, f)


if __name__ == "__main__":
    generate_test_data()

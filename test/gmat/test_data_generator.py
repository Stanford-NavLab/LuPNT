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
    generate_coord_conversions()


###############################################################
# Coordinate Conversion
###############################################################
def generate_coord_conversions():

    # parameter to be used
    epoch_str = "20 Jul 2020 12:00:00.000"
    cart_from = np.array([5102.5096, 6123.01152, 6378.1368, -4.7432196, 0.7905366, 5.55337561])

    # genererate gmat epoch
    time_sys_converter = gmat.TimeSystemConverter.Instance()
    epoch_gmat = (
        time_sys_converter.Convert(
            time_sys_converter.ConvertGregorianToMjd(epoch_str),
            gmat.TimeSystemConverter.UTC,
            gmat.TimeSystemConverter.TAI,
        )
    )
    print(epoch_gmat)

    # perform coordinate conversion
    coord_froms = ["GCRF"]
    coord_tos = ["ITRF", "ICRF", "MI", "PA"]

    # filename
    filename = "data/coord_conversions.pkl"
    cart_to_store = np.zeros((len(coord_froms), len(coord_tos), 6))

    for i, coord_from in enumerate(coord_froms):
        for j, coord_to in enumerate(coord_tos):
            cart_to = gmat_helpers.convert_coord(
                epoch_gmat, cart_from, coord_from, coord_to
            )
            # store in file
            cart_to_store[i, j, :] = cart_to

    # store in file
    data = {}
    data["epoch"] = epoch_gmat  
    data["cart_from"] = cart_from
    data["coord_froms"] = coord_froms
    data["coord_tos"] = coord_tos
    data["cart_to"] = cart_to_store

    with open(filename, "wb") as f:
        pickle.dump(data, f)

def load_coord_conversions():
    filename = "data/coord_conversions.pkl"
    with open(filename, "rb") as f:
        data = pickle.load(f)
    
    epoch_gmat = data["epoch"]

    # convert gmat epoch to pylupnt epoch
    data["epoch"] = gmat_helpers.convert_gmat_to_pylupnt_epoch(epoch_gmat)

    return data



if __name__ == "__main__":
    generate_test_data()


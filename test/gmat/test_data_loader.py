import pickle

def load_coord_conversions():
    filename = "data/coord_conversions.pkl"
    with open(filename, "rb") as f:
        data = pickle.load(f)
    return data


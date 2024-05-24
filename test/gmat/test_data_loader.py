import pickle


def load_frame_conversions():
    filename = "data/frame_conversions.pkl"
    with open(filename, "rb") as f:
        data = pickle.load(f)
    return data

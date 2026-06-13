import os
import pickle
import pylupnt as pnt
import pandas as pd


def pickle_dir_to_csvs(input_dir, filename_keyword="short"):
    """
    Load all pickle files containing `filename_keyword` in their name
    and save each DataFrame to a separate CSV file.
    """

    for fname in sorted(os.listdir(input_dir)):
        if filename_keyword in fname and fname.endswith(".pkl"):
            pkl_path = os.path.join(input_dir, fname)
            csv_name = fname.replace(".pkl", ".csv")
            csv_path = os.path.join(input_dir, csv_name)

            print(f"Loading: {pkl_path}")

            with open(pkl_path, "rb") as f:
                obj = pickle.load(f)

            if not isinstance(obj, pd.DataFrame):
                raise TypeError(f"{fname} does not contain a pandas DataFrame")

            print(f"Saving: {csv_path}")
            obj.to_csv(csv_path, index=False)


if __name__ == "__main__":

    n_orbit = 6
    dt = 1
    dtrt = 120

    datapath = os.path.join(pnt.get_output_dir(), "iono_delay", "labels_csv")
    orbit_dt_str = "norbit_{}_dt_{}s_dtrt_{}s".format(int(n_orbit), int(dt), int(dtrt))

    input_dir = os.path.join(datapath, orbit_dt_str)

    pickle_dir_to_csvs(input_dir, filename_keyword="short")

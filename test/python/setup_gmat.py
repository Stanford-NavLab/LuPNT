import sys
import os

gmat_path = os.environ["GMAT_PATH"]
assert gmat_path, "Please set the GMAT_PATH environment variable"

startup = os.path.join(gmat_path, "bin", "api_startup_file.txt")
assert os.path.exists(startup), "Cannot find " + startup

sys.path.insert(1, os.path.join(gmat_path, "bin"))
import gmatpy as gmat

gmat.Setup(startup)

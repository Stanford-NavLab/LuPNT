import os
import requests
import zipfile
import shutil

DATA_URL = "https://bit.ly/LuPNT_data"
DATA_FILENAME = "LuPNT_data.zip"
DATA_FOLDERNAME = "LuPNT_data"


if "LUPNT_DATA_PATH" not in os.environ or not os.path.isdir(
    os.path.join(os.environ["LUPNT_DATA_PATH"], "ephemeris")
):
    print("Downloading required data from", DATA_URL)
    response = requests.get(DATA_URL, stream=True)
    with open(DATA_FILENAME, "wb") as f:
        shutil.copyfileobj(response.raw, f)
    with zipfile.ZipFile(DATA_FILENAME, "r") as zip_ref:
        zip_ref.extractall()
    os.remove(DATA_FILENAME)
    os.environ["LUPNT_DATA_PATH"] = os.path.join(os.getcwd(), DATA_FOLDERNAME)
    print("Downloaded data to", os.path.relpath(os.environ["LUPNT_DATA_PATH"]))
else:
    print("Found required data at", os.path.relpath(os.environ["LUPNT_DATA_PATH"]))

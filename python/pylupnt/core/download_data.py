import os
import requests
import zipfile
import shutil

# The LuPNT_data reference bundle (ephemeris, GNSS antenna/clock products, plasma
# coefficients, gravity fields, DEMs, TLEs) is hosted as a GitHub Release asset on
# the public repo. Downloads of release assets on public repos are unauthenticated
# and unmetered, so no token or manual setup is required.
#
# To publish a NEW bundle, create a release and attach LuPNT_data.zip (whose root is
# a single LuPNT_data/ folder), then bump DATA_TAG to match:
#     cd data && zip -r ../LuPNT_data.zip LuPNT_data && cd ..
#     gh release create data-YYYY-MM-DD LuPNT_data.zip \
#         --repo Stanford-NavLab/LuPNT --title "LuPNT_data YYYY-MM-DD" \
#         --notes "LuPNT reference-data bundle"
DATA_TAG = "data-2026-07-30"
DATA_URL = f"https://github.com/Stanford-NavLab/LuPNT/releases/download/{DATA_TAG}/LuPNT_data.zip"
DATA_FILENAME = "LuPNT_data.zip"
DATA_FOLDERNAME = "LuPNT_data"


if "LUPNT_DATA_PATH" not in os.environ or not os.path.isdir(
    os.path.join(os.environ["LUPNT_DATA_PATH"], "ephemeris")
):
    print("Fetching data from", DATA_URL)
    response = requests.get(DATA_URL, stream=True)
    # Fail loudly on a missing/renamed release asset instead of writing the 404 page
    # into the zip and failing later with an opaque "not a zip file" error.
    response.raise_for_status()
    with open(DATA_FILENAME, "wb") as f:
        shutil.copyfileobj(response.raw, f)
    with zipfile.ZipFile(DATA_FILENAME, "r") as zip_ref:
        zip_ref.extractall()
    os.remove(DATA_FILENAME)
    os.environ["LUPNT_DATA_PATH"] = os.path.join(os.getcwd(), DATA_FOLDERNAME)
    print("Downloaded data to", os.path.relpath(os.environ["LUPNT_DATA_PATH"]))

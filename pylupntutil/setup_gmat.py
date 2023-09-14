import sys
from os import path

# username = "keidai"
username = "Guillem"

apistartup = "api_startup_file.txt"

if username == "keidai":
    GmatInstall = "/Users/keidaiiiyama/Documents/GMAT_R2022a"
elif username == "Guillem":
    GmatInstall = "/Users/guillemcv/Applications/GMAT R2022a"

GmatBinPath = GmatInstall + "/bin"
Startup = GmatBinPath + "/" + apistartup

if path.exists(Startup):
    sys.path.insert(1, GmatBinPath)

    import gmatpy as gmat

    gmat.Setup(Startup)

else:
    print("Cannot find ", Startup)
    print()
    print(
        "Please set up a GMAT startup file named ",
        apistartup,
        " in the ",
        GmatBinPath,
        " folder.",
    )

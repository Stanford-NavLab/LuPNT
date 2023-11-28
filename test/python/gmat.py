import sys
import os
import numpy as np

gmat_path = os.environ["GMAT_PATH"]
assert gmat_path, "Please set the GMAT_PATH environment variable"

startup = os.path.join(gmat_path, "bin", "api_startup_file.txt")
assert os.path.exists(startup), "Cannot find " + startup

sys.path.insert(1, os.path.join(gmat_path, "bin"))
import gmatpy as gmat

gmat.Setup(startup)


def create_spacecraft(coe, *args, **kwargs):
    # Spacecraft configuration preliminaries
    earthorb = gmat.Construct("Spacecraft", "EarthOrbiter")
    earthorb.SetField("DateFormat", "UTCGregorian")
    earthorb.SetField("Epoch", "20 Jul 2020 12:00:00.000")

    earthorb.SetField("CoordinateSystem", "EarthMJ2000Eq")
    earthorb.SetField("DisplayStateType", "Keplerian")

    # Orbital state
    earthorb.SetField("SMA", coe[0])
    earthorb.SetField("ECC", coe[1])
    earthorb.SetField("INC", np.rad2deg(coe[2]))
    earthorb.SetField("RAAN", np.rad2deg(coe[3]))
    earthorb.SetField("AOP", np.rad2deg(coe[4]))
    earthorb.SetField("MA", np.rad2deg(coe[5]))

    # Spacecraft ballistic properties for the SRP and Drag models
    if "SRPArea" in kwargs:
        earthorb.SetField("SRPArea", 2.5)
    if "Cr" in kwargs:
        earthorb.SetField("Cr", 1.75)
    if "DragArea" in kwargs:
        earthorb.SetField("DragArea", 1.8)
    if "Cd" in kwargs:
        earthorb.SetField("Cd", 2.1)

    earthorb.SetField("DryMass", 80)


def create_force_model():
    # Force model settings
    fm = gmat.Construct("ForceModel", "FM")
    fm.SetField("CentralBody", "Earth")

    # An 8x8 JGM-3 Gravity Model
    earthgrav = gmat.Construct("GravityField")
    earthgrav.SetField("BodyName", "Earth")
    earthgrav.SetField("PotentialFile", "../data/gravity/earth/JGM3.cof")
    earthgrav.SetField("Degree", 8)
    earthgrav.SetField("Order", 8)

    # Add forces into the ODEModel container
    fm.AddForce(earthgrav)

    gmat.Initialize()

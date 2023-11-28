import pylupnt as pnt
import numpy as np
import utils
import os
from gmat import gmat, gmat_path

kwargs = {}
coe = utils.coe_array_elfo
earth_potential_file_path = os.path.join(
    gmat_path, "data", "gravity", "earth", "JGM3.cof"
)
luna_potential_file_path = os.path.join(
    gmat_path, "data", "gravity", "luna", "grgm900c.cof"
)

EarthMJ2000Eq = gmat.Construct("CoordinateSystem", "EarthMJ2000Eq", "Earth", "MJ2000Eq")
LunaMJ2000Eq = gmat.Construct("CoordinateSystem", "LunaMJ2000Eq", "Luna", "MJ2000Ec")

# Spacecraft configuration preliminaries
sc = gmat.Construct("Spacecraft", "LunaOrbiter")
sc.SetField("DateFormat", "UTCGregorian")
sc.SetField("Epoch", "20 Jul 2020 12:00:00.000")
sc.SetField("CoordinateSystem", "LunaMJ2000Eq")
# sc.SetField("DisplayStateType", "Keplerian")

# Orbital state
sc.SetField("SMA", coe[0])
sc.SetField("ECC", coe[1])
sc.SetField("INC", np.rad2deg(coe[2]))
sc.SetField("RAAN", np.rad2deg(coe[3]))
sc.SetField("AOP", np.rad2deg(coe[4]))
sc.SetField("TA", np.rad2deg(pnt.mean_to_true_anomaly(coe[5], coe[1])))

# Spacecraft ballistic properties for the SRP and Drag models
if "SRPArea" in kwargs:
    sc.SetField("SRPArea", 2.5)
if "Cr" in kwargs:
    sc.SetField("Cr", 1.75)
if "DragArea" in kwargs:
    sc.SetField("DragArea", 1.8)
if "Cd" in kwargs:
    sc.SetField("Cd", 2.1)

sc.SetField("DryMass", 80)

# Force model settings
fm = gmat.Construct("ForceModel", "FM")
fm.SetField("CentralBody", "Luna")

# An 8x8 JGM-3 Gravity Model
grav = gmat.Construct("GravityField")
grav.SetField("BodyName", "Luna")
grav.SetField("PotentialFile", luna_potential_file_path)
grav.SetField("Degree", 0)
grav.SetField("Order", 0)

# Add forces into the ODEModel container
fm.AddForce(grav)

gmat.Initialize()

# Build the propagation container class
pdprop = gmat.Construct("Propagator", "PDProp")

# Create and assign a numerical integrator for use in the propagation
gator = gmat.Construct("PrinceDormand78", "Gator")
pdprop.SetReference(gator)

# Set some of the fields for the integration
pdprop.SetField("InitialStepSize", 60.0)
pdprop.SetField("Accuracy", 1.0e-12)
pdprop.SetField("MinStep", 0.0)

# Perform top level initialization
gmat.Initialize()

# Setup the spacecraft that is propagated
pdprop.AddPropObject(sc)
pdprop.PrepareInternals()

# Refresh the 'gator reference
fm = pdprop.GetODEModel()
gator = pdprop.GetPropagator()

# Take a 60 second step, showing the state before and after propagation
print(sc.GetEpoch())
print(gator.GetState())
print(sc.GetState().GetState())
print(sc.GetCartesianState())
gator.Step(30 * 60)
print(sc.GetEpoch())
print(gator.GetState())
print(sc.GetState().GetState())
print(sc.GetCartesianState())
gator.Step(30 * 60)
print(sc.GetEpoch())
print(gator.GetState())
print(sc.GetState().GetState())
print(sc.GetCartesianState())

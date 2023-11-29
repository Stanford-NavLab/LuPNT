import numpy as np
import pylupnt as pnt

try:
    from .gmat import gmat
except ImportError:
    from gmat import gmat


def unpack_rvector(rvector):
    return np.array([rvector.Get(i) for i in range(6)])


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


# pylupnt.CoordSystem {
#   ITRF = 0,  // International Terrestrial Reference Frame
#   GCRF,      // Geocentric Reference System
#   ICRF,      // International Celestial Reference System
#   SER,       // Sun-Earth Rotating Frame
#   GSE,       // Geocentric Solar Ecliptic
#   EME,       // Earth-Centered mean equator and equinox at J2000 epoch
#   MOD,       // Mean of date equatorial system
#   TOD,       // True of date equatorial system
#   EMR,       // Earth-Moon Rotating Frame
#   MI,        // Moon-centered Inertial Frame  (Axis aligened with ICRF)
#   PA,        // Moon-Fixed with principal axes
#   ME,        // Moon-Fixed with mean-Earth / polar axes
#   RTN,       // Radial-Tangential-Normal
#   CoordSystemCount,
#   NONE,
# };

# Gmat axes
# creatables = [
#     "MJ2000Eq",
#     "MJ2000Ec",
#     "TOEEq",
#     "TOEEc",
#     "MOEEq",
#     "MOEEc",
#     "TODEq",
#     "TODEc",
#     "MODEq",
#     "MODEc",
#     "ObjectReferenced",
#     "Equator",
#     "BodyFixed",
#     "BodyInertial",
#     "GSE",
#     "GSM",
#     "Topocentric",
#     "LocalAlignedConstrained",
#     #"ITRF", # This one is commented out in the original code
#     "ICRF",
#     "BodySpinSun"
# ]


def get_coordinate_system(name):
    if name == pnt.CoordSystem.ITRF:
        return gmat.Construct("CoordinateSystem", "ITRF", "Earth", "BodyFixed")
    elif name == pnt.CoordSystem.GCRF:
        return gmat.Construct("CoordinateSystem", "GCRF", "Earth", "MJ2000Eq")
    elif name == pnt.CoordSystem.ICRF:
        return gmat.Construct("CoordinateSystem", "ICRF", "Earth", "ICRF")
    elif name == pnt.CoordSystem.SER:
        assert False, "Not implemented"
    elif name == pnt.CoordSystem.GSE:
        return gmat.Construct("CoordinateSystem", "GSE", "Earth", "GSE")
    elif name == pnt.CoordSystem.EME:
        assert False, "Not implemented"
    elif name == pnt.CoordSystem.MOD:
        return gmat.Construct("CoordinateSystem", "MOD", "Earth", "MODEq")
    elif name == pnt.CoordSystem.TOD:
        return gmat.Construct("CoordinateSystem", "TOD", "Earth", "TODEq")
    elif name == pnt.CoordSystem.MI:
        return gmat.Construct("CoordinateSystem", "MI", "Moon", "MJ2000Eq")
    elif name == pnt.CoordSystem.PA:
        return gmat.Construct("CoordinateSystem", "PA", "Moon", "BodyFixed")
    elif name == pnt.CoordSystem.ME:
        assert False, "Not implemented"
    else:
        assert False, "Coordinate system not found"

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


# pylupnt.CoordSystem = {ITRF, GCRF, ICRF, SER, GSE, EME, MOD, TOD, MI, PA, ME, NONE}
# gmat_axes = [
# "MJ2000Eq", "MJ2000Ec", "TOEEq", "TOEEc", "MOEEq", "MOEEc", "TODEq", "TODEc", "MODEq", "MODEc",
# "ObjectReferenced", "Equator", "BodyFixed", "BodyInertial", "GSE", "GSM", "Topocentric",
# "LocalAlignedConstrained", "ICRF", "BodySpinSun"]


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
    elif name == pnt.CoordSystem.EMR:
        # GMAT Code Broken. Not working
        coord_sys = gmat.Construct("CoordinateSystem", "EMR")
        coord_sys.SetField("Origin", "Earth")
        coord_sys.SetField("Axes", "ObjectReferenced")
        coord_sys.SetField("XAxis", "R")
        coord_sys.SetField("ZAxis", "N")
        coord_sys.SetField("Primary", "Earth")
        coord_sys.SetField("Secondary", "Luna")
    elif name == pnt.CoordSystem.MI:
        return gmat.Construct("CoordinateSystem", "MI", "Luna", "MJ2000Eq")
    elif name == pnt.CoordSystem.PA:
        return gmat.Construct("CoordinateSystem", "PA", "Luna", "BodyFixed")
    elif name == pnt.CoordSystem.ME:
        assert False, "Not implemented"
    else:
        assert False, "Coordinate system not found"


def convert_coord(epoch, rv, coord_sys_from, coord_sys_to):
    epoch_gmat = gmat.GmatTime(convert_pylupnt_to_gmat_epoch(epoch))
    coord_from = get_coordinate_system(coord_sys_from)
    coord_to = get_coordinate_system(coord_sys_to)
    state_from = gmat.Rvector6(rv)
    state_to = gmat.Rvector6()
    gmat.Initialize()

    converter = gmat.CoordinateConverter()
    converter.Convert(epoch_gmat, state_from, coord_from, state_to, coord_to)
    return np.array([state_to[i] for i in range(state_to.GetSize())])


def convert_gmat_to_pylupnt_epoch(epoch):
    # pylupnt: TAI seconds since 2000/01/01 12:00:00.000
    # gmat: A1MJD since 1941/
    return epoch * gmat.SECS_PER_DAY - gmat.MJD_OF_J2000 * gmat.SECS_PER_DAY


def convert_pylupnt_to_gmat_epoch(epoch):
    # pylupnt: TAI seconds since J2000
    # gmat: A1MJD since 1941
    return (epoch + gmat.MJD_OF_J2000 * gmat.SECS_PER_DAY) / gmat.SECS_PER_DAY


def convert_time(epoch, time_sys_from, time_sys_to):
    time_converter = gmat.TimeSystemConverter.Instance()
    epoch_gmat_from = convert_pylupnt_to_gmat_epoch(epoch)
    epoh_gmat_to = time_converter.Convert(epoch_gmat_from, time_sys_from, time_sys_to)
    epoch_to = convert_gmat_to_pylupnt_epoch(epoh_gmat_to)
    return epoch_to

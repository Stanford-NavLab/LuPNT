import numpy as np
import pylupnt as pnt


def get_gnss_str(gnss_const):
    if gnss_const == "GPS":
        return "G"
    elif gnss_const == "GALILEO":
        return "E"
    elif gnss_const == "QZSS":
        return "J"
    else:
        raise ValueError(f"Unknown GNSS constellation: {gnss_const}")


def get_freq_str(gnss_const, signal_family):
    if gnss_const in ["GPS", "QZSS"]:
        if signal_family == 1:
            return "L1"
        elif signal_family == 2:
            return "L2"
        elif signal_family == 5:
            return "L5"
    elif gnss_const == "GALILEO":
        if signal_family == 1:
            return "E1"
        elif signal_family == 5:
            return "E5a"
    raise ValueError(f"Unknown frequency for {gnss_const} signal family {signal_family}")


def unit(v, eps=1e-12):
    norm = np.linalg.norm(v)
    if norm < eps:
        raise ValueError("Cannot normalize zero-length vector.")
    return v / norm


def enu_to_ecef_rot(t_tai, r_sat_ecef):
    r = np.asarray(r_sat_ecef, dtype=float).reshape(3)

    # Up: radial outward (geocentric up)
    u = unit(r, eps=1e-12)

    # Earth spin axis in ECEF (approx): z-axis
    z = np.array([0.0, 0.0, 1.0])

    # East: tangent direction
    # e = z x u (points east unless near poles)
    e_raw = np.cross(z, u)
    if np.linalg.norm(e_raw) < 1e-12:
        # Satellite nearly over pole; choose a different reference axis
        # Use x-axis as fallback
        x = np.array([1.0, 0.0, 0.0])
        e_raw = np.cross(x, u)
        if np.linalg.norm(e_raw) < 1e-12:
            raise ValueError("Cannot define East vector (degenerate geometry).")
    e = unit(e_raw, eps=1e-12)

    # North: complete right-handed frame (N = U x E)
    n = unit(np.cross(u, e), eps=1e-12)
    Cneu = np.column_stack((n, e, u))  # columns are NEU axes in ECEF

    return Cneu


def ijk_to_ecef_rot(t_tai, r_sat_ecef):
    r = np.asarray(r_sat_ecef, dtype=float).reshape(3)

    # z-axis: radial outward (geocentric up)
    kvec = -unit(r, eps=1e-12)

    # x-axis: projection of ECEF x-axis onto plane perpendicular to z
    r_sun = pnt.get_body_pos_vel(t_tai, pnt.EARTH, pnt.SUN, pnt.ECEF)[:3]
    jvec = unit(r_sun - r_sat_ecef, eps=1e-12)  # approximate Earth-Sun direction

    ivec = np.cross(jvec, kvec)

    Cijk = np.column_stack((ivec, jvec, kvec))  # columns are IJK axes in ECEF

    return Cijk

# ---- spice_adapter.py ----
import numpy as np
from datetime import datetime, timezone
import pylupnt as pnt


def moon_state_spice(t_dt: datetime):
    """
    Returns (rM_eci[km], vM_eci[km/s]) of the Moon wrt Earth center in J2000.
    """
    t_tai = pnt.datetime_to_tai(t_dt)
    # 'MOON' state relative to 'EARTH' in 'J2000'; no aberration correction for dynamics
    # spkezr returns (state, light_time); state in km and km/s
    state = pnt.get_body_pos_vel(t_tai, pnt.EARTH, pnt.MOON, pnt.ECI)
    r = np.asarray(state[:3], dtype=float)
    v = np.asarray(state[3:], dtype=float)
    return r, v


def sun_state_spice(t_dt: datetime):
    """
    Returns (rE_eci[km], vE_eci[km/s]) of the Earth center in J2000.
    """
    t_tai = pnt.datetime_to_tai(t_dt)
    # 'EARTH' state relative to 'EARTH' in 'J2000'; no aberration correction for dynamics
    state = pnt.get_body_pos_vel(t_tai, pnt.EARTH, pnt.SUN, pnt.ECI)
    r = np.asarray(state[:3], dtype=float)
    v = np.asarray(state[3:], dtype=float)
    return r, v

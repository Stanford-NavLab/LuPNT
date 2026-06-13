#!/usr/bin/env python3
"""
GNSSAttitude: build a NEU<->ECEF rotation at a satellite position.

NEU is defined locally at the satellite:
  U: outward radial direction (from Earth's center to satellite)
  E: East (tangent, increasing longitude direction)
  N: North (tangent, increasing latitude direction)

This is a *local topocentric* frame attached to the satellite position.
It is not the satellite body frame (yaw-steering).

Author: you + ChatGPT
"""

from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from typing import Callable, Optional


def _unit(v: np.ndarray, eps: float = 1e-15) -> np.ndarray:
    v = np.asarray(v, dtype=float).reshape(3)
    n = np.linalg.norm(v)
    if n < eps:
        raise ValueError("Cannot normalize near-zero vector.")
    return v / n


@dataclass
class GNSSAttitude:
    """
    Provide rotation matrices between local NEU at the satellite and ECEF.

    Parameters
    ----------
    get_sun_vector : callable
        Function get_sun_vector(t_tai) -> (3,) ndarray in ECEF.
        Not required for NEU frame, but accepted so you can extend this class
        to yaw-steering later without changing API.
    """

    get_sun_vector: Optional[Callable[[float], np.ndarray]] = None
    eps: float = 1e-12

    def C_ecef_from_neu(self, t_tai: float, r_sat_ecef: np.ndarray) -> np.ndarray:
        """
        Return C_ecef_neu (3x3) such that:
            v_ecef = C_ecef_neu @ v_neu

        NEU basis vectors (as columns) expressed in ECEF:
            C = [n_ecef, e_ecef, u_ecef]
        """
        r = np.asarray(r_sat_ecef, dtype=float).reshape(3)

        # Up: radial outward (geocentric up)
        u = _unit(r, eps=self.eps)

        # Earth spin axis in ECEF (approx): z-axis
        z = np.array([0.0, 0.0, 1.0])

        # East: tangent direction
        # e = z x u (points east unless near poles)
        e_raw = np.cross(z, u)
        if np.linalg.norm(e_raw) < self.eps:
            # Satellite nearly over pole; choose a different reference axis
            # Use x-axis as fallback
            x = np.array([1.0, 0.0, 0.0])
            e_raw = np.cross(x, u)
            if np.linalg.norm(e_raw) < self.eps:
                raise ValueError("Cannot define East vector (degenerate geometry).")
        e = _unit(e_raw, eps=self.eps)

        # North: complete right-handed frame (N = U x E)
        n = _unit(np.cross(u, e), eps=self.eps)

        C = np.column_stack((n, e, u))  # columns are NEU axes in ECEF
        return C

    def C_neu_from_ecef(self, t_tai: float, r_sat_ecef: np.ndarray) -> np.ndarray:
        """
        Return C_neu_ecef (3x3) such that:
            v_neu = C_neu_ecef @ v_ecef

        This is the transpose of C_ecef_from_neu because it's orthonormal.
        """
        C = self.C_ecef_from_neu(t_tai, r_sat_ecef)
        return C.T

    def sun_unit_neu(self, t_tai: float, r_sat_ecef: np.ndarray) -> np.ndarray:
        """
        OPTIONAL helper: return Sun direction unit vector expressed in NEU
        (useful if you later implement yaw-steered body axes).
        Requires get_sun_vector.
        """
        if self.get_sun_vector is None:
            raise ValueError("get_sun_vector(t_tai) was not provided.")
        s_ecef = np.asarray(self.get_sun_vector(t_tai), dtype=float).reshape(3)
        s_hat_ecef = _unit(s_ecef, eps=self.eps)
        C_neu_ecef = self.C_neu_from_ecef(t_tai, r_sat_ecef)
        return C_neu_ecef @ s_hat_ecef


# ---------------- example usage ----------------
if __name__ == "__main__":
    # Dummy Sun vector function (ECEF), replace with your real one
    def get_sun_vector(t_tai: float) -> np.ndarray:
        # Example: fixed direction (NOT real)
        return np.array([1.0, 0.2, 0.1])

    att = GNSSAttitude(get_sun_vector=get_sun_vector)

    # Example satellite position in ECEF (m)
    r_sat = np.array([15600e3, 0.0, 20180e3])
    t = 0.0  # your TAI seconds, MJD, etc.

    C_ecef_neu = att.C_ecef_from_neu(t, r_sat)
    print("C_ecef_from_neu:\n", C_ecef_neu)

    # Example: convert a NEU vector to ECEF
    v_neu = np.array([1.0, 2.0, 3.0])
    v_ecef = C_ecef_neu @ v_neu
    print("v_ecef:", v_ecef)

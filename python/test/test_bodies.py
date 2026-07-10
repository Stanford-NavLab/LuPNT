"""Unit tests for `pylupnt` body definitions and physical constants:
the `Body` factory helpers, `get_physical_constants`, and the module-level
gravitational-parameter / radius / misc constants.

All values are SI (metres, seconds, m^3/s^2) per the top-level pylupnt API.
"""

import numpy as np
import pylupnt as pnt
import pytest


# ---------------------------------------------------------------------------
# Body factory helpers
# ---------------------------------------------------------------------------


def test_body_gm_and_radius():
    moon = pnt.Body.Moon()
    earth = pnt.Body.Earth()
    sun = pnt.Body.Sun()

    assert moon.name == "MOON"
    # GM [m^3/s^2]
    assert moon.GM == pytest.approx(4.9028e12, rel=1e-3)
    assert earth.GM == pytest.approx(3.986004e14, rel=1e-3)
    assert sun.GM == pytest.approx(1.327e20, rel=1e-3)
    # Mean radius [m]
    assert moon.R == pytest.approx(1.7374e6, rel=1e-3)
    assert earth.R == pytest.approx(6.378e6, rel=1e-2)

    # Sun is far more massive than Earth, which is more massive than the Moon.
    assert sun.GM > earth.GM > moon.GM


def test_body_gm_matches_module_constants():
    assert pnt.Body.Moon().GM == pytest.approx(pnt.GM_MOON, rel=1e-6)
    assert pnt.Body.Earth().GM == pytest.approx(pnt.GM_EARTH, rel=1e-6)


# ---------------------------------------------------------------------------
# Module-level physical constants (SI)
# ---------------------------------------------------------------------------


def test_speed_of_light_si():
    assert pnt.C == pytest.approx(299_792_458.0, rel=1e-9)


def test_gravitational_parameters_si():
    # GM in m^3/s^2 -> these magnitudes confirm SI (not km-based) units.
    assert pnt.GM_EARTH == pytest.approx(3.986004e14, rel=1e-3)
    assert pnt.GM_MOON == pytest.approx(4.9028e12, rel=1e-3)


def test_body_radii_si():
    assert pnt.R_EARTH == pytest.approx(6.378e6, rel=1e-2)
    assert pnt.R_MOON == pytest.approx(1.7374e6, rel=1e-3)


def test_angle_constants():
    assert pnt.RAD == pytest.approx(np.pi / 180.0, rel=1e-12)
    assert pnt.DEG == pytest.approx(180.0 / np.pi, rel=1e-12)
    assert 180.0 * pnt.RAD == pytest.approx(np.pi, rel=1e-12)


def test_time_unit_constants():
    assert pnt.SECS_DAY == pytest.approx(86400.0, rel=1e-12)
    assert pnt.SECS_HOUR == pytest.approx(3600.0, rel=1e-12)


def test_get_physical_constants_returns_data():
    pc = pnt.get_physical_constants()
    assert pc is not None

"""Unit tests for the `pylupnt` time-conversion bindings
(`gregorian_to_time`, `convert_time`, `time_to_mjd`/`mjd_to_time`,
`mjd_to_gregorian`, `gregorian_to_mjd`, `time_to_gregorian_string`).

These exercise the C++ time system through its Python bindings: known
time-scale offsets, round-trip invertibility (both across time scales and
between the Gregorian/MJD/internal representations), and vectorized calls.
"""

import numpy as np
import pylupnt as pnt
import pytest

# Internal LuPNT epochs are seconds since J2000; conversions are done in double
# precision, so a few hundred ns of round-trip slack is expected.
ABS_TOL_S = 1e-5

ALL_SCALES = [
    pnt.Time.TAI,
    pnt.Time.UTC,
    pnt.Time.TT,
    pnt.Time.TDB,
    pnt.Time.GPS,
]


def _tai(y, mo, d, h=0, mi=0, s=0.0):
    """A TAI epoch for a UTC-ish calendar date (via the internal ~TDB scale)."""
    t = pnt.gregorian_to_time(y, mo, d, h, mi, s)
    return pnt.convert_time(t, pnt.Time.TDB, pnt.Time.TAI)


# ---------------------------------------------------------------------------
# Known constant time-scale offsets
# ---------------------------------------------------------------------------


def test_tt_tai_offset_is_exact():
    """TT = TAI + 32.184 s exactly (definitional)."""
    tai = _tai(2026, 1, 14)
    tt = pnt.convert_time(tai, pnt.Time.TAI, pnt.Time.TT)
    assert float(tt - tai) == pytest.approx(32.184, abs=ABS_TOL_S)


def test_gps_tai_offset_is_exact():
    """GPS = TAI - 19 s exactly (definitional, no leap seconds since 1980)."""
    tai = _tai(2026, 1, 14)
    gps = pnt.convert_time(tai, pnt.Time.TAI, pnt.Time.GPS)
    assert float(gps - tai) == pytest.approx(-19.0, abs=ABS_TOL_S)


def test_tai_utc_offset_is_integer_leap_seconds():
    """TAI - UTC is a positive integer number of leap seconds (37 s in 2026)."""
    tai = _tai(2026, 1, 14)
    utc = pnt.convert_time(tai, pnt.Time.TAI, pnt.Time.UTC)
    offset = float(tai - utc)
    assert offset == pytest.approx(round(offset), abs=1e-6)
    assert offset >= 32.0  # at least the 2006-era count; 37 as of 2017+


# ---------------------------------------------------------------------------
# Round trips across time scales
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("to_scale", ALL_SCALES)
def test_convert_time_round_trip(to_scale):
    """convert_time(t, TAI, X) then X -> TAI recovers the original epoch."""
    tai = _tai(2026, 3, 21, 6, 30, 12.25)
    other = pnt.convert_time(tai, pnt.Time.TAI, to_scale)
    back = pnt.convert_time(other, to_scale, pnt.Time.TAI)
    assert float(back - tai) == pytest.approx(0.0, abs=ABS_TOL_S)


def test_convert_time_identity():
    """Converting within the same scale is a no-op."""
    tai = _tai(2026, 1, 14)
    same = pnt.convert_time(tai, pnt.Time.TAI, pnt.Time.TAI)
    assert float(same - tai) == pytest.approx(0.0, abs=1e-9)


def test_convert_time_vectorized():
    """convert_time is elementwise over a numpy array and matches scalar calls."""
    t0 = pnt.gregorian_to_time(2026, 1, 14, 0, 0, 0)
    arr = t0 + np.arange(5) * 3600.0
    tai_arr = np.asarray(pnt.convert_time(arr, pnt.Time.TDB, pnt.Time.TAI))
    assert tai_arr.shape == (5,)
    for k in range(5):
        scalar = float(pnt.convert_time(arr[k], pnt.Time.TDB, pnt.Time.TAI))
        assert tai_arr[k] == pytest.approx(scalar, abs=ABS_TOL_S)


# ---------------------------------------------------------------------------
# Gregorian / MJD / internal representation round trips
# ---------------------------------------------------------------------------


def test_gregorian_mjd_round_trip():
    """gregorian -> time -> MJD -> gregorian recovers the calendar fields."""
    fields = (2026, 1, 14, 12, 30, 15.5)
    t = pnt.gregorian_to_time(*fields)
    mjd = pnt.time_to_mjd(t)
    y, mo, d, h, mi, s = pnt.mjd_to_gregorian(mjd)
    assert (y, mo, d, h, mi) == fields[:5]
    assert s == pytest.approx(fields[5], abs=1e-3)


def test_mjd_time_round_trip():
    """time -> MJD -> time is the identity."""
    t = pnt.gregorian_to_time(2026, 7, 9, 3, 14, 15.9)
    t2 = pnt.mjd_to_time(pnt.time_to_mjd(t))
    assert float(t2 - t) == pytest.approx(0.0, abs=ABS_TOL_S)


def test_gregorian_to_mjd_matches_time_to_mjd():
    """gregorian_to_mjd is consistent with time_to_mjd(gregorian_to_time(...))."""
    fields = (2026, 1, 14, 6, 0, 0.0)
    mjd_direct = pnt.gregorian_to_mjd(*fields)
    mjd_via_time = pnt.time_to_mjd(pnt.gregorian_to_time(*fields))
    assert float(mjd_direct) == pytest.approx(float(mjd_via_time), abs=1e-9)


def test_time_to_gregorian_string_is_wellformed():
    """The formatted timestamp contains the calendar date it was built from."""
    t = pnt.gregorian_to_time(2026, 1, 14, 0, 0, 0)
    s = pnt.time_to_gregorian_string(t)
    assert "2026" in s
    assert isinstance(s, str)


def test_time_ordering_is_monotonic():
    """Later calendar epochs map to larger internal time values."""
    t_earlier = pnt.gregorian_to_time(2026, 1, 14, 0, 0, 0)
    t_later = pnt.gregorian_to_time(2026, 1, 14, 0, 0, 1)
    assert float(t_later) > float(t_earlier)
    assert float(t_later - t_earlier) == pytest.approx(1.0, abs=1e-6)

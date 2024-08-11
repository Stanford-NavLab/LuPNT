import pylupnt as pnt
import numpy as np
import pytest


def test_constants():
    assert pnt.GM_EARTH == pytest.approx(398600.435507)

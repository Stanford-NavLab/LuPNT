import pylupnt as pnt
import numpy as np
import pytest


def test_frame_converter():
    rv_in = np.array([1000, 2000, 3000, 0.1, 0.2, 0.3])
    frame_in = pnt.Frame.MOON_CI
    frame_out = pnt.Frame.MOON_PA

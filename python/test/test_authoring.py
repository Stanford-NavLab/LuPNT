"""Tests for authoring simulation components in pure Python: subclassing ``pnt.Measurement`` and
``pnt.Application`` (pybind11 trampolines) and ``pnt.register_application``."""

import numpy as np
import pytest

import pylupnt as pnt


def test_measurement_subclass_evaluate_dispatches_to_python():
    class Bearing(pnt.Measurement):
        def __init__(self, sigma):
            pnt.Measurement.__init__(self)
            self.sigma = sigma
            self.target = np.zeros(3)

        def compute(self, x):
            d = self.target - x[:3]
            rn = np.linalg.norm(d)
            u = d / rn
            H = np.zeros((3, 6))
            H[:, :3] = -(np.eye(3) - np.outer(u, u)) / rn
            return u, H, self.sigma**2 * np.eye(3)

    m = Bearing(1e-5)
    m.target = np.array([1.0e6, 5.0e5, 3.0e5])
    x = np.array([2.0e5, 0.0, 0.0, 0.0, 1500.0, 300.0])
    z, H, R = m.evaluate(x)  # runs through the C++ Measurement base -> Python compute
    assert z.shape == (3,)
    np.testing.assert_allclose(np.linalg.norm(z), 1.0, atol=1e-9)  # unit line-of-sight
    assert H.shape == (3, 6)
    np.testing.assert_allclose(R, (1e-5) ** 2 * np.eye(3))


def test_measurement_without_compute_raises():
    class Bad(pnt.Measurement):
        pass

    with pytest.raises(Exception):
        Bad().evaluate(np.zeros(6))


def test_application_subclass_name_and_frequency():
    class MyApp(pnt.Application):
        def __init__(self):
            pnt.Application.__init__(self)

        def step(self, t):
            pass

    a = MyApp()
    a.set_name("my_app")
    a.set_frequency(2.0)
    assert a.get_name() == "my_app"
    assert a.get_frequency() == pytest.approx(2.0)


def test_register_application_accepts_python_subclass():
    class RegApp(pnt.Application):
        def __init__(self, config):
            pnt.Application.__init__(self)

        def step(self, t):
            pass

    # Registering a Python Application subclass into the C++ asset factory must not raise.
    pnt.register_application("TestRegApp_unit", RegApp)

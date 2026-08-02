"""Tests for pylupnt.plot._mpl (matplotlib plotting helpers)."""

import matplotlib

matplotlib.use("Agg")  # headless backend BEFORE pyplot / pylupnt.plot import

import matplotlib.pyplot as plt
import numpy as np
import pytest

import pylupnt as pnt
from pylupnt.plot import _mpl as M


@pytest.fixture(autouse=True)
def _close_figs():
    yield
    plt.close("all")


@pytest.fixture(scope="module")
def antenna():
    return pnt.Antenna("Block-IIF_ACE")


# --------------------------------------------------------------------------- #
# Plot3D                                                                       #
# --------------------------------------------------------------------------- #
def test_plot3d_init_defaults():
    p = M.Plot3D()
    assert isinstance(p.fig, plt.Figure)
    assert p.scatters == []
    assert p.plots == []


def test_plot3d_init_custom_view():
    p = M.Plot3D(azim=10, elev=45, figsize=(4, 4))
    assert p.fig.get_size_inches()[0] == pytest.approx(4)


def test_plot_surface_moon_adjusts_axis():
    p = M.Plot3D()
    p.plot_surface(pnt.MOON, scale=20)
    assert p.name == pnt.MOON


def test_plot_surface_earth_no_adjust_with_limit_offset():
    p = M.Plot3D()
    p.plot_surface(
        pnt.EARTH,
        offset=np.array([100.0, 0.0, 0.0]),
        adjust_axis=False,
        limit=5000.0,
        scale=30,
    )
    assert p.name == pnt.EARTH


def test_plot_surface_earth_with_limit_adjust():
    p = M.Plot3D()
    p.plot_surface(pnt.EARTH, adjust_axis=True, limit=8000.0, scale=30)
    assert p.name == pnt.EARTH


# --------------------------------------------------------------------------- #
# occultation / visibility                                                     #
# --------------------------------------------------------------------------- #
def _surface_plot(body=pnt.MOON):
    p = M.Plot3D()
    p.plot_surface(body, scale=20)
    return p


def test_check_occultation():
    p = _surface_plot()
    data = np.random.default_rng(0).standard_normal((6, 3)) * 3000
    cond, alphas = p.check_occultation(data)
    assert cond.shape == (6,)
    assert alphas.shape == (6,)


def test_plot_visible():
    p = _surface_plot()
    (line,) = p.ax.plot([0.0], [0.0], [0.0])
    p.points = [line]
    p.data = [np.random.default_rng(1).standard_normal((5, 3)) * 3000]
    p.plot_visible(p.azim, p.elev)


def test_rotate_ignores_other_axes():
    p = _surface_plot()

    class _Evt:
        inaxes = None

    # inaxes != p.ax -> no-op branch
    p.rotate(_Evt())


def test_rotate_matching_axis_calls_plot_visible():
    p = _surface_plot()
    (line,) = p.ax.plot([0.0], [0.0], [0.0])
    p.points = [line]
    p.data = [np.random.default_rng(8).standard_normal((5, 3)) * 3000]

    class _Evt:
        pass

    evt = _Evt()
    evt.inaxes = p.ax  # matches -> triggers plot_visible branch
    p.rotate(evt)


# --------------------------------------------------------------------------- #
# scatter / plot                                                               #
# --------------------------------------------------------------------------- #
def test_scatter_no_mask():
    p = _surface_plot()
    p.scatter(np.random.default_rng(2).standard_normal((4, 3)) * 3000)


def test_scatter_with_mask_and_kwargs():
    p = _surface_plot()
    data = np.random.default_rng(3).standard_normal((4, 3)) * 3000
    p.scatter(data, mask=True, color="red", s=5)


def test_plot_2d_no_mask():
    p = _surface_plot()
    p.plot(np.random.default_rng(4).standard_normal((5, 3)) * 3000, color="k")


def test_plot_2d_with_mask():
    p = _surface_plot()
    p.plot(np.random.default_rng(5).standard_normal((5, 3)) * 3000, mask=True)


def test_plot_3d_stack_no_mask():
    p = _surface_plot()
    p.plot(np.random.default_rng(6).standard_normal((3, 5, 3)) * 3000)


def test_plot_3d_stack_with_mask():
    p = _surface_plot()
    p.plot(np.random.default_rng(7).standard_normal((3, 5, 3)) * 3000, mask=True)


# --------------------------------------------------------------------------- #
# setters                                                                      #
# --------------------------------------------------------------------------- #
def test_setters():
    p = M.Plot3D()
    p.set_labels("X", "Y", "Z")
    assert p.ax.get_xlabel() == "X"

    p.set_pane_color((0.9, 0.9, 0.9, 1.0))
    p.set_labelpad(1, 2, 3)
    assert p.ax.xaxis.labelpad == 1
    p.set_tickpad(0)

    ticks = [-1000, 0, 1000]
    p.set_ticks(ticks, ticks, ticks)
    p.set_tick_multiplier(2)

    lims = (-1000, 1000)
    p.set_lims(lims, lims, lims, equal=True)
    assert p.ax.get_xlim() == pytest.approx(lims)


def test_set_lims_not_equal():
    p = M.Plot3D()
    lims = (-500, 500)
    p.set_lims(lims, lims, lims, equal=False)
    assert p.ax.get_zlim() == pytest.approx(lims)


# --------------------------------------------------------------------------- #
# antenna gain pattern (2D)                                                    #
# --------------------------------------------------------------------------- #
def test_antenna_gain_2d_defaults(antenna):
    fig, ax = plt.subplots()
    M.plot_antenna_gain_patter_2D(ax=ax, antenna=antenna)
    assert ax.get_xlabel() != ""
    assert len(ax.lines) >= 1


def test_antenna_gain_2d_custom_inputs(antenna):
    M.plot_antenna_gain_patter_2D(
        antenna=antenna,
        phi=np.linspace(-180, 180, 50),
        theta=np.array([0.0]),  # single -> no legend branch
        show_max=False,
        name="MyAntenna",
    )
    assert plt.gca().get_title() == "MyAntenna"


def test_antenna_gain_2d_multi_theta_legend(antenna):
    M.plot_antenna_gain_patter_2D(
        antenna=antenna,
        theta=np.array([0.0, 30.0, 60.0]),
    )
    assert plt.gca().get_legend() is not None


# --------------------------------------------------------------------------- #
# covariance ellipse                                                           #
# --------------------------------------------------------------------------- #
def test_plot_2d_ellipse():
    fig, ax = plt.subplots()
    plt.sca(ax)
    mu = np.array([1.0, 2.0])
    cov = np.array([[4.0, 1.0], [1.0, 9.0]])
    patch = M.plot_2d_ellipse(mu, cov, n_std=2.0, edgecolor="blue")
    assert patch in ax.patches


def test_plot_2d_ellipse_defaults():
    fig, ax = plt.subplots()
    plt.sca(ax)
    mu = np.zeros(2)
    cov = np.eye(2)
    patch = M.plot_2d_ellipse(mu, cov)
    assert patch in ax.patches


# --------------------------------------------------------------------------- #
# namespace export                                                             #
# --------------------------------------------------------------------------- #
def test_exported_via_namespace():
    assert pnt.plot.Plot3D is M.Plot3D
    assert pnt.plot.plot_2d_ellipse is M.plot_2d_ellipse

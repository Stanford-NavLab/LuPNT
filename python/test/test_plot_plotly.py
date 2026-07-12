"""Tests for pylupnt.plot._plotly (plotly plotting helpers)."""

import numpy as np
import plotly.graph_objs as go
import pytest

import pylupnt as pnt
from pylupnt.plot import _plotly as P


@pytest.fixture(scope="module")
def antenna():
    return pnt.Antenna("Block-IIF_ACE")


# --------------------------------------------------------------------------- #
# camera / limits / aspect                                                     #
# --------------------------------------------------------------------------- #
def test_set_view():
    fig = go.Figure()
    P.set_view(fig, azimuth=45, elevation=30, zoom=2.0)
    eye = fig.layout.scene.camera.eye
    assert eye.x is not None and eye.z is not None


def test_set_lims_1d():
    fig = go.Figure()
    P.set_lims(fig, [-5, 5])  # 1D -> tiled to 3 axes
    assert tuple(fig.layout.scene.xaxis.range) == (-5, 5)
    assert tuple(fig.layout.scene.zaxis.range) == (-5, 5)


def test_set_lims_2d_array():
    fig = go.Figure()
    lims = np.array([[-1, 1], [-2, 2], [-3, 3]])
    P.set_lims(fig, lims)
    assert tuple(fig.layout.scene.yaxis.range) == (-2, 2)


def test_set_equal_aspect_ratio():
    fig = go.Figure()
    P.set_lims(fig, np.array([[-1, 1], [-2, 2], [-4, 4]]))
    P.set_equal_aspect_ratio(fig)
    ar = fig.layout.scene.aspectratio
    assert fig.layout.scene.aspectmode == "manual"
    assert ar.z == pytest.approx(1.0)  # largest range -> ratio 1


# --------------------------------------------------------------------------- #
# meshgrid / bodies                                                            #
# --------------------------------------------------------------------------- #
def test_create_sphere_meshgrid():
    x, y, z = P.create_sphere_meshgrid(8, 12, radius=2.0)
    assert x.shape == (8, 12)
    assert np.allclose(np.sqrt(x**2 + y**2 + z**2), 2.0)


def test_plot_body_basic():
    fig = go.Figure()
    out = P.plot_body(fig, pnt.MOON, size_factor=40, n_colors=4, n_training_pixels=200)
    assert out is fig
    assert len(fig.data) == 1
    assert isinstance(fig.data[0], go.Mesh3d)


def test_plot_body_with_lighting_and_rotation():
    fig = go.Figure()
    P.plot_body(
        fig,
        pnt.EARTH,
        size_factor=40,
        n_colors=4,
        n_training_pixels=200,
        R_b2frame=np.eye(3),
        r_b2s=np.array([1.0, 0.0, 0.0]),
        r_body=np.array([1e6, 0.0, 0.0]),
        alpha=0.3,
    )
    assert len(fig.data) == 1


# --------------------------------------------------------------------------- #
# frames / arrows                                                              #
# --------------------------------------------------------------------------- #
def test_plot_frame_scalar_length():
    fig = go.Figure()
    P.plot_frame(fig, np.zeros(3), np.eye(3), length=2.0)
    # 3 axes * (line + cone) = 6 traces
    assert len(fig.data) == 6


def test_plot_frame_array_length():
    fig = go.Figure()
    P.plot_frame(fig, np.zeros(3), np.eye(3), length=np.array([1.0, 2.0, 3.0]))
    assert len(fig.data) == 6


def test_plot_arrow3_single():
    fig = go.Figure()
    out = P.plot_arrow3(fig, origin=np.zeros(3), direction=np.array([1.0, 0.0, 0.0]))
    assert out is fig
    assert len(fig.data) == 2  # line + cone


def test_plot_arrow3_multiple_with_color_list():
    fig = go.Figure()
    origins = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    directions = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    P.plot_arrow3(fig, origin=origins, direction=directions, color=["red", "green"])
    assert len(fig.data) == 4


def test_plot_arrow3_requires_fig():
    with pytest.raises(AssertionError):
        P.plot_arrow3(None, origin=np.zeros(3))


# --------------------------------------------------------------------------- #
# orbits                                                                       #
# --------------------------------------------------------------------------- #
def _orbit_data(n_sat=2, n_t=10):
    rng = np.random.default_rng(0)
    return rng.standard_normal((n_sat, n_t, 6)) * 1e7


def test_plot_orbits_2d_input():
    fig = go.Figure()
    rv = _orbit_data(1, 8)[0]  # 2D
    out = P.plot_orbits(fig, rv)
    assert out is fig
    assert len(fig.data) == 1


def test_plot_orbits_3d_with_marker():
    fig = go.Figure()
    rv = _orbit_data(3, 8)
    P.plot_orbits(fig, rv, t=2, marker_size=6)
    # 3 line traces + 1 marker trace
    assert len(fig.data) == 4


def test_plot_orbits_single_color_string():
    fig = go.Figure()
    rv = _orbit_data(2, 6)
    P.plot_orbits(fig, rv, color="black")
    assert len(fig.data) == 2


def test_plot_orbits_relabels_preset_tickvals():
    fig = go.Figure()
    # Preset tickvals so the tick-relabel branch (kkm) is exercised.
    ticks = [-1e6, 0.0, 1e6]
    fig.update_layout(
        scene=dict(
            xaxis=dict(tickvals=ticks),
            yaxis=dict(tickvals=ticks),
            zaxis=dict(tickvals=ticks),
        )
    )
    rv = _orbit_data(1, 6)
    P.plot_orbits(fig, rv, t=1)
    assert fig.layout.scene.xaxis.ticktext is not None


# --------------------------------------------------------------------------- #
# image quantization / triangulation                                          #
# --------------------------------------------------------------------------- #
def test_image2zvals_normalizes_uint8():
    rng = np.random.default_rng(1)
    img = (rng.random((8, 8, 3)) * 255).astype(np.uint8)
    z, colorscale = P.image2zvals(img, n_colors=4, n_training_pixels=50)
    assert z.shape == (8, 8)
    assert len(colorscale) == 4


def test_image2zvals_float_no_normalize():
    rng = np.random.default_rng(2)
    img = rng.random((6, 6, 3))  # already in [0,1], range0 <= 1
    z, colorscale = P.image2zvals(img, n_colors=3, n_training_pixels=20)
    assert z.shape == (6, 6)


def test_image2zvals_bad_ndim():
    with pytest.raises(ValueError):
        P.image2zvals(np.zeros((4, 4)), n_colors=2)


def test_image2zvals_bad_depth():
    with pytest.raises(ValueError):
        P.image2zvals(np.zeros((4, 4, 2)), n_colors=2)


def test_regular_tri():
    tri = P.regular_tri(3, 4)
    assert tri.shape[1] == 3
    assert tri.shape[0] == 2 * (3 - 1) * (4 - 1)


def test_mesh_data():
    rng = np.random.default_rng(3)
    img = (rng.random((6, 6, 3)) * 255).astype(np.uint8)
    I, J, K, intensity, colorscale = P.mesh_data(img, n_colors=4, n_training_pixels=50)
    assert len(I) == len(J) == len(K) == len(intensity)


# --------------------------------------------------------------------------- #
# scatter                                                                      #
# --------------------------------------------------------------------------- #
def test_scatter_1d():
    fig = go.Figure()
    out = P.scatter(fig, np.array([1e6, 2e6, 3e6]))
    assert out is fig
    assert len(fig.data) == 1


def test_scatter_2d():
    fig = go.Figure()
    P.scatter(fig, np.random.default_rng(4).standard_normal((5, 3)) * 1e6)
    assert len(fig.data) == 1


def test_scatter_3d_stack_with_color_list():
    fig = go.Figure()
    xyz = np.random.default_rng(5).standard_normal((3, 5, 3)) * 1e6
    P.scatter(fig, xyz, mode="lines", marker_size=2, color=["red", "green", "blue"])
    assert len(fig.data) == 3


# --------------------------------------------------------------------------- #
# antenna gain pattern (3D)                                                    #
# --------------------------------------------------------------------------- #
def test_plot_antenna_gain_pattern_from_antenna(antenna):
    fig = P.plot_antenna_gain_pattern(antenna=antenna)
    assert isinstance(fig, go.Figure)
    # 1 surface + frame (6 traces)
    assert len(fig.data) == 7
    assert any(isinstance(t, go.Surface) for t in fig.data)


def test_plot_antenna_gain_pattern_explicit_arrays(antenna):
    theta = antenna.get_theta_vector()
    phi = antenna.get_phi_vector()
    gain = antenna.get_gain_matrix()
    fig = P.plot_antenna_gain_pattern(theta=theta, phi=phi, gain=gain)
    assert isinstance(fig, go.Figure)
    assert len(fig.data) == 7


# --------------------------------------------------------------------------- #
# namespace export                                                             #
# --------------------------------------------------------------------------- #
def test_exported_via_namespace():
    assert pnt.plot.plot_orbits is P.plot_orbits
    assert pnt.plot.scatter is P.scatter

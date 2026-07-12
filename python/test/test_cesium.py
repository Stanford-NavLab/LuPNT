"""Tests for ``pylupnt.plot.cesium`` — CZML / HTML scene serialization.

No browser and no Cesium ion token are involved: the module only builds JSON /
HTML strings, so everything is asserted on the serialized output.
"""

import json

import numpy as np
import pytest

import pylupnt as pnt
from pylupnt import _pylupnt as _pnt
from pylupnt.plot import cesium
from pylupnt.plot.cesium import BODY_PRESETS, CesiumScene, frame_orientation_quaternions


# --------------------------------------------------------------------- helpers
def _circular_orbit(body_radius, alt, n=24):
    """A simple planar circular ring of positions [n, 3] in metres."""
    r = body_radius + alt
    th = np.linspace(0.0, 2 * np.pi, n)
    return np.c_[r * np.cos(th), r * np.sin(th), np.zeros(n)]


def test_iso_and_rgba_helpers():
    from datetime import datetime, timezone

    dt = datetime(2027, 3, 1, 4, 5, 6, tzinfo=timezone.utc)
    assert cesium._iso(dt) == "2027-03-01T04:05:06Z"
    assert cesium._rgba((10, 20, 30)) == [10, 20, 30, 255]
    assert cesium._rgba((10, 20, 30, 40)) == [10, 20, 30, 40]
    assert cesium._rgba_with_alpha((10, 20, 30), 99) == [10, 20, 30, 99]


def test_matrix_to_unit_quat_identity_and_norm():
    q = cesium._matrix_to_unit_quat(np.eye(3))
    np.testing.assert_allclose(q, [0, 0, 0, 1], atol=1e-12)
    # 90 deg about z
    rz = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)
    q = cesium._matrix_to_unit_quat(rz)
    np.testing.assert_allclose(np.linalg.norm(q), 1.0)


def test_matrix_to_unit_quat_all_branches():
    # Rotations chosen so each diagonal-dominant branch is exercised.
    rx = np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]], dtype=float)  # trace<0, m00 largest
    ry = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]], dtype=float)  # m11 largest
    rz = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 1]], dtype=float)  # m22 largest
    for m in (rx, ry, rz):
        q = cesium._matrix_to_unit_quat(m)
        np.testing.assert_allclose(np.linalg.norm(q), 1.0, atol=1e-9)


def test_matrix_to_unit_quat_rejects_degenerate():
    # A non-finite matrix yields a non-finite quaternion norm -> guarded ValueError.
    with pytest.raises(ValueError):
        cesium._matrix_to_unit_quat(np.full((3, 3), np.nan))


def test_frame_orientation_quaternions_shapes():
    q = frame_orientation_quaternions([0.0, 100.0, 200.0], pnt.MOON_PA, pnt.MOON_CI)
    assert q.shape == (3, 4)
    np.testing.assert_allclose(np.linalg.norm(q, axis=1), 1.0, atol=1e-9)
    # scalar input still yields a (1, 4) result
    q1 = frame_orientation_quaternions(0.0, pnt.MOON_PA, pnt.MOON_CI)
    assert q1.shape == (1, 4)


# --------------------------------------------------------------------- texture data URI
def test_texture_data_uri_missing_file_returns_none(monkeypatch, tmp_path):
    import pylupnt.core.pylupnt_utils as pu

    monkeypatch.setattr(pu, "LUPNT_DATA_PATH", str(tmp_path))
    assert cesium._texture_data_uri("does_not_exist.png") is None


def test_texture_data_uri_handles_decode_error(monkeypatch, tmp_path):
    pytest.importorskip("PIL")
    import pylupnt.core.pylupnt_utils as pu

    topo = tmp_path / "topo"
    topo.mkdir()
    (topo / "bad.png").write_text("this is not a valid image")  # PIL will raise
    monkeypatch.setattr(pu, "LUPNT_DATA_PATH", str(tmp_path))
    assert cesium._texture_data_uri("bad.png") is None


# --------------------------------------------------------------------- construction
def test_invalid_body_raises():
    with pytest.raises(ValueError):
        CesiumScene(body="PLUTO")


def test_epoch_variants_and_defaults():
    s_default = CesiumScene(body="MOON")
    assert s_default.epoch.year == 2000
    assert s_default.name == "Moon scene"
    s_str = CesiumScene(body="EARTH", epoch="2027-03-01T00:00:00Z", name="X")
    assert s_str.epoch.year == 2027 and s_str.name == "X"
    # naive datetime gets a UTC tzinfo attached
    from datetime import datetime

    s_naive = CesiumScene(body="MOON", epoch=datetime(2025, 1, 2))
    assert s_naive.epoch.tzinfo is not None
    assert s_default.fixed_frame == BODY_PRESETS["MOON"]["frame"]


# --------------------------------------------------------------------- add_satellite
def test_add_satellite_builds_packet_and_defaults():
    s = CesiumScene(body="MOON")
    pos = _circular_orbit(_pnt.R_MOON, 2000e3, n=10)
    out = s.add_satellite("SV-1", pos, label=True)
    assert out is s  # chainable
    pk = s._entities[-1]
    assert pk["id"] == "sat/SV-1"
    assert pk["position"]["referenceFrame"] == "FIXED"
    # 10 samples * 4 (t,x,y,z)
    assert len(pk["position"]["cartesian"]) == 40
    assert "label" in pk
    assert "solidColor" in pk["path"]["material"]
    assert s._satellite_ids["SV-1"] == "sat/SV-1"


def test_add_satellite_dashed_and_inertial_and_offsets():
    s = CesiumScene(body="EARTH")
    pos = _circular_orbit(_pnt.R_EARTH, 500e3, n=8)
    offs = np.linspace(0.0, 5400.0, 8)
    s.add_satellite(
        "LEO", pos, offsets_s=offs, color=(1, 2, 3), full_orbit=False,
        dashed=True, reference_frame="INERTIAL",
    )
    pk = s._entities[-1]
    assert pk["position"]["referenceFrame"] == "INERTIAL"
    assert "polylineDash" in pk["path"]["material"]
    assert pk["path"]["leadTime"] == 0.0  # full_orbit=False
    assert s._max_offset == pytest.approx(5400.0)


@pytest.mark.parametrize(
    "bad,kwargs,msg",
    [
        (np.zeros((4, 2)), {}, "shape"),
        (np.zeros((3, 3)), {"offsets_s": np.zeros(2)}, "offsets_s"),
    ],
)
def test_add_satellite_shape_errors(bad, kwargs, msg):
    s = CesiumScene(body="MOON")
    with pytest.raises(ValueError):
        s.add_satellite("bad", bad, **kwargs)


def test_add_satellite_rejects_nonfinite_and_bad_frame():
    s = CesiumScene(body="MOON")
    pos = _circular_orbit(_pnt.R_MOON, 1000e3, n=5)
    with pytest.raises(ValueError):
        s.add_satellite("nan", np.full((5, 3), np.nan))
    with pytest.raises(ValueError):
        s.add_satellite("f", pos, reference_frame="GALACTIC")


def test_auto_color_cycles_palette():
    s = CesiumScene(body="MOON")
    pos = _circular_orbit(_pnt.R_MOON, 1000e3, n=4)
    for i in range(len(cesium._PALETTE) + 2):
        s.add_satellite(f"S{i}", pos)
    # color index wrapped past the palette length
    assert s._color_idx == len(cesium._PALETTE) + 2


# --------------------------------------------------------------------- add_trajectory
def test_add_trajectory_converts_frame():
    s = CesiumScene(body="MOON")
    t = np.linspace(0.0, 3600.0, 6)
    rv = np.tile(np.array([_pnt.R_MOON + 3000e3, 0, 0, 0, 1500.0, 0]), (6, 1))
    s.add_trajectory("relay", t, rv, frame_in=pnt.MOON_CI)
    pk = s._entities[-1]
    assert pk["id"] == "sat/relay"
    assert len(pk["position"]["cartesian"]) == 24  # 6 * 4


# --------------------------------------------------------------------- add_station
def test_add_station_from_latlon_and_options():
    s = CesiumScene(body="MOON")
    out = s.add_station("Shackleton", lat=-89.9, lon=0.0, always_visible=True)
    assert out is s
    pk = s._entities[-1]
    assert pk["id"] == "station/Shackleton"
    assert len(pk["position"]["cartesian"]) == 3
    assert pk["point"]["disableDepthTestDistance"] == 1.0e12
    assert "label" in pk and pk["label"]["disableDepthTestDistance"] == 1.0e12


def test_add_station_from_cartesian_no_label():
    s = CesiumScene(body="EARTH")
    s.add_station("GS", pos_m=[_pnt.R_EARTH, 0.0, 0.0], label=False)
    pk = s._entities[-1]
    assert "label" not in pk
    assert pk["position"]["cartesian"][0] == pytest.approx(_pnt.R_EARTH)


def test_add_station_requires_position():
    s = CesiumScene(body="MOON")
    with pytest.raises(ValueError):
        s.add_station("bad")


# --------------------------------------------------------------------- add_link
def test_add_link_builds_reference_polyline():
    s = CesiumScene(body="MOON")
    pos = _circular_orbit(_pnt.R_MOON, 2000e3, n=6)
    s.add_satellite("SV-1", pos)
    s.add_station("HQ", lat=0.0, lon=0.0)
    s.add_link("SV-1", "HQ")
    pk = s._entities[-1]
    assert pk["id"] == "link/SV-1/HQ"
    refs = pk["polyline"]["positions"]["references"]
    assert refs == ["sat/SV-1#position", "station/HQ#position"]
    assert "depthFailMaterial" in pk["polyline"]  # always_visible default
    # solid variant
    s.add_link("SV-1", "HQ", dashed=False, always_visible=False)
    pk2 = s._entities[-1]
    assert "solidColor" in pk2["polyline"]["material"]
    assert "depthFailMaterial" not in pk2["polyline"]


def test_add_link_unknown_endpoints():
    s = CesiumScene(body="MOON")
    pos = _circular_orbit(_pnt.R_MOON, 2000e3, n=6)
    s.add_satellite("SV-1", pos)
    s.add_station("HQ", lat=0.0, lon=0.0)
    with pytest.raises(ValueError):
        s.add_link("ghost", "HQ")
    with pytest.raises(ValueError):
        s.add_link("SV-1", "ghost")


# --------------------------------------------------------------------- add_body_trajectory
def test_add_body_trajectory_full_options():
    s = CesiumScene(body="EARTH")
    pos = _circular_orbit(_pnt.R_EARTH, 380000e3, n=8)  # Moon-ish distance
    offs = np.linspace(0.0, 2360000.0, 8)
    quat = np.tile(np.array([0.0, 0.0, 0.0, 1.0]), (8, 1))
    s.add_body_trajectory(
        "MOON", pos, offsets_s=offs, radius=_pnt.R_MOON, orientation_quat=quat,
        texture=True, path=True, reference_frame="INERTIAL",
    )
    pk = s._entities[-1]
    assert pk["id"] == "body/MOON"
    assert "ellipsoid" in pk
    assert "orientation" in pk and pk["orientation"]["unitQuaternion"]
    assert "path" in pk and "label" in pk
    # texture=True with a real moon texture -> image material (falls back to
    # solidColor only if the texture is missing).
    assert "material" in pk["ellipsoid"]


def test_add_body_trajectory_minimal_defaults():
    s = CesiumScene(body="MOON")
    pos = _circular_orbit(0.0, 5000e3, n=5)
    s.add_body_trajectory("Blob", pos, path=False, label=False, texture=None)
    pk = s._entities[-1]
    assert "path" not in pk
    assert "label" not in pk
    assert "orientation" not in pk
    assert "solidColor" in pk["ellipsoid"]["material"]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"pos_m": np.zeros((3, 2))},
        {"pos_m": np.zeros((3, 3)), "offsets_s": np.zeros(2)},
        {"pos_m": np.full((3, 3), np.nan)},
        {"pos_m": np.zeros((3, 3)), "orientation_quat": np.zeros((3, 3))},
        {"pos_m": np.zeros((3, 3)), "orientation_quat": np.full((3, 4), np.nan)},
        {"pos_m": np.zeros((3, 3)), "reference_frame": "WARP"},
    ],
)
def test_add_body_trajectory_validation(kwargs):
    s = CesiumScene(body="MOON")
    with pytest.raises(ValueError):
        s.add_body_trajectory("B", **kwargs)


# --------------------------------------------------------------------- output: czml / html / save
def _populated_scene():
    s = CesiumScene(body="MOON", name="Relay net")
    pos = _circular_orbit(_pnt.R_MOON, 2000e3, n=12)
    s.add_satellite("SV-1", pos)
    s.add_station("HQ", lat=0.0, lon=0.0)
    s.add_link("SV-1", "HQ")
    return s


def test_to_czml_structure_with_wireframe():
    s = _populated_scene()
    czml = s.to_czml()
    assert isinstance(czml, list)
    doc = czml[0]
    assert doc["id"] == "document" and doc["name"] == "Relay net"
    assert "clock" in doc
    ids = {p["id"] for p in czml}
    assert "sat/SV-1" in ids and "station/HQ" in ids and "link/SV-1/HQ" in ids
    # wireframe polylines present (parallels + meridians)
    assert any(str(p["id"]).startswith("body/par") for p in czml)
    assert any(str(p["id"]).startswith("body/mer") for p in czml)


def test_to_czml_without_wireframe():
    s = CesiumScene(body="EARTH", wireframe=False)
    pos = _circular_orbit(_pnt.R_EARTH, 500e3, n=6)
    s.add_satellite("A", pos)
    czml = s.to_czml()
    assert not any(str(p["id"]).startswith("body/par") for p in czml)


def test_to_html_is_valid_and_substitutes_tokens():
    s = _populated_scene()
    html = s.to_html()
    assert html.lstrip().startswith("<!doctype html>")
    # all template tokens must be replaced
    for tok in ("__VER__", "__CZML__", "__CAMERA_RANGE__", "__BODY_REACH__",
                "__BODY_TEXTURE__", "__BASE_LAYER__", "__GLOBE__",
                "__BASE_COLOR__", "__ATMOSPHERE__", "__INERTIAL_VIEW__"):
        assert tok not in html
    assert cesium.CESIUM_VERSION in html
    # embedded CZML round-trips back to JSON
    marker = "const czml="
    start = html.index(marker) + len(marker)
    end = html.index(";", start)
    parsed = json.loads(html[start:end])
    assert parsed[0]["id"] == "document"


def test_to_html_moon_uses_moon_ellipsoid_globe():
    html = CesiumScene(body="MOON").to_html()
    assert "Cesium.Ellipsoid.MOON" in html


def test_to_html_no_texture_disables_base_layer():
    s = CesiumScene(body="MOON", texture=False)
    html = s.to_html()
    # base layer JS becomes the literal `false` and texture is JSON null
    assert "baseLayer:false" in html
    assert "const BODY_TEXTURE=null" in html


def test_to_html_inertial_view_flag():
    s = CesiumScene(body="EARTH", inertial_view=True)
    assert "const INERTIAL_VIEW=true" in s.to_html()


def test_save_writes_html(tmp_path):
    s = _populated_scene()
    out = s.save(tmp_path / "nested" / "scene.html")
    assert out.is_file()
    assert out.read_text().lstrip().startswith("<!doctype html>")


def test_show_writes_and_returns_iframe(tmp_path, capsys):
    pytest.importorskip("IPython")
    s = _populated_scene()
    iframe = s.show(out_dir=tmp_path)
    # default filename derived from the scene name
    assert (tmp_path / "relay_net.html").is_file()
    assert "Saved" in capsys.readouterr().out
    assert hasattr(iframe, "src")


def test_show_default_out_dir_and_explicit_filename(tmp_path, monkeypatch):
    pytest.importorskip("IPython")
    monkeypatch.chdir(tmp_path)
    s = CesiumScene(body="MOON", name="!!!")  # name sanitizes to empty -> 'scene'
    pos = _circular_orbit(_pnt.R_MOON, 1000e3, n=4)
    s.add_satellite("A", pos)
    s.show()  # out_dir defaults to cwd/cesium_scenes
    assert (tmp_path / "cesium_scenes" / "scene.html").is_file()

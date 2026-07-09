"""Interactive Cesium (CesiumJS) visualization of constellations and surface assets.

This module turns LuPNT trajectories into a self-contained `CesiumJS
<https://cesium.com/platform/cesiumjs/>`_ web scene. Trajectories are serialized as
`CZML <https://github.com/AnalyticalGraphicsInc/czml-writer/wiki/CZML-Guide>`_ (Cesium's
time-dynamic JSON format) and embedded in a small HTML page that loads CesiumJS from a CDN.

No Cesium ion account or access token is required: the central body (Earth or Moon) is a
native Cesium globe with an embedded surface texture and a lat/lon wireframe overlay, so the
only network dependency is the CDN that serves the CesiumJS library itself.

Typical use::

    scene = CesiumScene(body="MOON", name="Lunar relays", epoch=datetime(2027, 3, 1))
    scene.add_trajectory("SV-1", t_tdb, rv_mci, frame_in=pnt.MOON_CI)   # inertial -> fixed
    scene.add_station("Shackleton", lat=-89.9, lon=0.0)
    scene.show()                       # inline in Jupyter; also writes a standalone .html

Positions passed to :meth:`CesiumScene.add_satellite` are body-fixed Cartesian metres
(ECEF for Earth, ``MOON_PA`` for the Moon) — the same rotating frame Cesium renders — so
orbits trace the natural rosette and surface stations stay locked to the terrain.
"""

import base64
import io
import json
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

from pylupnt import _pylupnt as _pnt

__all__ = ["CesiumScene", "BODY_PRESETS", "frame_orientation_quaternions"]

# Default CesiumJS version served from jsDelivr. Any recent release (>= 1.100) works.
CESIUM_VERSION = "1.123"

# Per-body rendering presets: radius, body-fixed frame, texture, base/wireframe colours.
# `base_color` is the globe's solid colour (shown wherever/until the texture drapes), kept
# light so the body is clearly visible even if the embedded imagery is slow or unavailable.
BODY_PRESETS = {
    "EARTH": dict(
        radius=_pnt.R_EARTH,
        frame=_pnt.ECEF,
        body=_pnt.EARTH,
        texture="earth_surface.jpg",
        base_color=(90, 130, 180),
        texture_flip_x=True,
        texture_flip_y=True,
        wire_color=(200, 220, 255),
        moon_ellipsoid=False,
        atmosphere=True,
    ),
    "MOON": dict(
        radius=_pnt.R_MOON,
        frame=_pnt.MOON_PA,
        body=_pnt.MOON,
        texture="moon_surface.jpeg",
        base_color=(165, 165, 170),
        texture_flip_x=False,
        texture_flip_y=False,
        wire_color=(210, 220, 235),
        moon_ellipsoid=True,
        atmosphere=False,
    ),
}

# Qualitative palette (Plotly D3) used to auto-colour satellites.
_PALETTE = [
    (31, 119, 180),
    (255, 127, 14),
    (44, 160, 44),
    (214, 39, 40),
    (148, 103, 189),
    (140, 86, 75),
    (227, 119, 194),
    (127, 127, 127),
    (188, 189, 34),
    (23, 190, 207),
]


def _iso(dt):
    return dt.strftime("%Y-%m-%dT%H:%M:%SZ")


def _rgba(color, alpha=255):
    """Coerce an (r,g,b[,a]) tuple/list to a CZML [r,g,b,a] byte list."""
    c = list(color)
    if len(c) == 3:
        c = c + [alpha]
    return [int(v) for v in c]


def _rgba_with_alpha(color, alpha):
    """Coerce a colour tuple/list to CZML bytes with a specific alpha."""
    c = _rgba(color)
    c[3] = int(alpha)
    return c


def _texture_data_uri(filename, flip_x=False, flip_y=False):
    """Return a base64 ``data:`` URI for a normalized LuPNT topo texture."""
    try:
        from pylupnt.core.pylupnt_utils import LUPNT_DATA_PATH

        path = Path(LUPNT_DATA_PATH) / "topo" / filename
        if not path.is_file():
            return None
        from PIL import Image, ImageOps

        with Image.open(path) as img:
            img = ImageOps.exif_transpose(img).convert("RGB")
            if flip_x:
                img = ImageOps.mirror(img)
            if flip_y:
                img = ImageOps.flip(img)
            buf = io.BytesIO()
            img.save(buf, format="PNG", optimize=True)
        b64 = base64.b64encode(buf.getvalue()).decode("ascii")
        return f"data:image/png;base64,{b64}"
    except Exception:
        return None


def _body_texture_data_uri(body):
    preset = BODY_PRESETS[body.upper()]
    return _texture_data_uri(
        preset["texture"],
        flip_x=preset.get("texture_flip_x", False),
        flip_y=preset.get("texture_flip_y", False),
    )


def _matrix_to_unit_quat(mat):
    """Return a Cesium ``[x, y, z, w]`` unit quaternion from a 3x3 rotation matrix."""
    m = np.asarray(mat, dtype=float)
    tr = float(np.trace(m))
    if tr > 0.0:
        s = np.sqrt(tr + 1.0) * 2.0
        qw = 0.25 * s
        qx = (m[2, 1] - m[1, 2]) / s
        qy = (m[0, 2] - m[2, 0]) / s
        qz = (m[1, 0] - m[0, 1]) / s
    elif m[0, 0] > m[1, 1] and m[0, 0] > m[2, 2]:
        s = np.sqrt(1.0 + m[0, 0] - m[1, 1] - m[2, 2]) * 2.0
        qw = (m[2, 1] - m[1, 2]) / s
        qx = 0.25 * s
        qy = (m[0, 1] + m[1, 0]) / s
        qz = (m[0, 2] + m[2, 0]) / s
    elif m[1, 1] > m[2, 2]:
        s = np.sqrt(1.0 + m[1, 1] - m[0, 0] - m[2, 2]) * 2.0
        qw = (m[0, 2] - m[2, 0]) / s
        qx = (m[0, 1] + m[1, 0]) / s
        qy = 0.25 * s
        qz = (m[1, 2] + m[2, 1]) / s
    else:
        s = np.sqrt(1.0 + m[2, 2] - m[0, 0] - m[1, 1]) * 2.0
        qw = (m[1, 0] - m[0, 1]) / s
        qx = (m[0, 2] + m[2, 0]) / s
        qy = (m[1, 2] + m[2, 1]) / s
        qz = 0.25 * s
    q = np.array([qx, qy, qz, qw], dtype=float)
    n = np.linalg.norm(q)
    if n == 0.0 or not np.isfinite(n):
        raise ValueError("rotation matrix cannot be converted to a finite unit quaternion")
    return q / n


def frame_orientation_quaternions(t_tdb, frame_in, frame_out):
    """Sample LuPNT frame rotations as Cesium ``[x, y, z, w]`` quaternions.

    The returned quaternion rotates axes fixed in ``frame_in`` into ``frame_out`` at each
    supplied TDB epoch. It is useful for orienting textured body ellipsoids in inertial scenes.
    """
    times = np.atleast_1d(np.asarray(t_tdb, dtype=float))
    qs = []
    for t in times:
        rot, _ = _pnt.get_frame_rotation_translation(float(t), frame_in, frame_out)
        qs.append(_matrix_to_unit_quat(rot))
    return np.vstack(qs)


class CesiumScene:
    """Builds a CZML document and renders it as an interactive CesiumJS scene.

    Parameters
    ----------
    body : str
        Central body, ``"EARTH"`` or ``"MOON"`` (see :data:`BODY_PRESETS`).
    name : str, optional
        Scene title (shown in the Cesium info box / used as a default filename).
    epoch : datetime or str, optional
        UTC reference epoch for the animation clock. Trajectory time offsets are measured
        from this instant. Defaults to 2000-01-01. The absolute date only labels the
        timeline; it does not affect the (body-fixed) geometry.
    multiplier : float
        Initial animation speed (simulated seconds per real second).
    texture : bool
        If ``True`` (default) drape the body with LuPNT's embedded surface image. If the
        texture is unavailable, Cesium falls back to the body's light solid colour.
    wireframe : bool
        If ``True`` (default) overlay a lat/lon graticule wireframe on the body. Built from
        polylines, it renders in every environment and reads clearly as a sphere.
    cesium_version : str
        CesiumJS version to load from the CDN.
    """

    def __init__(
        self,
        body="MOON",
        name=None,
        epoch=None,
        multiplier=200.0,
        texture=True,
        wireframe=True,
        cesium_version=CESIUM_VERSION,
        inertial_view=False,
    ):
        body = body.upper()
        if body not in BODY_PRESETS:
            raise ValueError(f"body must be one of {list(BODY_PRESETS)}, got {body!r}")
        self.body = body
        self.preset = BODY_PRESETS[body]
        self.name = name or f"{body.title()} scene"
        self.multiplier = float(multiplier)
        self.cesium_version = cesium_version
        self._use_texture = bool(texture)
        self._wireframe = bool(wireframe)
        self._inertial_view = bool(inertial_view)

        if epoch is None:
            epoch = datetime(2000, 1, 1, tzinfo=timezone.utc)
        elif isinstance(epoch, str):
            epoch = datetime.fromisoformat(epoch.replace("Z", "+00:00"))
        if epoch.tzinfo is None:
            epoch = epoch.replace(tzinfo=timezone.utc)
        self.epoch = epoch

        self._entities = []  # list of CZML packet dicts
        self._max_offset = 0.0
        self._max_radius = 0.0  # largest object distance from body centre [m], for camera framing
        self._color_idx = 0
        self._satellite_ids = {}
        self._station_ids = {}

    # ------------------------------------------------------------------ helpers
    def _next_color(self):
        c = _PALETTE[self._color_idx % len(_PALETTE)]
        self._color_idx += 1
        return c

    @property
    def fixed_frame(self):
        """The body-fixed :class:`Frame` this scene renders in (ECEF or MOON_PA)."""
        return self.preset["frame"]

    # ------------------------------------------------------------------ builders
    def add_satellite(
        self,
        name,
        pos_m,
        offsets_s=None,
        color=None,
        full_orbit=True,
        label=False,
        width=1.5,
        pixel_size=7,
        dashed=False,
        dash_length=18,
        reference_frame="FIXED",
    ):
        """Add a moving satellite with a trailing (or full-orbit) path.

        Parameters
        ----------
        name : str
            Entity name (shown on hover/click and as the optional label).
        pos_m : array_like, shape (N, 3)
            Body-fixed Cartesian positions in metres (ECEF for Earth, MOON_PA for the Moon).
        offsets_s : array_like, shape (N,), optional
            Seconds from the scene epoch for each sample. Defaults to ``0, 1, ..., N-1``.
        color : tuple, optional
            ``(r, g, b[, a])`` bytes. Auto-assigned from a palette if omitted.
        full_orbit : bool
            If ``True`` the whole sampled arc is always drawn (set the sample span to one
            orbital period to show a complete closed orbit); if ``False`` a comet-like trail
            grows behind the satellite.
        label : bool
            Draw the satellite name next to its marker.
        dashed : bool
            Draw the sampled trajectory path as a dashed line.
        reference_frame : str
            CZML reference frame for ``pos_m`` (``"FIXED"`` or ``"INERTIAL"``).
        """
        reference_frame = reference_frame.upper()
        if reference_frame not in {"FIXED", "INERTIAL"}:
            raise ValueError("reference_frame must be 'FIXED' or 'INERTIAL'")
        pos = np.asarray(pos_m, dtype=float)
        if pos.ndim != 2 or pos.shape[1] != 3:
            raise ValueError(f"pos_m must have shape (N, 3), got {pos.shape}")
        n = pos.shape[0]
        if offsets_s is None:
            offsets_s = np.arange(n, dtype=float)
        offsets_s = np.asarray(offsets_s, dtype=float)
        if offsets_s.shape != (n,):
            raise ValueError("offsets_s must have shape (N,) matching pos_m")
        if not np.isfinite(pos).all():
            raise ValueError(f"pos_m for {name!r} contains non-finite values")

        color = color or self._next_color()
        self._max_offset = max(self._max_offset, float(offsets_s[-1]))
        self._max_radius = max(self._max_radius, float(np.linalg.norm(pos, axis=1).max()))

        cart = []
        for k in range(n):
            cart += [float(offsets_s[k]), float(pos[k, 0]), float(pos[k, 1]), float(pos[k, 2])]
        start = _iso(self.epoch)
        end = _iso(self.epoch + timedelta(seconds=float(offsets_s[-1])))
        # full_orbit -> keep the entire arc lit at all times via a large lead+trail window.
        trail = float(offsets_s[-1]) if full_orbit else float(offsets_s[-1]) * 0.35
        lead = float(offsets_s[-1]) if full_orbit else 0.0

        sat_id = f"sat/{name}"
        self._satellite_ids[name] = sat_id
        material = (
            {
                "polylineDash": {
                    "color": {"rgba": _rgba(color)},
                    "gapColor": {"rgba": _rgba_with_alpha(color, 35)},
                    "dashLength": float(dash_length),
                }
            }
            if dashed
            else {"solidColor": {"color": {"rgba": _rgba(color)}}}
        )

        pk = {
            "id": sat_id,
            "name": name,
            "availability": f"{start}/{end}",
            "position": {
                "interpolationAlgorithm": "LAGRANGE",
                "interpolationDegree": 5,
                "referenceFrame": reference_frame,
                "epoch": start,
                "cartesian": cart,
            },
            "path": {
                "material": material,
                "width": width,
                "leadTime": lead,
                "trailTime": trail,
                "resolution": 120,
                "show": True,
            },
            "point": {
                "pixelSize": pixel_size,
                "color": {"rgba": _rgba(color)},
                "outlineColor": {"rgba": [0, 0, 0, 255]},
                "outlineWidth": 1,
            },
        }
        if label:
            pk["label"] = {
                "text": name,
                "font": "12px sans-serif",
                "fillColor": {"rgba": _rgba(color)},
                "pixelOffset": {"cartesian2": [10, 0]},
                "scale": 0.9,
            }
        self._entities.append(pk)
        return self

    def add_trajectory(self, name, t_tdb, rv, frame_in, **kwargs):
        """Add a satellite from an inertial (or any-frame) trajectory, converting to fixed.

        Convenience wrapper that rotates ``rv`` from ``frame_in`` into this scene's body-fixed
        frame with :func:`pylupnt.convert_frame` and derives the time offsets from ``t_tdb``.

        Parameters
        ----------
        t_tdb : array_like, shape (N,)
            TDB seconds (LuPNT internal time) for each sample.
        rv : array_like, shape (N, 6) or (N, 3)
            State (or position) in ``frame_in``.
        frame_in : Frame
            Input frame, e.g. ``pnt.MOON_CI`` or ``pnt.ECI``.
        """
        t_tdb = np.asarray(t_tdb, dtype=float)
        rv = np.asarray(rv, dtype=float)
        rv_fixed = np.asarray(_pnt.convert_frame(t_tdb, rv, frame_in, self.fixed_frame))
        pos = rv_fixed[:, :3]
        offsets = t_tdb - t_tdb[0]
        return self.add_satellite(name, pos, offsets_s=offsets, **kwargs)

    def add_station(
        self,
        name,
        pos_m=None,
        lat=None,
        lon=None,
        alt=0.0,
        color=(255, 255, 255),
        label=True,
        pixel_size=10,
        always_visible=False,
        label_offset=(12, 0),
    ):
        """Add a static surface station.

        Provide either a body-fixed Cartesian position ``pos_m`` (metres), or geodetic
        ``lat``/``lon`` in **degrees** (placed on a sphere of the body radius via
        :func:`pylupnt.lat_lon_alt_to_cart`).
        """
        if pos_m is None:
            if lat is None or lon is None:
                raise ValueError("provide pos_m, or both lat and lon (degrees)")
            pos_m = _pnt.lat_lon_alt_to_cart(
                np.array([lat, lon, alt], dtype=float), self.preset["radius"]
            )
        pos_m = np.asarray(pos_m, dtype=float).reshape(3)
        self._max_radius = max(self._max_radius, float(np.linalg.norm(pos_m)))
        station_id = f"station/{name}"
        self._station_ids[name] = station_id
        pk = {
            "id": station_id,
            "name": name,
            "position": {"cartesian": [float(pos_m[0]), float(pos_m[1]), float(pos_m[2])]},
            "point": {
                "pixelSize": pixel_size,
                "color": {"rgba": _rgba(color)},
                "outlineColor": {"rgba": [255, 255, 255, 255]},
                "outlineWidth": 1.5,
            },
        }
        if always_visible:
            pk["point"]["disableDepthTestDistance"] = 1.0e12
        if label:
            pk["label"] = {
                "text": name,
                "font": "13px sans-serif",
                "fillColor": {"rgba": [255, 255, 255, 255]},
                "pixelOffset": {"cartesian2": [float(label_offset[0]), float(label_offset[1])]},
                "showBackground": True,
                "backgroundColor": {"rgba": [0, 0, 0, 150]},
                "scale": 0.9,
            }
            if always_visible:
                pk["label"]["disableDepthTestDistance"] = 1.0e12
        self._entities.append(pk)
        return self

    def add_link(
        self,
        satellite,
        station,
        color=(255, 255, 255, 130),
        width=1.0,
        dashed=True,
        dash_length=14,
        always_visible=True,
    ):
        """Add a dynamic line from a satellite entity to a station entity.

        ``satellite`` and ``station`` are the names previously passed to
        :meth:`add_satellite` / :meth:`add_trajectory` and :meth:`add_station`.
        """
        if satellite not in self._satellite_ids:
            raise ValueError(f"unknown satellite {satellite!r}")
        if station not in self._station_ids:
            raise ValueError(f"unknown station {station!r}")

        material = (
            {
                "polylineDash": {
                    "color": {"rgba": _rgba(color)},
                    "gapColor": {"rgba": _rgba_with_alpha(color, 20)},
                    "dashLength": float(dash_length),
                }
            }
            if dashed
            else {"solidColor": {"color": {"rgba": _rgba(color)}}}
        )
        pk = {
            "id": f"link/{satellite}/{station}",
            "name": f"{satellite} to {station}",
            "polyline": {
                "positions": {
                    "references": [
                        f"{self._satellite_ids[satellite]}#position",
                        f"{self._station_ids[station]}#position",
                    ]
                },
                "arcType": "NONE",
                "width": width,
                "material": material,
            },
        }
        if always_visible:
            pk["polyline"]["depthFailMaterial"] = material
        self._entities.append(pk)
        return self

    def add_body_trajectory(
        self,
        name,
        pos_m,
        offsets_s=None,
        radius=None,
        color=(180, 180, 185),
        label=True,
        path=True,
        path_color=None,
        path_width=1.0,
        reference_frame="FIXED",
        texture=None,
        orientation_quat=None,
    ):
        """Add a moving spherical body, such as the Moon in an Earth-centred scene."""
        reference_frame = reference_frame.upper()
        if reference_frame not in {"FIXED", "INERTIAL"}:
            raise ValueError("reference_frame must be 'FIXED' or 'INERTIAL'")
        pos = np.asarray(pos_m, dtype=float)
        if pos.ndim != 2 or pos.shape[1] != 3:
            raise ValueError(f"pos_m must have shape (N, 3), got {pos.shape}")
        n = pos.shape[0]
        if offsets_s is None:
            offsets_s = np.arange(n, dtype=float)
        offsets_s = np.asarray(offsets_s, dtype=float)
        if offsets_s.shape != (n,):
            raise ValueError("offsets_s must have shape (N,) matching pos_m")
        if not np.isfinite(pos).all():
            raise ValueError(f"pos_m for {name!r} contains non-finite values")
        if orientation_quat is not None:
            orientation_quat = np.asarray(orientation_quat, dtype=float)
            if orientation_quat.shape != (n, 4):
                raise ValueError("orientation_quat must have shape (N, 4) matching pos_m")
            if not np.isfinite(orientation_quat).all():
                raise ValueError(f"orientation_quat for {name!r} contains non-finite values")
        if radius is None:
            radius = float(self.preset["radius"])
        radius = float(radius)

        self._max_offset = max(self._max_offset, float(offsets_s[-1]))
        self._max_radius = max(self._max_radius, float(np.linalg.norm(pos, axis=1).max() + radius))
        cart = []
        for k in range(n):
            cart += [float(offsets_s[k]), float(pos[k, 0]), float(pos[k, 1]), float(pos[k, 2])]
        quat = []
        if orientation_quat is not None:
            for k in range(n):
                quat += [
                    float(offsets_s[k]),
                    float(orientation_quat[k, 0]),
                    float(orientation_quat[k, 1]),
                    float(orientation_quat[k, 2]),
                    float(orientation_quat[k, 3]),
                ]
        start = _iso(self.epoch)
        end = _iso(self.epoch + timedelta(seconds=float(offsets_s[-1])))

        material = {"solidColor": {"color": {"rgba": _rgba(color)}}}
        if texture is not None:
            tex_body = name if texture is True else str(texture)
            tex = _body_texture_data_uri(tex_body)
            if tex:
                material = {"image": {"image": tex, "repeat": {"cartesian2": [1.0, 1.0]}}}

        body_id = f"body/{name}"
        pk = {
            "id": body_id,
            "name": name,
            "availability": f"{start}/{end}",
            "position": {
                "interpolationAlgorithm": "LAGRANGE",
                "interpolationDegree": 5,
                "referenceFrame": reference_frame,
                "epoch": start,
                "cartesian": cart,
            },
            "ellipsoid": {
                "radii": {"cartesian": [radius, radius, radius]},
                "material": material,
                "outline": True,
                "outlineColor": {"rgba": _rgba_with_alpha(color, 120)},
                "slicePartitions": 48,
                "stackPartitions": 32,
            },
        }
        if quat:
            pk["orientation"] = {
                "interpolationAlgorithm": "LINEAR",
                "epoch": start,
                "unitQuaternion": quat,
            }
        if path:
            pc = path_color if path_color is not None else _rgba_with_alpha(color, 160)
            pk["path"] = {
                "material": {"solidColor": {"color": {"rgba": _rgba(pc)}}},
                "width": path_width,
                "leadTime": float(offsets_s[-1]),
                "trailTime": float(offsets_s[-1]),
                "resolution": 600,
                "show": True,
            }
        if label:
            pk["label"] = {
                "text": name,
                "font": "14px sans-serif",
                "fillColor": {"rgba": [255, 255, 255, 255]},
                "pixelOffset": {"cartesian2": [14, 0]},
                "showBackground": True,
                "backgroundColor": {"rgba": [0, 0, 0, 150]},
            }
        self._entities.append(pk)
        return self

    # ------------------------------------------------------------------ output
    def _wireframe_packets(self):
        """Lat/lon graticule of the body as polyline entities (a wireframe globe).

        This is drawn *in addition* to the textured Cesium globe (see :meth:`to_html`). The
        globe gives the solid textured planet; the wireframe is built from polylines — the one
        primitive that renders in every environment — so a clear spherical body is always
        visible even where the globe's tiled imagery has not (yet) loaded.
        """
        r = float(self.preset["radius"])
        col = _rgba(self.preset["wire_color"], 170)
        n = 72  # samples per line
        packets = []

        def line(pid, pts, width):
            cart = [float(v) for p in pts for v in p]
            packets.append(
                {
                    "id": pid,
                    "polyline": {
                        "positions": {"cartesian": cart},
                        "width": width,
                        "arcType": "NONE",
                        "material": {"solidColor": {"color": {"rgba": col}}},
                    },
                }
            )

        lons = np.linspace(0.0, 2 * np.pi, n)
        for lat_deg in range(-75, 90, 15):  # parallels
            lat = np.deg2rad(lat_deg)
            pts = (
                r
                * np.c_[
                    np.cos(lat) * np.cos(lons),
                    np.cos(lat) * np.sin(lons),
                    np.full_like(lons, np.sin(lat)),
                ]
            )
            line(f"body/par{lat_deg}", pts, 2.5 if lat_deg == 0 else 1.0)
        lats = np.linspace(-np.pi / 2, np.pi / 2, n)
        for lon_deg in range(0, 360, 15):  # meridians
            lon = np.deg2rad(lon_deg)
            pts = r * np.c_[np.cos(lats) * np.cos(lon), np.cos(lats) * np.sin(lon), np.sin(lats)]
            line(f"body/mer{lon_deg}", pts, 1.0)
        return packets

    def to_czml(self):
        """Return the scene as a list of CZML packets (document, body wireframe, entities).

        The solid textured body is the Cesium globe itself (configured in :meth:`to_html`); a
        polyline wireframe of the body is added here as a robust always-visible outline.
        (A planet-sized ellipsoid *entity* at the frame origin does not render reliably even on
        real GPUs, which is why the body is the native globe plus this wireframe.)
        """
        span = max(self._max_offset, 60.0)
        start, end = _iso(self.epoch), _iso(self.epoch + timedelta(seconds=span))
        doc = {
            "id": "document",
            "name": self.name,
            "version": "1.0",
            "clock": {
                "interval": f"{start}/{end}",
                "currentTime": start,
                "multiplier": self.multiplier,
                "range": "LOOP_STOP",
                "step": "SYSTEM_CLOCK_MULTIPLIER",
            },
        }
        body = self._wireframe_packets() if self._wireframe else []
        return [doc] + body + self._entities

    def to_html(self):
        """Return a standalone HTML document (CesiumJS from CDN, texture embedded).

        The body is rendered as the native Cesium globe: for the Moon the globe uses
        ``Ellipsoid.MOON``; the surface texture is draped via a single embedded imagery tile
        (no Cesium ion token). Globe lighting is disabled so the whole body is evenly visible.
        """
        # Frame the whole scene from a fixed distance so nothing depends on async zoom-to.
        reach = max(self._max_radius, float(self.preset["radius"]))
        camera_range = 3.0 * reach

        texture = _body_texture_data_uri(self.body) if self._use_texture else None
        if texture:
            base_layer_js = (
                "Cesium.ImageryLayer.fromProviderAsync("
                "Cesium.SingleTileImageryProvider.fromUrl(BODY_TEXTURE))"
            )
        else:
            base_layer_js = "false"
        globe_js = (
            "new Cesium.Globe(Cesium.Ellipsoid.MOON)"
            if self.preset["moon_ellipsoid"]
            else "undefined"
        )

        html = _HTML_TEMPLATE
        html = html.replace("__VER__", self.cesium_version)
        html = html.replace("__CZML__", json.dumps(self.to_czml()))
        html = html.replace("__CAMERA_RANGE__", repr(float(camera_range)))
        html = html.replace("__BODY_REACH__", repr(float(reach)))
        html = html.replace("__BODY_TEXTURE__", json.dumps(texture))
        html = html.replace("__BASE_LAYER__", base_layer_js)
        html = html.replace("__GLOBE__", globe_js)
        html = html.replace("__BASE_COLOR__", json.dumps(_rgba(self.preset["base_color"])))
        html = html.replace("__ATMOSPHERE__", "true" if self.preset["atmosphere"] else "false")
        html = html.replace("__INERTIAL_VIEW__", "true" if self._inertial_view else "false")
        return html

    def save(self, path):
        """Write the standalone HTML scene to ``path`` and return the :class:`Path`."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.to_html())
        return path

    def show(self, filename=None, out_dir=None, height=520):
        """Render inline in Jupyter and write a standalone ``.html`` alongside.

        Returns an :class:`IPython.display.IFrame`. If your notebook front-end sandboxes
        iframes and the inline view is blank, open the saved ``.html`` in a browser.
        """
        from IPython.display import IFrame

        if out_dir is None:
            out_dir = Path.cwd() / "cesium_scenes"
        out_dir = Path(out_dir)
        if filename is None:
            safe = "".join(c if c.isalnum() else "_" for c in self.name).strip("_").lower()
            filename = f"{safe or 'scene'}.html"
        path = self.save(out_dir / filename)
        try:
            rel = path.relative_to(Path.cwd())
            src = str(rel)
        except ValueError:
            src = str(path)
        print(f"Saved {path}  ({len(self._entities)} entities, {path.stat().st_size/1024:.0f} kB)")
        return IFrame(src=src, width="100%", height=height)


_HTML_TEMPLATE = r"""<!doctype html><html><head><meta charset="utf-8">
<script>window.CESIUM_BASE_URL='https://cdn.jsdelivr.net/npm/cesium@__VER__/Build/Cesium/';</script>
<script src="https://cdn.jsdelivr.net/npm/cesium@__VER__/Build/Cesium/Cesium.js"></script>
<link href="https://cdn.jsdelivr.net/npm/cesium@__VER__/Build/Cesium/Widgets/widgets.css" rel="stylesheet">
<style>html,body,#c{width:100%;height:100%;margin:0;padding:0;overflow:hidden;background:#000}</style>
</head><body><div id="c"></div><script>
Cesium.Ion.defaultAccessToken='';                       // no ion token needed
const czml=__CZML__;
const CAMERA_RANGE=__CAMERA_RANGE__;
const BODY_REACH=__BODY_REACH__;
const BODY_TEXTURE=__BODY_TEXTURE__;
const INERTIAL_VIEW=__INERTIAL_VIEW__;
// The central body is the native Cesium globe (Moon uses the Moon ellipsoid); its surface is
// draped with one embedded imagery tile so no ion imagery/terrain is fetched.
const viewer=new Cesium.Viewer('c',{
  baseLayer:__BASE_LAYER__,
  globe:__GLOBE__,
  baseLayerPicker:false,geocoder:false,homeButton:false,sceneModePicker:false,
  navigationHelpButton:false,fullscreenButton:true,
  animation:true,timeline:true,infoBox:true,selectionIndicator:true,
  automaticallyTrackDataSourceClocks:true,shouldAnimate:true});
window.viewer=viewer;                                  // useful when debugging saved scenes
viewer.scene.backgroundColor=Cesium.Color.BLACK;
viewer.scene.skyAtmosphere.show=__ATMOSPHERE__;
viewer.scene.globe.showGroundAtmosphere=__ATMOSPHERE__;
viewer.scene.globe.enableLighting=false;      // evenly-lit body (no day/night terminator)
// Always give the globe a light solid colour so the body is visible even before/without the
// embedded texture (Cesium's default baseColor is near-black and reads as "no body").
viewer.scene.globe.baseColor=Cesium.Color.fromBytes(...__BASE_COLOR__);
Cesium.CzmlDataSource.load(czml).then(function(ds){
  viewer.dataSources.add(ds);
  if (Cesium.defined(ds.clock)) {
    ds.clock.getValue(viewer.clock);           // adopt the CZML interval, speed and loop mode
    if (viewer.timeline) {
      viewer.timeline.zoomTo(viewer.clock.startTime,viewer.clock.stopTime);
    }
  }
  viewer.clock.shouldAnimate=true;             // press-play works / auto-animates
  if (INERTIAL_VIEW) {
    Cesium.Transforms.preloadIcrfFixed(
      new Cesium.TimeInterval({start:viewer.clock.startTime,stop:viewer.clock.stopTime})
    ).catch(function(e){console.warn('ICRF preload failed; using pseudo-fixed fallback:',e);});
    // Keep the camera in an inertial frame; the fixed globe texture then rotates underneath
    // inertial trajectories instead of appearing frozen to the screen.
    viewer.scene.postUpdate.addEventListener(function(scene,time){
      const icrfToFixed=Cesium.Transforms.computeIcrfToFixedMatrix(time) ||
                        Cesium.Transforms.computeTemeToPseudoFixedMatrix(time);
      if (Cesium.defined(icrfToFixed)) {
        const offset=Cesium.Cartesian3.clone(viewer.camera.position);
        const transform=Cesium.Matrix4.fromRotationTranslation(icrfToFixed);
        viewer.camera.lookAtTransform(transform,offset);
      }
    });
  }
  // Frame the body (centred at the origin) from a fixed distance — deterministic and keeps the
  // body centred (unlike zoomTo on the asymmetric orbit cloud).
  viewer.camera.flyToBoundingSphere(
    new Cesium.BoundingSphere(Cesium.Cartesian3.ZERO,BODY_REACH),
    {duration:0.0,offset:new Cesium.HeadingPitchRange(0.0,Cesium.Math.toRadians(-20.0),CAMERA_RANGE)});
}).catch(function(e){console.error('CesiumScene load failed:',e);});
</script></body></html>"""

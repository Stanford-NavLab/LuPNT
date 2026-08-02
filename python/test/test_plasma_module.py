"""Tests for the pure-Python parts of ``pylupnt.plasma`` (the package ``__init__``).

Covers the re-exported symbols / ``__all__`` wiring and the base-path auto-config
block that runs when the C++ plasma base path is not yet set. The C++ setter is
monkeypatched so no real global state is mutated.
"""

import importlib
import os

import pytest

import pylupnt.plasma as plasma
from pylupnt import _pylupnt as _pnt


# --------------------------------------------------------------------- re-exports
def test_kp_loader_functions_reexported():
    for name in ("update_kp", "update_kp_table", "convert_to_csv"):
        assert hasattr(plasma, name)
        assert callable(getattr(plasma, name))


def test_all_lists_kp_and_freq_symbols():
    for name in ("update_kp", "update_kp_table", "convert_to_csv", "freq_L1", "freq_L2", "freq_L5"):
        assert name in plasma.__all__


def test_km_constants_present_but_not_in_all():
    # RE/C/GM_EARTH etc. are the km-based pecsim values; deliberately excluded
    # from __all__ so they don't shadow the SI core constants.
    for name in ("RE", "C", "PI", "SECS_DAY", "RAD2DEG", "DEG2RAD", "TECU", "GM_EARTH"):
        assert hasattr(plasma, name)
        assert name not in plasma.__all__
    # C in km/s is ~299792.458 (not the SI 299792458).
    assert abs(plasma.C - 299792.458) < 1e-3


def test_pybind_symbols_reexported():
    for name in (
        "trace_ray",
        "compute_ne",
        "gcpm_v24",
        "Satellite",
        "get_plasma_base_path",
        "set_plasma_base_path",
    ):
        assert hasattr(plasma, name)
    for name in ("trace_ray", "compute_ne", "gcpm_v24", "Satellite"):
        assert name in plasma.__all__


# --------------------------------------------------------------------- base-path autoconfig
@pytest.fixture
def _restore_plasma():
    """Restore the real C++ getters/setters and a clean module after a reload."""
    orig_get = _pnt.get_plasma_base_path
    orig_set = _pnt.set_plasma_base_path
    orig_pecsim = os.environ.get("PECSIMPY_BASE_PATH")
    orig_lupnt = os.environ.get("LUPNT_DATA_PATH")
    yield
    _pnt.get_plasma_base_path = orig_get
    _pnt.set_plasma_base_path = orig_set
    for key, val in (("PECSIMPY_BASE_PATH", orig_pecsim), ("LUPNT_DATA_PATH", orig_lupnt)):
        if val is None:
            os.environ.pop(key, None)
        else:
            os.environ[key] = val
    importlib.reload(plasma)


def _reload_with_fake_base(monkeypatch, recorded):
    _pnt.get_plasma_base_path = lambda: ""  # force the auto-config branch
    _pnt.set_plasma_base_path = lambda p: recorded.append(p)
    importlib.reload(plasma)


def test_autoconfig_uses_lupnt_data_path(_restore_plasma, tmp_path, monkeypatch):
    plasma_dir = tmp_path / "plasma"
    plasma_dir.mkdir()
    monkeypatch.setenv("LUPNT_DATA_PATH", str(tmp_path))
    recorded = []
    _reload_with_fake_base(monkeypatch, recorded)
    assert recorded == [str(plasma_dir)]
    assert os.environ["PECSIMPY_BASE_PATH"] == str(plasma_dir)


def test_autoconfig_falls_back_to_file_relative_path(_restore_plasma, monkeypatch):
    # LUPNT_DATA_PATH invalid -> walk up to <root>/data/LuPNT_data/plasma.
    monkeypatch.setenv("LUPNT_DATA_PATH", "/definitely/not/a/dir")
    recorded = []
    _reload_with_fake_base(monkeypatch, recorded)
    assert len(recorded) == 1
    assert recorded[0].endswith(os.path.join("data", "LuPNT_data", "plasma"))
    assert os.path.isdir(recorded[0])

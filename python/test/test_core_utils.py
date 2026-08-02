"""Unit tests for pylupnt's pure-Python helper modules ``pylupnt.core.base`` and
``pylupnt.core.pylupnt_utils`` (previously exercised only indirectly, so uncovered)."""

import re

import numpy as np
import pytest
import yaml

import pylupnt as pnt
from pylupnt.core import base
from pylupnt.core import pylupnt_utils as pu


# --------------------------------------------------------------------------- base.get_hash
def test_get_hash_is_deterministic_8_hex():
    h = base.get_hash("a", 1, 2.5)
    assert h == base.get_hash("a", 1, 2.5)
    assert re.fullmatch(r"[0-9a-f]{8}", h)


def test_get_hash_varies_with_args():
    assert base.get_hash("a", 1) != base.get_hash("a", 2)


# --------------------------------------------------------------------------- base.convert_title_case
@pytest.mark.parametrize(
    "s,expected",
    [
        ("hello_world", "Hello World"),
        ("lunar_gnss_odts", "Lunar Gnss Odts"),
        ("single", "Single"),
        ("", ""),
    ],
)
def test_convert_title_case(s, expected):
    assert base.convert_title_case(s) == expected


# --------------------------------------------------------------------------- base.get_timestamp
def test_get_timestamp_format():
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}", base.get_timestamp())


# --------------------------------------------------------------------------- base.set_seed
def test_set_seed_makes_numpy_reproducible():
    base.set_seed(123)
    a = np.random.rand(5)
    base.set_seed(123)
    np.testing.assert_array_equal(a, np.random.rand(5))


# --------------------------------------------------------------------------- base.enum_to_list / dict
def test_enum_to_list_and_dict_skip_underscores():
    class E:
        A = 1
        B = 2
        _hidden = 3

    assert set(base.enum_to_list(E)) == {1, 2}
    d = base.enum_to_dict(E)
    assert d == {"A": 1, "B": 2}
    assert "_hidden" not in d


def test_enum_to_list_on_bound_pnt_enum():
    assert len(base.enum_to_list(pnt.Frame)) > 0


# --------------------------------------------------------------------------- base.BracketedListDumper
def test_bracketed_list_dumper_uses_flow_style():
    out = yaml.dump({"x": [1, 2, 3]}, Dumper=base.BracketedListDumper)
    assert "[1, 2, 3]" in out


# --------------------------------------------------------------------------- pylupnt_utils.normalize
def test_normalize_1d_unit_vector():
    n = pu.normalize(np.array([3.0, 4.0]))
    np.testing.assert_allclose(np.linalg.norm(n), 1.0)
    np.testing.assert_allclose(n, [0.6, 0.8])


def test_normalize_normalizes_each_row():
    n = pu.normalize(np.array([[3.0, 4.0], [0.0, 2.0]]))
    np.testing.assert_allclose(np.linalg.norm(n, axis=1), [1.0, 1.0])


# --------------------------------------------------------------------------- pylupnt_utils.timed
def test_timed_returns_result_and_nonnegative_elapsed():
    res, dt = pu.timed(lambda x: x * 2, 21)
    assert res == 42
    assert dt >= 0.0


# --------------------------------------------------------------------------- pylupnt_utils.is_notebook
def test_is_notebook_false_under_pytest():
    assert pu.is_notebook() is False


# --------------------------------------------------------------------------- pylupnt_utils.get_output_dir
def test_get_output_dir_creates_nested_dir(tmp_path, monkeypatch):
    monkeypatch.setattr(pu, "LUPNT_OUTPUT_PATH", tmp_path)
    d = pu.get_output_dir("sub/dir")
    assert d.exists() and d == tmp_path / "sub" / "dir"


# --------------------------------------------------------------------------- pylupnt_utils.find_file
def test_find_file_hit_and_miss(tmp_path):
    (tmp_path / "a").mkdir()
    f = tmp_path / "a" / "target.txt"
    f.write_text("x")
    assert pu.find_file("target.txt", path=tmp_path) == str(f)
    assert pu.find_file("nope.txt", path=tmp_path) is None

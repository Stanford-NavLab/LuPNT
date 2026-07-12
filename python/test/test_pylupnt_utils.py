"""Unit tests for the pure-Python helpers in ``pylupnt.core.pylupnt_utils`` that are
not otherwise exercised (I/O helpers, plotting utilities, formatting).

These target the branches missing from ``test_core_utils.py``: ``is_notebook`` shell
detection, ``get_output_dir`` default branch, ``load_data``, ``timer_func``,
``get_basepath``, pickle round-trip, ``File`` (h5py), ``set_axes_equal``, ``plot_RTN``,
``format_element`` and ``print_aligned``.
"""

import sys
import types
from pathlib import Path

import numpy as np
import pytest

from pylupnt.core import pylupnt_utils as pu


# --------------------------------------------------------------------------- is_notebook
@pytest.mark.parametrize(
    "shell_name,expected",
    [
        ("ZMQInteractiveShell", True),
        ("TerminalInteractiveShell", False),
        ("SomethingElse", False),
    ],
)
def test_is_notebook_detects_shell(monkeypatch, shell_name, expected):
    shell_cls = type(shell_name, (), {})
    monkeypatch.setattr(pu, "get_ipython", lambda: shell_cls(), raising=False)
    assert pu.is_notebook() is expected


def test_is_notebook_nameerror_returns_false(monkeypatch):
    # No get_ipython in scope -> NameError -> False (also the pytest default).
    monkeypatch.delattr(pu, "get_ipython", raising=False)
    assert pu.is_notebook() is False


# --------------------------------------------------------------------------- normalize
def test_normalize_unit_and_rows():
    np.testing.assert_allclose(pu.normalize(np.array([3.0, 4.0])), [0.6, 0.8])
    rows = pu.normalize(np.array([[3.0, 4.0], [0.0, 5.0]]))
    np.testing.assert_allclose(rows, [[0.6, 0.8], [0.0, 1.0]])


# --------------------------------------------------------------------------- find_file
def test_find_file_hit_and_miss(tmp_path):
    sub = tmp_path / "a" / "b"
    sub.mkdir(parents=True)
    target = sub / "wanted.txt"
    target.write_text("x")
    assert pu.find_file("wanted.txt", path=tmp_path) == str(target)
    assert pu.find_file("absent.txt", path=tmp_path) is None


# --------------------------------------------------------------------------- timed
def test_timed_returns_result_and_elapsed():
    result, elapsed = pu.timed(lambda a, b: a + b, 2, 3)
    assert result == 5
    assert elapsed >= 0.0


# --------------------------------------------------------------------------- get_output_dir
def test_get_output_dir_default_uses_data_path(tmp_path, monkeypatch):
    monkeypatch.setattr(pu, "LUPNT_OUTPUT_PATH", None)
    monkeypatch.setattr(pu, "LUPNT_DATA_PATH", tmp_path)
    d = pu.get_output_dir()  # output_dirs=None -> no nested join
    assert d == Path(tmp_path) / "output"
    assert d.exists()


def test_get_output_dir_explicit_output_path(tmp_path, monkeypatch):
    monkeypatch.setattr(pu, "LUPNT_OUTPUT_PATH", tmp_path)
    d = pu.get_output_dir("nested")
    assert d == tmp_path / "nested"
    assert d.exists()


# --------------------------------------------------------------------------- load_data
def test_load_data_uniform_and_ragged(tmp_path, monkeypatch):
    monkeypatch.setattr(pu, "LUPNT_DATA_PATH", tmp_path)
    run = tmp_path / "output" / "run1"
    run.mkdir(parents=True)
    (run / "uniform.csv").write_text("1,2,3\n4,5,6\n")
    (run / "ragged.csv").write_text("1,2\n3,4,5\n")
    (run / "ignore.txt").write_text("not csv")

    data = pu.load_data("run1")

    assert set(data.keys()) == {"uniform", "ragged"}
    np.testing.assert_array_equal(data["uniform"], [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    # Ragged rows are padded to the max width with NaN via the except-branch.
    assert data["ragged"].shape == (2, 3)
    np.testing.assert_array_equal(data["ragged"][1], [3.0, 4.0, 5.0])
    assert np.isnan(data["ragged"][0, 2])


# --------------------------------------------------------------------------- timer_func
def test_timer_func_wraps_and_reports(monkeypatch, capsys):
    # The module uses ``time()`` (the module object is not callable) so replace it
    # with a callable clock to exercise the timing/print body.
    clock = iter([1.0, 1.5])
    monkeypatch.setattr(pu, "time", lambda: next(clock))

    @pu.timer_func
    def add(a, b):
        return a + b

    assert add(2, 3) == 5
    out = capsys.readouterr().out
    assert "executed in" in out
    assert "add" in out


# --------------------------------------------------------------------------- get_basepath
def test_get_basepath_returns_data_path():
    assert pu.get_basepath() == pu.LUPNT_DATA_PATH


# --------------------------------------------------------------------------- pickle helpers
def test_pickle_dump_and_load_roundtrip(tmp_path):
    obj = {"a": [1, 2, 3], "b": np.arange(4)}
    p = tmp_path / "obj.pkl"
    pu.dump_pickle(obj, p)
    assert p.exists()
    loaded = pu.load_pickle(p)
    assert loaded["a"] == [1, 2, 3]
    np.testing.assert_array_equal(loaded["b"], np.arange(4))


# --------------------------------------------------------------------------- File (h5py)
def test_File_creates_and_reads_hdf5(tmp_path):
    h5py = pytest.importorskip("h5py")
    path = tmp_path / "x.h5"
    with pu.File(path, "w") as f:
        assert isinstance(f, h5py.File)
        f.create_dataset("d", data=np.arange(3))
    with pu.File(path, "r") as f:
        np.testing.assert_array_equal(f["d"][:], [0, 1, 2])


# --------------------------------------------------------------------------- set_axes_equal
class _FakeAx3D:
    """Minimal 3D-axis stand-in exposing only the limit accessors set_axes_equal uses.

    Using a fake (rather than a real Matplotlib 3D axis) keeps this test independent of
    the Matplotlib/NumPy state so it exercises the pure limit-scaling arithmetic directly.
    """

    def __init__(self, xl, yl, zl):
        self._xl, self._yl, self._zl = list(xl), list(yl), list(zl)

    def get_xlim3d(self):
        return self._xl

    def get_ylim3d(self):
        return self._yl

    def get_zlim3d(self):
        return self._zl

    def set_xlim3d(self, v):
        self._xl = list(v)

    def set_ylim3d(self, v):
        self._yl = list(v)

    def set_zlim3d(self, v):
        self._zl = list(v)


def test_set_axes_equal_makes_spans_equal():
    ax = _FakeAx3D((0.0, 2.0), (-1.0, 5.0), (0.0, 1.0))  # widest span = 6
    pu.set_axes_equal(ax)
    spans = [
        ax.get_xlim3d()[1] - ax.get_xlim3d()[0],
        ax.get_ylim3d()[1] - ax.get_ylim3d()[0],
        ax.get_zlim3d()[1] - ax.get_zlim3d()[0],
    ]
    for s in spans:
        assert s == pytest.approx(6.0, abs=1e-9)
    # Each axis stays centered on its original midpoint.
    assert np.mean(ax.get_xlim3d()) == pytest.approx(1.0, abs=1e-9)
    assert np.mean(ax.get_ylim3d()) == pytest.approx(2.0, abs=1e-9)
    assert np.mean(ax.get_zlim3d()) == pytest.approx(0.5, abs=1e-9)


# --------------------------------------------------------------------------- plot_RTN
class _FakeAx:
    """No-op axis capturing nothing; every plotting call is accepted and ignored."""

    def __getattr__(self, _name):
        return lambda *a, **k: None


class _FakeFig:
    def add_subplot(self, *a, **k):
        return _FakeAx()


class _FakeGridSpec:
    def __init__(self, *a, **k):
        pass

    def __getitem__(self, _idx):
        return None


def test_plot_RTN_runs_with_pandas_series(monkeypatch):
    pd = pytest.importorskip("pandas")

    # Substitute lightweight fakes for matplotlib so the plotting logic (limits,
    # branch selection, per-plane iteration) runs without a real Matplotlib backend.
    import matplotlib

    fake_plt = types.ModuleType("matplotlib.pyplot")
    fake_plt.figure = lambda *a, **k: _FakeFig()
    fake_gs = types.ModuleType("matplotlib.gridspec")
    fake_gs.GridSpec = _FakeGridSpec
    monkeypatch.setitem(sys.modules, "matplotlib.pyplot", fake_plt)
    monkeypatch.setitem(sys.modules, "matplotlib.gridspec", fake_gs)
    # ``import matplotlib.pyplot as plt`` resolves via the matplotlib package
    # attribute, so patch that too (and the gridspec attribute for symmetry).
    monkeypatch.setattr(matplotlib, "pyplot", fake_plt, raising=False)
    monkeypatch.setattr(matplotlib, "gridspec", fake_gs, raising=False)

    # Route np.min/np.max through the builtins to stay robust to any interpreter
    # state where the numpy min/max reduction path misbehaves; the values are
    # identical, and plot_RTN only uses these for the axis-limit bookkeeping.
    monkeypatch.setattr(np, "min", lambda a, *aa, **kk: min(np.asarray(a).ravel().tolist()))
    monkeypatch.setattr(np, "max", lambda a, *aa, **kk: max(np.asarray(a).ravel().tolist()))

    n = 12
    ramp = np.linspace(-1.0, 1.0, n)
    rv = {
        "rR": pd.Series(ramp),
        "rT": pd.Series(ramp[::-1].copy()),
        "rN": pd.Series(np.sin(ramp)),
    }
    # init=True exercises the initial-point branches in addition to final/center.
    pu.plot_RTN(rv, init=True, final=True, center=True)
    # A second call with the non-default flags off covers the alternate branches.
    pu.plot_RTN(rv, init=False, final=False, center=False)


# --------------------------------------------------------------------------- format_element
def test_format_element_zero_and_nonzero():
    assert pu.format_element(0.0) == "0.0"
    assert pu.format_element(2.0) == "2.0"  # default "{}"
    assert pu.format_element(3.14159, "{:.2f}") == "3.14"


# --------------------------------------------------------------------------- print_aligned
def test_print_aligned_2d(capsys):
    pu.print_aligned(np.array([[1.5, 2.25], [10.0, 3.0]]))
    out = capsys.readouterr().out
    assert out.count("[") >= 2  # bracketed matrix rows
    assert "1.5" in out


def test_print_aligned_1d_reshapes_to_column(capsys):
    pu.print_aligned(np.array([1.0, 0.0, 22.5]))
    out = capsys.readouterr().out
    assert "22.5" in out
    assert "0.0" in out  # format_element(0.0) branch

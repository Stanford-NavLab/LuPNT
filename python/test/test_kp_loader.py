"""Unit tests for ``pylupnt.plasma.kp_loader`` (Kp geomagnetic-index loader).

Network access is always monkeypatched — no real HTTP requests are made. The
loader's file-parsing / CSV-conversion logic is exercised against tiny synthetic
JSON tables written into ``tmp_path``.
"""

import json
import os

import pytest

from pylupnt.plasma import kp_loader as kl


# --------------------------------------------------------------------- _get_base_dir
def test_get_base_dir_prefers_pecsimpy(monkeypatch, tmp_path):
    monkeypatch.setenv("PECSIMPY_BASE_PATH", str(tmp_path))
    assert kl._get_base_dir() == str(tmp_path)


def test_get_base_dir_falls_back_to_lupnt_plasma(monkeypatch, tmp_path):
    monkeypatch.delenv("PECSIMPY_BASE_PATH", raising=False)
    (tmp_path / "plasma").mkdir()
    monkeypatch.setenv("LUPNT_DATA_PATH", str(tmp_path))
    assert kl._get_base_dir() == os.path.join(str(tmp_path), "plasma")


def test_get_base_dir_raises_when_nothing_available(monkeypatch, tmp_path):
    monkeypatch.delenv("PECSIMPY_BASE_PATH", raising=False)
    # LUPNT_DATA_PATH points at a dir without a `plasma` subdir.
    monkeypatch.setenv("LUPNT_DATA_PATH", str(tmp_path))
    with pytest.raises(RuntimeError):
        kl._get_base_dir()


def test_get_base_dir_raises_when_no_env(monkeypatch):
    monkeypatch.delenv("PECSIMPY_BASE_PATH", raising=False)
    monkeypatch.delenv("LUPNT_DATA_PATH", raising=False)
    with pytest.raises(RuntimeError):
        kl._get_base_dir()


# --------------------------------------------------------------------- load_table
class _FakeResp:
    def __init__(self, status_code=200, content=b"{}"):
        self.status_code = status_code
        self.content = content


def _read_json(path):
    with open(path) as f:
        return json.load(f)


def test_load_table_downloads_and_updates_last_date(monkeypatch, tmp_path, capsys):
    calls = {}

    def fake_get(url):
        calls["url"] = url
        return _FakeResp(200, b'{"ok": 1}')

    monkeypatch.setattr(kl.requests, "get", fake_get)

    kl.load_table(str(tmp_path), 2020, 12, 31)

    json_dir = tmp_path / "data" / "kp" / "json"
    # Downloaded file written with the response content.
    target = json_dir / "kp_2020.json"
    assert target.is_file()
    assert target.read_bytes() == b'{"ok": 1}'
    # last_date advanced, first_date left at sentinel (end year not < sentinel).
    assert _read_json(json_dir / "last_date.json") == {"yyyy": 2020, "mm": 12, "dd": 31}
    assert _read_json(json_dir / "first_date.json") == {"yyyy": 0, "mm": 0, "dd": 0}
    # URL includes the requested date range.
    assert "2020-01-01" in calls["url"] and "2020-12-31" in calls["url"]
    assert "Downloaded Kp index for 2020" in capsys.readouterr().out


def test_load_table_updates_first_date_branch(monkeypatch, tmp_path):
    json_dir = tmp_path / "data" / "kp" / "json"
    json_dir.mkdir(parents=True)
    # Pre-seed a first_date later than the requested year to trigger the
    # update_first_date branch, and a last_date already ahead so only first
    # date needs updating.
    (json_dir / "first_date.json").write_text(json.dumps({"yyyy": 2025, "mm": 6, "dd": 1}))
    (json_dir / "last_date.json").write_text(json.dumps({"yyyy": 2030, "mm": 12, "dd": 31}))

    monkeypatch.setattr(kl.requests, "get", lambda url: _FakeResp(200, b"{}"))
    kl.load_table(str(tmp_path), 2020, 12, 31)

    assert _read_json(json_dir / "first_date.json") == {"yyyy": 2020, "mm": 1, "dd": 1}


def test_load_table_handles_http_error(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(kl.requests, "get", lambda url: _FakeResp(404, b""))
    kl.load_table(str(tmp_path), 2021, 6, 15)

    json_dir = tmp_path / "data" / "kp" / "json"
    assert not (json_dir / "kp_2021.json").exists()
    assert "Download failed (HTTP 404)" in capsys.readouterr().out


def test_load_table_skips_when_up_to_date(monkeypatch, tmp_path, capsys):
    json_dir = tmp_path / "data" / "kp" / "json"
    json_dir.mkdir(parents=True)
    (json_dir / "first_date.json").write_text(json.dumps({"yyyy": 1990, "mm": 1, "dd": 1}))
    (json_dir / "last_date.json").write_text(json.dumps({"yyyy": 2030, "mm": 12, "dd": 31}))

    def boom(url):  # must never be called
        raise AssertionError("requests.get should not run when up to date")

    monkeypatch.setattr(kl.requests, "get", boom)
    kl.load_table(str(tmp_path), 2020, 12, 31)
    assert "already up to date" in capsys.readouterr().out


# --------------------------------------------------------------------- update_kp_table
def test_update_kp_table_iterates_years(monkeypatch):
    recorded = []
    monkeypatch.setattr(kl, "load_table", lambda base, y, m, d: recorded.append((y, m, d)))

    class _FixedNow:
        year, month, day = 1997, 3, 4

    monkeypatch.setattr(kl, "datetime", type("D", (), {"now": staticmethod(lambda: _FixedNow)}))
    kl.update_kp_table("/base", start_year=1995)

    # Full years 1995, 1996 then the partial current year.
    assert recorded == [(1995, 12, 31), (1996, 12, 31), (1997, 3, 4)]


def test_update_kp_table_resolves_base_dir(monkeypatch):
    monkeypatch.setattr(kl, "_get_base_dir", lambda: "/resolved")
    seen = []
    monkeypatch.setattr(kl, "load_table", lambda base, y, m, d: seen.append(base))

    class _FixedNow:
        year, month, day = 1996, 1, 1

    monkeypatch.setattr(kl, "datetime", type("D", (), {"now": staticmethod(lambda: _FixedNow)}))
    kl.update_kp_table(base_dir=None, start_year=1995)
    assert seen and all(b == "/resolved" for b in seen)


# --------------------------------------------------------------------- convert_to_csv
def test_convert_to_csv_parses_json(tmp_path, capsys):
    json_dir = tmp_path / "data" / "kp" / "json"
    json_dir.mkdir(parents=True)
    data = {
        "datetime": ["2020-01-01T00:00:00", "2020-01-01T03:00:00"],
        "Kp": [1.3, 2.7],
    }
    (json_dir / "kp_2020.json").write_text(json.dumps(data))
    # A non-matching file that must be skipped by the loop's continue branch.
    (json_dir / "not_kp.json").write_text(json.dumps({"x": 1}))
    (json_dir / "readme.txt").write_text("ignore me")

    kl.convert_to_csv(str(tmp_path))

    csv_path = tmp_path / "data" / "kp" / "csv" / "kp_2020.csv"
    assert csv_path.is_file()
    lines = csv_path.read_text().strip().splitlines()
    assert lines[0] == "yyyy,mm,dd,hh,kp"
    assert lines[1] == "2020,01,01,00,1.3"
    assert lines[2] == "2020,01,01,03,2.7"
    # The skipped file must not produce a CSV.
    assert not (tmp_path / "data" / "kp" / "csv" / "not_kp.csv").exists()
    assert "Converted Kp JSON files to CSV." in capsys.readouterr().out


def test_convert_to_csv_resolves_base_dir(monkeypatch, tmp_path):
    json_dir = tmp_path / "data" / "kp" / "json"
    json_dir.mkdir(parents=True)
    (json_dir / "kp_1999.json").write_text(
        json.dumps({"datetime": ["1999-12-31T21:00:00"], "Kp": [0.0]})
    )
    monkeypatch.setattr(kl, "_get_base_dir", lambda: str(tmp_path))
    kl.convert_to_csv(base_dir=None)
    assert (tmp_path / "data" / "kp" / "csv" / "kp_1999.csv").is_file()


# --------------------------------------------------------------------- update_kp
def test_update_kp_orchestrates(monkeypatch, capsys):
    order = []
    monkeypatch.setattr(kl, "update_kp_table", lambda base, sy: order.append(("table", base, sy)))
    monkeypatch.setattr(kl, "convert_to_csv", lambda base: order.append(("csv", base)))

    kl.update_kp(base_dir="/my/base", start_year=2001)
    assert order == [("table", "/my/base", 2001), ("csv", "/my/base")]
    out = capsys.readouterr().out
    assert "Plasma data directory: /my/base" in out
    assert "updated and converted to CSV successfully" in out


def test_update_kp_resolves_base_dir(monkeypatch):
    monkeypatch.setattr(kl, "_get_base_dir", lambda: "/auto/base")
    seen = {}
    monkeypatch.setattr(kl, "update_kp_table", lambda base, sy: seen.setdefault("table", base))
    monkeypatch.setattr(kl, "convert_to_csv", lambda base: seen.setdefault("csv", base))
    kl.update_kp(base_dir=None)
    assert seen == {"table": "/auto/base", "csv": "/auto/base"}

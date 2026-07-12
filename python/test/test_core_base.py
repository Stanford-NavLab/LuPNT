"""Tests for the remaining uncovered logic in ``pylupnt.core.base`` (``wait_for_key``)
and the data-download bootstrap in ``pylupnt.core.download_data``.

No real key press and no real network download ever occur — terminal I/O and HTTP
are fully monkeypatched.
"""

import importlib
import io
import os
import shutil
import sys
import types
import zipfile

import requests

from pylupnt.core import base


# --------------------------------------------------------------------- wait_for_key (Windows branch)
def test_wait_for_key_windows_branch(monkeypatch):
    got = {}
    fake_msvcrt = types.ModuleType("msvcrt")
    fake_msvcrt.getch = lambda: got.setdefault("pressed", True) or b"x"
    monkeypatch.setitem(sys.modules, "msvcrt", fake_msvcrt)

    base.wait_for_key("press please", name="T")
    assert got.get("pressed") is True


# --------------------------------------------------------------------- wait_for_key (POSIX branch)
def test_wait_for_key_posix_branch(monkeypatch):
    # Force `import msvcrt` to fail so the termios/tty fallback runs.
    monkeypatch.setitem(sys.modules, "msvcrt", None)

    import termios
    import tty

    events = []
    monkeypatch.setattr(termios, "tcgetattr", lambda fd: "OLD")
    monkeypatch.setattr(
        termios, "tcsetattr", lambda fd, when, old: events.append(("restore", old))
    )
    monkeypatch.setattr(tty, "setraw", lambda fd: events.append(("setraw", fd)))

    class _FakeStdin:
        def fileno(self):
            return 0

        def read(self, n):
            events.append(("read", n))
            return "y"

    monkeypatch.setattr(sys, "stdin", _FakeStdin())

    base.wait_for_key("hit a key", name="T")
    kinds = [e[0] for e in events]
    assert kinds == ["setraw", "read", "restore"]  # raw -> read -> restore in finally


# --------------------------------------------------------------------- download_data bootstrap
def test_download_data_triggers_and_sets_env(tmp_path, monkeypatch):
    """Reload the module with a data path lacking `ephemeris` so the download
    branch executes; all I/O is faked."""
    monkeypatch.chdir(tmp_path)
    # LUPNT_DATA_PATH exists but has no `ephemeris` subdir -> triggers the branch.
    monkeypatch.setenv("LUPNT_DATA_PATH", str(tmp_path))

    class _FakeResp:
        def __init__(self):
            self.raw = io.BytesIO(b"fake-zip-bytes")

    monkeypatch.setattr(requests, "get", lambda url, stream=False: _FakeResp())
    monkeypatch.setattr(shutil, "copyfileobj", lambda src, dst: dst.write(src.read()))

    class _FakeZip:
        def __init__(self, name, mode):
            self.name = name

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def extractall(self):
            os.makedirs(os.path.join(os.getcwd(), "LuPNT_data"), exist_ok=True)

    monkeypatch.setattr(zipfile, "ZipFile", _FakeZip)

    import pylupnt.core.download_data as dd

    removed = {}
    monkeypatch.setattr(os, "remove", lambda p: removed.setdefault("path", p))

    importlib.reload(dd)

    # The extracted data folder was created and the env var points at it.
    assert (tmp_path / "LuPNT_data").is_dir()
    assert os.environ["LUPNT_DATA_PATH"] == os.path.join(str(tmp_path), "LuPNT_data")
    # The temporary zip download was cleaned up.
    assert removed["path"] == "LuPNT_data.zip"

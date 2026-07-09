import os

# ----------------------------------------------------------------------
# 1. Core C++ module import (after gtsam)
# ----------------------------------------------------------------------
from ._pylupnt import Logger

Logger.info("Initializing", name="PyLuPNT")

if "DISPLAY" not in os.environ:
    for candidate in (":0", ":1", ":2"):
        try:
            import subprocess

            # Check if X server is running by trying xdpyinfo
            ret = subprocess.run(
                ["xdpyinfo", "-display", candidate],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            if ret.returncode == 0:
                os.environ["DISPLAY"] = candidate
                # Export DISPLAY so subprocesses inherit it
                os.putenv("DISPLAY", candidate)
                Logger.info(
                    f"No DISPLAY environment variable found, using DISPLAY={candidate}", "PyLuPNT"
                )
                break
        except Exception:
            pass
    else:
        pass

# Version information
__version__ = "0.1.0"

import warnings

Config = dict

# C++ bindings
from . import _pylupnt
from ._pylupnt import *


# Python modules (optional — fail gracefully if deps are missing)
def _try_import(module):
    import importlib

    try:
        mod = importlib.import_module(f".{module}", package=__name__)
        # pull exported names into this namespace
        names = getattr(mod, "__all__", [n for n in dir(mod) if not n.startswith("_")])
        globals().update({n: getattr(mod, n) for n in names})
    except Exception as e:
        Logger.warn(f"Skipping .{module}: {e}", name="PyLuPNT")


# `core` is cheap (utilities used throughout pylupnt) and stays eager.
_try_import("core")

# interfaces/measurements/plot/plasma pull in heavy optional deps (plotly, sklearn,
# scipy, pandas) that most users of the core C++ bindings never touch. Import them
# lazily, on first attribute access, instead of paying their cost at `import pylupnt`.
_lazy_submodules = ["interfaces", "measurements", "plot", "plasma"]


def __getattr__(name):
    while _lazy_submodules:
        _try_import(_lazy_submodules.pop(0))
        if name in globals():
            return globals()[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


# Lander descent-guidance trajectory generators (pure numpy; scipy/cvxpy are imported only
# inside the convex solver). Exposed as a submodule: `pylupnt.lander_guidance.zem_zev_trajectory`.
try:
    from . import lander_guidance  # noqa: F401
except Exception as e:  # pragma: no cover - keep `import pylupnt` robust
    Logger.warn(f"Skipping .lander_guidance: {e}", name="PyLuPNT")

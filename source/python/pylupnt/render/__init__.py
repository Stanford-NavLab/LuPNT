import logging

try:
    from ._blender import *
except ImportError:
    try:
        from _blender import *
    except ImportError:
        logging.warning("Could not import _blender module")

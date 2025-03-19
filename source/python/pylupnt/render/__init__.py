import logging

try:
    from ._blender import *
except ImportError:
    try:
        from _blender import *
    except ImportError:
        pass

try:
    # py2 py3 compatible
    from . import download_data
    from .utils import *
    from . import render
    from .math_utils import *
    from . import scenarios

    # Order matters
    from ._pylupnt import *
    from . import plot

except ImportError as e:
    # this was installed with as a python wheel
    from utils import *
    import render
    from math_utils import *
    import scenarios

    # Order matters
    from _pylupnt import *
    import plot

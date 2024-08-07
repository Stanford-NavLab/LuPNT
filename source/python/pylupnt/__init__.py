try:
    from ._pylupnt import *  # py2 py3 compatible

    from .utils import *
    from . import plot
    from . import render
    from .math_utils import *
    from . import scenarios

except ImportError as e:
    print(e)
    print("Import failed. Trying as python wheel.")

    # this was installed with as a python wheel
    from _pylupnt import *

    from utils import *
    import plot
    import render
    from math_utils import *
    import scenarios

.. LuPNT documentation master file, created by
   sphinx-quickstart on Mon Apr  3 14:18:28 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   NOTE: The python modules documentation for PURE python modules must be explicitly annotated with meta data. See below Line 51.

.. image:: _static/LuPNT_background.png
    :alt: LuPNT background
    :width: 100%
    :align: center

-----------

LuPNT: A C++/Python Library for Lunar Communications, Positioning, Navigation, and Timing (PNT)
===============================================

.. only: not latex

    Contents:

.. toctree::
    :maxdepth: 1
    :caption: Getting Started

    introduction
    builddocs


.. _tutorial_index:

.. toctree::
    :maxdepth: 1
    :caption: Tutorial

    tutorial/Python/index
    tutorial/C++/index

.. _python_api_index:

.. toctree::
    :maxdepth: 1
    :caption: Python API

    python_api/pylupnt
    python_api/pylupnt.python_code

..
    This is a comment. Please duplicate the module and meta data you want here!
    Meta data can be a string at the end that says python_only. To be used for pure python module (no pybind).
    MAKE_DOCS/python_api/pylupnt
    MAKE_DOCS/python_api/pylupnt.python_code python_only

.. _cpp_api_index:

.. toctree::
    :maxdepth: 1
    :caption: C++ API

    cpp_api/cpp_library_root

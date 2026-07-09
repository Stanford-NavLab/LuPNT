.. LuPNT documentation master file, created by
   sphinx-quickstart on Mon Apr  3 14:18:28 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   NOTE: The python modules documentation for PURE python modules must be explicitly annotated with meta data. See below Line 51.

.. image::_static/LuPNT_background.png
    :alt: LuPNT background
    :width: 100%
    :align: center

LuPNT
================
A C++/Python Library for Lunar Communications, Positioning, Navigation, and Timing (PNT)


.. only: not latex

    Contents:

.. toctree::
    :maxdepth: 1
    :caption: Getting Started

    introduction
    development
    pages/new_simulation
    builddocs
    pages/sp3_download
    pages/cross_validation

.. _tutorial_index:

.. toctree::
    :maxdepth: 1
    :caption: Tutorial

    tutorial/Python/index
    tutorial/C++/index
    pages/related_publications

.. _math_specs_index:

.. toctree::
    :maxdepth: 1
    :caption: Math Specs

    math/index

.. _python_api_index:

.. toctree::
    :maxdepth: 1
    :caption: Python API

    python_api/pylupnt
    python_api/pylupnt.plot
    python_api/pylupnt.core.pylupnt_utils
    python_api/pylupnt.interfaces.gnss_utils
    python_api/pylupnt.plasma.kp_loader

..
    Please duplicate the module and meta data here.
    Meta data can be a string at the end that says python_only. To be used for pure python module (no pybind).
    MAKE_DOCS/python_api/pylupnt
    MAKE_DOCS/python_api/pylupnt.plot python_only
    MAKE_DOCS/python_api/pylupnt.core.pylupnt_utils python_only
    MAKE_DOCS/python_api/pylupnt.interfaces.gnss_utils python_only
    MAKE_DOCS/python_api/pylupnt.plasma.kp_loader python_only

.. _cpp_api_index:

.. toctree::
    :maxdepth: 1
    :caption: C++ API

    cpp_api/cpp_library_root

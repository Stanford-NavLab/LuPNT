----------------------
pylupnt
----------------------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

What
====

A simple example of how to use
`pybind11 <https://github.com/pybind/pybind11>`__ with
`numpy <https://numpy.org/>`__.

This C++/Python library creates a ``std::vector`` of 16-bit ints,
and provides a Python interface to the contents of this vector in a
few different ways:

-  a Python
   `List <https://docs.python.org/3/tutorial/datastructures.html#more-on-lists>`__
   (copy the data)
-  a NumPy
   `ndarray <https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html>`__
   (copy the data).
-  a NumPy
   `ndarray <https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html>`__
   (move the data).

Why
===

Python Lists are great!
However, when storing many small elements of the same type,
a Numpy array is much faster and uses a lot less memory:

|Memory used vs number of elements|

|Time used vs number of elements|

How
===

The pybind11 code is in
`src/pybind11\_numpy\_example\_python.cpp <https://github.com/ssciwr/pylupnt/blob/main/src/pylupnt_python.cpp>`__.

The python project is defined in `pyproject.toml <https://github.com/ssciwr/pylupnt/blob/main/pyproject.toml>`__
and uses `scikit-build-core <https://github.com/scikit-build/scikit-build-core>`__.

Each tagged commit triggers a `GitHub action job <https://github.com/ssciwr/pylupnt/actions/workflows/pypi.yml>`__
which uses `cibuildwheel <https://cibuildwheel.readthedocs.io/>`__ to build and upload wheels to `PyPI <https://pypi.org/project/pylupnt/>`__.

The scripts used to generate the above plots are in
`scripts <https://github.com/ssciwr/pylupnt/tree/main/scripts>`__.

This repo was quickly set up using the SSC `C++ Project
Cookiecutter <https://github.com/ssciwr/cookiecutter-cpp-project>`__.

.. |License: MIT| image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT
.. |GitHub Workflow Status| image:: https://img.shields.io/github/workflow/status/lkeegan/pylupnt/CI
   :target: https://github.com/lkeegan/pylupnt/actions?query=workflow%3ACI
.. |PyPI Release| image:: https://img.shields.io/pypi/v/pylupnt.svg
   :target: https://pypi.org/project/pylupnt
.. |Documentation Status| image:: https://readthedocs.org/projects/pylupnt/badge/
   :target: https://pylupnt.readthedocs.io/
.. |Memory used vs number of elements| image:: https://raw.githubusercontent.com/ssciwr/pylupnt/main/scripts/memory.png
.. |Time used vs number of elements| image:: https://raw.githubusercontent.com/ssciwr/pylupnt/main/scripts/time.png

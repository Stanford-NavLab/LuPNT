.. _getting_started:

Getting Started
###############

This page walks through a complete LuPNT setup from a fresh machine: installing
the toolchain, building the C++ library and Python bindings, verifying the
install, registering the Jupyter kernel, and configuring the NASA Earthdata
credentials needed by the data-driven GNSS examples and tests. For a first
overview of what LuPNT is, see :doc:`introduction`; for the developer workflow
(build variants, debugging, VS Code) see :doc:`development`.

All dependencies — the C++ compiler toolchain, CMake/Ninja, Python, and every
C++/Python package — are provisioned by `Pixi <https://pixi.sh>`_ from
conda-forge, so **you do not install any compiler or library by hand**.

1. Prerequisites
================

* **Pixi** — the only tool you install manually:

  .. code-block:: bash

     curl -fsSL https://pixi.sh/install.sh | bash

  Restart the shell (or ``source`` your profile) so ``pixi`` is on ``PATH``.

* **Git** — to clone the repository.

* **Windows** — build and run inside **WSL2 (Ubuntu)**, not native Windows;
  create the ``~/.netrc`` file (Step 6) inside WSL as well.

* **Disk** — the bundled data set (SPICE kernels, gravity fields, ephemerides,
  TLEs) is ~1 GB and is downloaded automatically on the first build.

2. Clone the repository
=======================

.. code-block:: bash

   git clone https://github.com/Stanford-NavLab/LuPNT.git
   cd LuPNT

.. note::

   The optical-navigation / stereo research code under ``python/thirdparty``
   lives in git submodules. The core library and all numbered examples build
   without them; fetch them only if you need those projects, with
   ``git submodule update --init --recursive``.

3. Install the environment
==========================

From the repository root, provision the isolated conda-forge environment (this
downloads the toolchain and all dependencies into ``.pixi/`` — it does not touch
your system Python or compilers):

.. code-block:: bash

   pixi install

Every subsequent command is run through Pixi, either as ``pixi run <task>`` or
from inside an activated shell:

.. code-block:: bash

   pixi shell        # activate the environment (exit with `exit`)
   # ...or prefix individual commands:
   pixi run <task>

4. Build LuPNT
==============

.. code-block:: bash

   pixi run build       # build the C++ library (target: lupnt)
   pixi run build-py    # build and deploy the Python bindings (_pylupnt)

.. note::

   The **first** build downloads the bundled ``LuPNT_data`` set (SPICE kernels,
   gravity fields, ephemerides, EOP/TAI-UTC tables, TLEs) via
   ``cmake/FetchLuPNTData.cmake``; later builds skip it. No credentials are
   needed for this bundled data — only the live GNSS products in Step 6 require
   an Earthdata login.

Optionally build the standalone C++ tutorial programs:

.. code-block:: bash

   pixi run build-examples          # builds all_examples into build-examples-pixi/
   ./build-examples-pixi/ex1_propagate_orbit

5. Verify the install
=====================

Inside ``pixi shell`` the ``pylupnt`` package is already importable
(``PYTHONPATH`` and the data paths are set by the environment):

.. code-block:: bash

   python -c "import pylupnt as pnt; print(pnt.R_MOON)"

Run the test suites and a first example:

.. code-block:: bash

   pixi run test-py                 # Python tests (no credentials needed)
   pixi run test-cpp                # C++ tests (some GNSS tests need Earthdata — Step 6)
   python python/examples/ex1_propagate_orbit.py

6. NASA Earthdata credentials (for the GNSS examples)
=====================================================

The Earth-GNSS examples — :doc:`ex3 <tutorial/C++/ex3_gnss_interface>`,
:doc:`ex5 <tutorial/C++/ex5_gnss_measurement_sim>`,
:doc:`ex6 <tutorial/C++/ex6_gnss_odts>`, and their Python notebooks — plus a few
``pixi run test-cpp`` cases download precise SP3/BRDC/ANTEX products from
`NASA CDDIS <https://cddis.nasa.gov/>`_, which requires a free **NASA Earthdata
Login**. All non-GNSS examples run from the bundled data with no credentials.

**Step 1 — create an account** at
`urs.earthdata.nasa.gov <https://urs.earthdata.nasa.gov/>`_ ("Register for a
profile").

**Step 2 — create a** ``~/.netrc`` **file** in your home directory (on Windows,
inside WSL). Note it is ``~/.netrc``, *not* a file in the repository:

.. code-block:: bash

   touch ~/.netrc
   chmod 0600 ~/.netrc
   echo "machine urs.earthdata.nasa.gov login <your_username> password <your_password>" >> ~/.netrc

.. warning::

   The ``chmod 0600`` step is required — ``curl``/``wget`` refuse to use a
   world-readable ``.netrc``. Keep this file private and **never commit it**; it
   lives outside the repository at ``~/.netrc``, so it cannot be checked into git.

**Step 3 (alternative)** — for short local tests both loaders also read the
credentials from environment variables instead of ``~/.netrc``:

.. code-block:: bash

   export EARTHDATA_USERNAME="<your_username>"
   export EARTHDATA_PASSWORD="<your_password>"

**Step 4** — on the first scripted CDDIS access, sign in through a browser once
and approve any Earthdata application-access prompt for CDDIS.

The full download workflow (cache location, running the Python/C++ SP3 examples,
and the ex6 precompute stages) is on the :doc:`SP3 Downloads <pages/sp3_download>`
page.

7. Jupyter / notebook setup
===========================

The Python example notebooks under ``python/examples`` run in a dedicated kernel
that bakes in ``PYTHONPATH``, ``LUPNT_DATA_PATH``, and the other runtime paths,
so ``import pylupnt`` and the data-backed examples work when the kernel is
launched directly by Jupyter or VS Code (outside Pixi activation). Register it
once per machine:

.. code-block:: bash

   pixi run install-kernel     # registers the "LuPNT (pixi)" kernel

Then select the **LuPNT (pixi)** kernel in Jupyter or VS Code.

8. Next steps
=============

* :doc:`tutorial/Python/index` and :doc:`tutorial/C++/index` — the worked
  examples (start with ``ex1_propagate_orbit``).
* :doc:`How to Create a New Simulation <pages/new_simulation>` — the agent-based
  architecture and the workflow for building your own scenario.
* :doc:`math/index` — the mathematical specifications for the dynamics,
  measurement, and estimation models.

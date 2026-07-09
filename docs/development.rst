.. _development:

Development with VSCode
=======================

This page collects the developer workflow for building LuPNT, wiring up the
Python bindings, and debugging the C++/Python code from Visual Studio Code.
Day-to-day builds are driven by `Pixi <https://pixi.sh>`_, which provisions the
compiler toolchain, Python, and all C++ dependencies from conda-forge — see the
:ref:`introduction` and the project ``README.md`` for the one-time setup
(``pixi install``, ``pixi run build``, ``pixi run build-py``).

Environment
-----------

LuPNT uses Pixi for reproducible environments; you do not need to create a conda
or mamba environment by hand. The Pixi manifest ``pixi.toml`` defines both the
dependencies and the task shortcuts (``pixi run build``, ``build-py``,
``test``, ``build-docs``, ...).

Activate the environment in a shell — this exports ``LUPNT_DATA_PATH``,
``LUPNT_OUTPUT_PATH``, ``PECSIMPY_BASE_PATH``, and ``PYTHONPATH`` (pointed at
``python/``) automatically:

.. code-block:: bash

   pixi shell

or run a single command without entering a shell:

.. code-block:: bash

   pixi run <command>

To point VS Code's Python and CMake extensions at the Pixi interpreter, select
the interpreter under ``.pixi/envs/default/bin/python`` (Command Palette →
*Python: Select Interpreter*).

Building the Python bindings
----------------------------

The Pixi task ``build-py`` runs the CMake ``pylupnt-dev`` target, which builds
LuPNT and the pybind11 bindings, copies the generated ``_pylupnt.*.so`` into
``python/pylupnt/``, and generates type stubs (``.pyi``) so IDEs like VS Code
can autocomplete the compiled package:

.. code-block:: bash

   pixi run build-py          # Release bindings
   pixi run build-py-reldbg   # RelWithDebInfo bindings
   pixi run build-py-debug    # Debug bindings

Inside ``pixi shell`` the package is already importable (``PYTHONPATH`` includes
``python/``). Verify with:

.. code-block:: bash

   python -c "import pylupnt as pnt; print(pnt.R_MOON)"

Outside of Pixi activation (for example a Jupyter/VS Code kernel launched
directly), register the LuPNT kernel once so ``import pylupnt`` and the runtime
data paths resolve:

.. code-block:: bash

   pixi run install-kernel     # bakes PYTHONPATH / LUPNT_DATA_PATH into the "LuPNT (pixi)" kernel

Debugging
---------

You can debug Python and C++ together by installing the
`Python C++ Debugger <https://marketplace.visualstudio.com/items?itemName=benjamin-simmonds.pythoncpp-debug>`_
extension; its website has examples for Windows and ``gdb``.

For Apple silicon (and generally on macOS), install the
`CodeLLDB <https://marketplace.visualstudio.com/items?itemName=vadimcn.vscode-lldb>`_
extension for C++ debugging and create ``.vscode/launch.json`` with the
following configurations. Change the path to ``eigenlldb.py`` (shipped at the
repo root) to enable pretty-printing of Eigen / LuPNT matrix types.

.. code-block:: python

        {
            "configurations": [
                {
                    "name": "* C++ Attach",
                    "type": "lldb",
                    "request": "attach",
                    "pid": "",
                    "initCommands": [
                        # ***** CHANGE THIS *****
                        "command script import \"YOUR-PATH-TO-LUPNT/LuPNT/eigenlldb.py\"",
                        # ***** CHANGE THIS *****
                    ],
                },
                {
                    "name": "* Python",
                    "type": "debugpy",
                    "request": "launch",
                    "program": "${file}",
                    "cwd": "${fileDirname}",
                    "console": "integratedTerminal"
                },
                {
                    "name": "* Python/C++ Debugger",
                    "type": "pythoncpp",
                    "request": "launch",
                    "pythonLaunchName": "* Python Debugger: Current File",
                    "cppAttachName": "* Attach",
                },
            ],
        }

The ``Python/C++ Debugger`` starts a normal Python debug session and passes the
process PID to the C++ debugger so you can step across the pybind11 boundary.

To debug a pure C++ target with the CMake extension, edit ``.vscode/settings.json``:

.. code-block:: python

    {
        "cmake.debugConfig": {
            "name": "* C++ Launch",
            "type": "lldb",
            "request": "launch",
            "initCommands": [
                // ***** CHANGE THIS *****
                "command script import \"YOUR-PATH-TO-LUPNT/LuPNT/eigenlldb.py\"",
                // ***** CHANGE THIS *****
            ],
        },
    }

Pretty printing
---------------

To use ``eigenlldb.py`` for pretty-printing with CodeLLDB, install NumPy into
LLDB's Python. Open the Command Palette (``Cmd/Ctrl + Shift + P``), select
``LLDB: Command Prompt``, and run ``pip install numpy``.

In the debug console, use ``p <variable>`` to print a value or ``? <variable>``
to inspect its raw contents. For example, given

.. code-block:: cpp

        Vec3 r(4338.99, -4338.99, -0.0757297);
        Body moon = Body::Moon();

the console shows

.. code-block:: bash

        p r
        (lupnt::Vec3) (3,1) (static,static)
        [[ 4338.99     ]
         [-4338.99     ]
         [   -0.0757297]]

        ? moon
        {id:MOON, name:"MOON", GM:4902.8001180000001, R:1737.4, ...}
            id = MOON
            name = "MOON"
            GM = 4902.8001180000001
            R = 1737.4000000000001
            fixed_frame = MOON_PA
            inertial_frame = MOON_CI
            use_gravity_field = true

Pre-commit hooks
----------------

LuPNT uses `pre-commit <https://pre-commit.com/>`_ to enforce formatting
(``clang-format`` / ``cmake-format``). The hooks run automatically on
``git commit``; run them manually with:

.. code-block:: bash

   pixi run pre-commit run --all-files

A failing formatting check rejects the commit; ``git commit --no-verify``
bypasses the hook if you need to commit work in progress.

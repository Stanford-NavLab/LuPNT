.. _development:

Development with VSCode
=======================

Debugging
---------

You can debug Python and C++ by installing the `Python C++ Debugger <https://marketplace.visualstudio.com/items?itemName=benjamin-simmonds.pythoncpp-debug>`_ extension. The website provides examples for debugging in Windows and with ``gdb``.

For Apple silicon, install the `CodeLLDB <https://marketplace.visualstudio.com/items?itemName=vadimcn.vscode-lldb>`_ extension and create the file ``.vscode/launch.json`` with the following configurations.
This extension allows debugging C++ code.


Launch Configurations
---------------------

Change the path to ``eigenlldb.py`` to enable pretty printing of C++ types.

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

The ``Python/C++ Debugger`` starts your typical debugging session in python and passes the PID of your process to the C++ debugger.

Settings
--------

Edit the ``.vscode/settings.json`` file to debug all targets using CodeLLDB.
Change the path to ``eigenlldb.py`` to enable pretty printing of C++ types.

.. code-block:: python

    {
        "cmake.sourceDirectory": "${workspaceFolder}/all",
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

Pretty Printing
---------------

To use the ``eigenlldb.py`` script for pretty printing with the CodeLLDB extension, we need to install numpy within its Python distribution.
Open the VSCode Command Palette by pressing ``Command + Shift + P`` (or ``Ctrl + Shift + P``), and select ``LLDB: Command Prompt``.
Once the LLDB command prompt opens, enter ``pip install numpy``.

Debugging Variables
-------------------

In the VSCode debug console use ``p <variable-name>`` to print the name of a variables or ``? <variable-name>`` to inspect its raw contents. For example, for position vector and velocity vector would be.

.. code-block:: cpp

        Vec3 r(4338.99, -4338.99, -0.0757297);
        Body moon = Body::Moon();

.. code-block:: bash

        p position
        (lupnt::Vec3) (3,1) (static,static)
        [[ 4338.99     ]
         [-4338.99     ]
         [   -0.0757297]]

        p body
        (const lupnt::Body &) 0x00000001674e0210

        ? body
        {id:MOON, name:"MOON", GM:4902.8001180000001, ...}
            id = MOON
            name = "MOON"
            GM = 4902.8001180000001
            R = 1737.4000000000001
            n = 20
            m = 1
            fixed_frame = MOON_PA
            inertial_frame = MOON_CI
            use_gravity_field = true
            gravity_field = {n_max:1, m_max:3, n:1803047168, m:1, ...}
            [raw] = const lupnt::Body &


Building the Python Bindings
----------------------------

Run the CMake target ``pylupnt-dev`` to automatically build ``pylupnt``, i.e., LuPNT and the python bindings, copy the generated ``_pylupnt.*.so`` file to ``source/python/pylupnt``, and generate stubs for the bindings. The stubs allow IDEs like VSCode to find the contents of a compiled package such us our bindings. Finally, add ``source/python/pylupnt`` to your python path by adding the environment variable

.. code-block:: bash

        export PYTHONPATH="${PYTHONPATH}:YOUR-PATH-TO-LUPNT/LuPNT/source/python"

You can check whether it was added by fully restarting your terminal session or VSCode and executing the following command

.. code-block:: bash

        python -c "import sys; print(sys.path)"

Now you can import the package and use it in your python scripts.

.. code-block:: bash

        python -c "import pylupnt as pnt; print(pnt.R_MOON)"

Pre-Commit Hook
----------------

To ensure that the code is formatted correctly, we use the `pre-commit <https://pre-commit.com/>`_ tool. The configuration file ``.pre-commit-config.yaml`` is located in the root directory of the repository. To install the pre-commit hook, run the following commands in the root directory of the repository.

.. code-block:: bash

        pip install pre-commit
        pre-commit install

This command installs the pre-commit hook in the repository. The hook is executed before each commit and checks the code formatting. If the code is not formatted correctly, the commit is rejected. To bypass the hook, use the ``--no-verify`` option when committing.

To test the pre-commit hook, run the following command in the root directory of the repository.

.. code-block:: bash

        pre-commit run --all-files

Python Version
--------------

You can specify the python version to use in the ``.vscode/settings.json`` file.

.. code-block:: python

        "cmake.configureArgs": [
            # ***** CHANGE THIS *****
            "-DPYTHON_EXECUTABLE=YOUR-PATH-TO-PYTHON/bin/python"
            # ***** CHANGE THIS *****
        ],

The path to the python executable can be found by running the following command in the terminal.

.. code-block:: bash

        which python

Make sure you use `Clean Reconfigure All Projects` in the CMake extension to apply the changes.
You should see the following message in the CMake output.

.. code-block:: bash

        [cmake] -- pybind11 v2.13.1
        [cmake] -- Found PythonInterp: YOUR-PATH-TO-PYTHON/bin/python (found suitable version "3.9.19", minimum required is "3.7")
        [cmake] -- Found PythonLibs: YOUR-PATH-TO-PYTHON/lib/libpython3.X.dylib

It is ok to see other versions of python in the output of other libraries like Format.cmake or autodiff, as long as pybind11 finds the correct version.


Python Environment
------------------

We recommend using `micromamba <https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html>`_  to create a new environment for the project.

Conda and mamba can be used for package management, but mamba is generally preferred due to its faster performance. Mamba is a partial re-implementation of conda in C++, utilizing multi-threading and libsolv for improved speed. It still relies on some conda code and has a dependency on Python. Alternatively, micromamba is a pure C++ re-implementation with no Python dependency, making it a better choice for creating lightweight images from scratch. It is possible to `use mamba's solver from within conda <https://conda.github.io/conda-libmamba-solver/user-guide/>`_.

To create a new environment with micromamba, run the following commands in the root directory of the repository.

.. code-block:: bash

        micromamba create -n lupnt python=3.11
        micromamba activate lupnt
        micromamba install --file requirements.yml

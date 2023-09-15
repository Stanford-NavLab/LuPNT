.. _install:

Installation
============

Prerequisites
---------------------------

**Operating System**: MacOS, Windows

**Python**: >=3.7 (for Python bindings)

Preliminaries
---------------------------

.. Todo: Create a bashfile to do the installation
1. Add user path to the base of this simulator in :code:`lupnt/core/UserFilePath.h`

2. If you are using VSCode (recommended), do the additional setups as listed :ref:`here <reference>`. 

3. Install the following third party requirements

.. TODO: make this work
.. [here](#working-with-vscode)

Third Party Requirements
^^^^^^^^^^^^^^^^^^^^^^^^

+-----------------------------------------------------------------------+------------------------------------------+----------+--------------------------------------------------------------------------------------------------------------------------------------+
| Name                                                                  | Usage                                    | Version  | Installation                                                                                                                         |
+=======================================================================+==========================================+==========+======================================================================================================================================+
| `autodiff <https://github.com/autodiff/autodiff>`__                   | automatic                                | 0.6.12   | - Rename the entire folder to "autodiff"                                                                                             |
|                                                                       | differentiation                          |          | - Place it under  :code:`lupnt/3rdparty`                                                                                             |
+-----------------------------------------------------------------------+------------------------------------------+----------+--------------------------------------------------------------------------------------------------------------------------------------+
| `cspice <https://naif.jpl.nasa.gov/naif/toolkit_C.html>`__            |planetary ephemeris                       | n/a      | - Name the folder :code:`cspice`lace it under :code:`lupnt/3rdparty`                                                                 |
|                                                                       |and frame conversion                      |          | - Place it under :code:`lupnt/3rdparty`                                                                                              |
|                                                                       |                                          |          | - Move :code:`cpsice.a` and :code:`csupport.a` under :code:`cspice/lib` to :code:`cspice/`                                           |
+-----------------------------------------------------------------------+------------------------------------------+----------+--------------------------------------------------------------------------------------------------------------------------------------+
| `Eigen <https://eigen.tuxfamily.org/index.php?title=Main_Page>`__     | vector and                               | 3.4.0    | - Rename the entire folder to "Eigen"                                                                                                |
|                                                                       | matrix computation                       |          | - Place it under :code:`lupnt/3rdparty`                                                                                              |
+-----------------------------------------------------------------------+------------------------------------------+----------+--------------------------------------------------------------------------------------------------------------------------------------+
| `boost <https://www.boost.org/users/download/>`__                     |                                          | n/a      |                                                                                                                                      |
+-----------------------------------------------------------------------+------------------------------------------+----------+--------------------------------------------------------------------------------------------------------------------------------------+
| `libInterpolate <https://github.com/CD3/libInterpolate>`__            |function interpolation                    | 2.6.2    | - Rename the entire folder to "libInterpolate"                                                                                       |
|                                                                       |                                          |          | - Place it under :code:`lupnt/3rdparty`.                                                                                             |
+-----------------------------------------------------------------------+------------------------------------------+----------+--------------------------------------------------------------------------------------------------------------------------------------+
| `pybind <https://pybind11.readthedocs.io/en/stable/installing.html>`__|Python bindings                           | n/a      |                                                                                                                                      |
+-----------------------------------------------------------------------+------------------------------------------+----------+--------------------------------------------------------------------------------------------------------------------------------------+
| Ephemeris Files                                                       | ephemeris data                           | n/a      | - See `here <data/ephemeris/readme.md>`__ for instructions                                                                           |
|                                                                       |                                          |          | - You can extract the kernel files from `here <https://www.dropbox.com/sh/npgjdndt9ma3tmr/AADnjjwIsdwQwsuarLrHRF76a?dl=0>`__ as well |
|                                                                       |                                          |          | - Place the files under :code:`/data/ephemeris`                                                                                      |
+-----------------------------------------------------------------------+------------------------------------------+----------+--------------------------------------------------------------------------------------------------------------------------------------+
 
Build LuPNT Library
---------------------------
.. note::
   After the build is completed, the generated library will be located at :code:`build/lupnt/liblupnt.a`

- **Option 1** Using VScode

    - From Extensions, download :code:`CMAKE Tools`

    - Inside the .vscode directory in this project, create "settings.json" and add the following CMAKE option (This required to robustly build pybind)

        .. code-block:: bash

            {
                // cmake settings
                "cmake.configureArgs": [
                    "-DPYTHON_INCLUDE_DIRS=path1string",
                    "-DPYTHON_LIBRARIES=path2string"
                    "-DBUILD_EXAMPLES=ON"
                ]
            }

        In the code block above, replace :code:`path1string` with the output of the following command as executed in the terminal:

        .. code-block:: bash

            python3 -c "import sysconfig; print(sysconfig.get_path('include'))"


        Similarly, replace :code:`path2string` with the output of the following command:

        .. code-block:: bash

            python3 -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))"

    - Configure and Build Project from the `CMAKE` tab


- **Option 2** From Terminal 

    - You can build the project by calling

        .. code-block:: bash

            $ mkdir build 
            $ cd build
            $ cmake .. -DPYTHON_INCLUDE_DIRS=$(python -c "import sysconfig; print(sysconfig.get_path('include'))")  \
            -DPYTHON_LIBRARIES=$(python -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))") \
            -DBUILD_EXAMPLES=ON
            $ make 

        The two cmake options will add the path to the python libraries which is required to build pybind.


Install the Python LuPNT Library (pylupnt)
------------------------------------------

For developers, see [here](bindings/readme.md) for details on how to add new python bindings.

1. Run CMake to build the python bindings

2. Create and activate your local venv environment

    .. code-block:: bash

        $ python3 -m venv venv
        $ . /venv/bin/activate

3. Install the lupnt library using pip (you will need to re-run this every time you add a new function to pybind)

    .. code-block:: bash

        $ pip3 install .

4. Now you can use the pylupnt library inside your project (see the codes under :code:`examples_python/`)

    .. code-block:: bash
        
        install pylupnt as lpt

Run Unit Tests
------------------------------------------

To run the tests for the c++ codes, run the following script in the project root

    .. code-block:: bash
        
        $ ./build/test/runUnitTests


To run the tests for the Python bindings, run the following script in the project root 

    .. code-block:: bash
        
        $ python3 -m pytest test_python
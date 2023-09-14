.. _reference:

References
====================================

References for Developers
-------------------------

- `Astrodynmaics Convention and Modeling Reference for Lunar, Cisluunar, and Librartion Point Orbits <https://www.colorado.edu/faculty/bosanac/sites/default/files/attached-files/nasa_tp_20220014814_final.pdf>`__.
- `Functions in Google Test <https://qiangbo-workspace.oss-cn-shanghai.aliyuncs.com/2018-12-05-gtest-and-coverage/PlainGoogleQuickTestReferenceGuide1.pdf>`__.
- Debugging with GDB

    - Follow the process in this `video <https://www.youtube.com/watch?v=-tGSO5-eRRg>`__.
    - If the error :code:`ERROR: Unable to start debugging. Unexpected GDB output from command "-exec-run". Unable to find Mach task port for process-id 5264: (os/kern) failure (0x5).` appears when running the gdb, follow the process shown `here <https://sourceware.org/gdb/wiki/PermissionsDarwin>`__.


Testing with GMAT
^^^^^^^^^^^^^^^^^

- Some of the dynamics functions are tested by comparing outputs with the GMAT library Python API.
- See `here <https://sourceforge.net/p/gmat/git/ci/GMAT-R2020a/tree/application/api/API_README.txt>`__ for instructions on how to setup the GMAT API for function.


Working with Visual Studio Code
--------------------------------

1. This project uses the `Google C++ Style <https://google.github.io/styleguide/cppguide.html>`__.

    - Set the setting :code:`C_Cpp: Clang_format_fallback Style` to :code:`Google`.
    - Set :code:`C_Cpp: Clang_format_style` to :code:`Google` if it is not set to :code:`file`.

2. Install `prerequisites <https://code.visualstudio.com/docs/cpp/cmake-linux>`__.
3. Install the extension CodeLLDB.
4. Set up run/debug targets in the CMake window, or directly run targets.
5. Setup useful shortcuts (:code:`Cmd+K` :code:`Cmd+S`) to run your run/debug targets easily.
    
    - CMake: Run Without Debugging
    - CMake: Debug
    - Debug: Start debugging

Debugging with Eigen and autodiff objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To view Eigen and autodiff objects when debugging, modify the debugger configuration:

    .. code-block:: bash

        {
            "version": "0.2.0",
            "configurations": [
                {
                    "name": "(lldb) Launch",
                    "type": "lldb",
                    "request": "launch",
                    // Resolved by CMake Tools:
                    "program": "${command:cmake.launchTargetPath}",
                    "args": [],
                    "cwd": "${workspaceFolder}",
                    "initCommands": [
                        "command script import \"/Users/guillemcv/AppData/Eigen/eigenlldb.py\"",
                    ]
                }
            ]
        }

- A useful command for the debug console is `p` to print a variable or expression. For example, :code:`$p rv_rx_gcrf->x_`, where :code:`rv_rx_gcrf` is a :code:`std::shared_ptr<State>`, results in:
    
    .. code-block:: bash

        (autodiff::VectorXreal) $2 = ([0] = 291587.67232231156, [1] = 269354.82986367267, [2] = 76112.184704362284, [3] = -1.3616218570959222, [4] = 0.66603497742054196, [5] = 1.8768878200960224)

# LuPNT
C++/Python library for Lunar Positioning, Navigation, and Timing Analysis (LuPNT).

## Step1: Install Required Packages and Files
Todo: Create a bashfile to do the installation

- [autodiff](https://github.com/autodiff/autodiff)
  - For automatic differentiation
  - Tested with v.0.6.12
  - Rename the entire folder to "autodiff", and place it under lupnt/3rdparty
- [cspice](https://naif.jpl.nasa.gov/naif/toolkit_C.html)
  - For planetary ephemris and frame conversion
  - Name the folder `cspice` and place it under lupnt/3rdparty
  - Then move `cpsice.a` and `csupport.a` under `cspice/lib` to under `cspice/`
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
  - For vector and matrix computation
  - Tested with v3.4.0
  - Place the entire folder to "Eigen", and place it under lupnt/3rdparty
- [boost](https://www.boost.org/users/download/)
- [libInterpolate](https://github.com/CD3/libInterpolate)
  - For function interpolation
  - Tested with version 2.6.2
  - Rename the entire folder to "libInterpolate", and place it under lupnt/3rdparty
- [pybind](https://pybind11.readthedocs.io/en/stable/installing.html)
  ```
  git submodule add -b stable ../../pybind/pybind11 pybind11
  git submodule update --init
  ```
- Ephemeris Files 
  - See [here](data/ephemeris/readme.md) for instructions
  - You can extract the kernel files from [here](https://www.dropbox.com/sh/npgjdndt9ma3tmr/AADnjjwIsdwQwsuarLrHRF76a?dl=0) as well
  - Place the files under `/data/ephemeris`
- Spherical Harmonics Files
  - Place the files under `/data/spherical_harmonics`

## Step2: Do Additional Setups
1. Add the path to the base of this Simulator to UserFilePath.h
2. In the project directory, execute following command to prohibit comitting your path changes
```
git update-index --assume-unchanged lupnt/core/UserFilePath.h
```
3. If you are using VSCode (recommended), do the additional setups as listed [here](#working-with-vscode)

##  Step3: Build the LuPNT Library
After the build is completed, the generated library will be located at `build/lupnt/liblupnt.a`

### Option1: From VScode
1. From Extensions, download `CMAKE Tools`
2. Inside the .vscode directory in this project, create "settings.json" and add the following CMAKE option (This required to robustly build pybind)
```{
      // cmake settings
      "cmake.configureArgs": [
          "-DPYTHON_INCLUDE_DIRS=path1string",
          "-DPYTHON_LIBRARIES=path2string"
          "-DBUILD_EXAMPLES=ON"
      ]
    }
```
In above, replace `path1string` with the output you get when typing in the following command in terminal
```
python3 -c "import sysconfig; print(sysconfig.get_path('include'))"       
```
Similarly, replace `path2string` with the output of the following command
```
python3 -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))"
```
3. Configure and Build Project from the `CMAKE` tab

### Option2: From Terminal 
You can build the project by calling
```
$ mkdir build 
$ cd build
$ cmake .. -DPYTHON_INCLUDE_DIRS=$(python -c "import sysconfig; print(sysconfig.get_path('include'))")  \
-DPYTHON_LIBRARIES=$(python -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))") \
-DBUILD_EXAMPLES=ON
$ make 
```
The two cmake options will add the path to the python libraries which is required to build pybind.

## Step4: Install the Python LuPNT Library (pylupnt)
For developers, see [here](bindings/readme.md) for details on how to add new python bindings.
1. Run CMake to build the python bindings
2. Create and activate your local venv environment
```
$ python3 -m venv venv
$ . /venv/bin/activate
```
3. install the lupnt library using pip (you will need to re-run this every time you add a new function to pybind)
```
$ pip3 install .
```
4. Now you can use the pylupnt library inside your project (see the codes under `examples_python/`)
```
install pylupnt as lpt
```

## Step5: Run Unit Tests
To run the tests for the c++ codes, run the following script in the project root
```
$ ./build/test/runUnitTests
```

To run the tests for the Python bindings, run the following script in the project root 
```
$ python3 -m pytest test_python
```


## Appendix A: Reference for Developers
- [Astrodynmaics Convention and Modeling Reference for Lunar, Cisluunar, and Librartion Point Orbits](https://www.colorado.edu/faculty/bosanac/sites/default/files/attached-files/nasa_tp_20220014814_final.pdf)
- [Functions in Google Test](https://qiangbo-workspace.oss-cn-shanghai.aliyuncs.com/2018-12-05-gtest-and-coverage/PlainGoogleQuickTestReferenceGuide1.pdf)
- Debugging with GDB
  - Follow the process in this [video](https://www.youtube.com/watch?v=-tGSO5-eRRg)
  - If the error 'ERROR: Unable to start debugging. Unexpected GDB output from command "-exec-run". Unable to find Mach task port for process-id 5264: (os/kern) failure (0x5).' appears when running the gdb, follow the process shown [here](https://sourceware.org/gdb/wiki/PermissionsDarwin)


### Testing with GMAT
- Some of the dynamics functions are tested by comparing outputs with the [GMAT](library) python API
- See [here](https://sourceforge.net/p/gmat/git/ci/GMAT-R2020a/tree/application/api/API_README.txt) for instructions on how to setup the GMAT API for function


## Appendix B: Working with VSCODE
- This project uses the [Google C++ Style](https://google.github.io/styleguide/cppguide.html).
  - Set the setting `C_Cpp: Clang_format_fallback Style` to `Google`.
  - Set `C_Cpp: Clang_format_style` to `Google` if it is not set to `file`.
- Install [prerequisites](https://code.visualstudio.com/docs/cpp/cmake-linux)
- Install the extension CodeLLDB
- Set up run/debug targets in the CMake window, or directly run targets
- Setup useful shortcuts (Cmd+K Cmd+S) to run your run/debug targets easily
  - CMake: Run Without Debugging
  - CMake: Debug
  - Debug: Start debugging
- To view Eigen and autodiff objects when debugging, modify the debugger configuration:
```
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
                "command script import \"absolute/path/to/your/eigenlldb.py\"",
            ]
        }
    ]
}
```
- A useful command for the debug console is `p` to print a variable or expression. For example, `$p rv_rx_gcrf->x_`, where `rv_rx_gcrf` is a `std::shared_ptr<State>`, results in:
```
(autodiff::VectorXreal) $2 = ([0] = 291587.67232231156, [1] = 269354.82986367267, [2] = 76112.184704362284, [3] = -1.3616218570959222, [4] = 0.66603497742054196, [5] = 1.8768878200960224)
```

# CMake-Based Project Using lupnt

## Introduction
This example demonstrates how CMake's command `find_package` can be used to
resolve the dependency of an executable `app` on **lupnt**, a header-only
C++17 library.

The source file `main.cpp` includes the header-file `lupnt/forward.hpp` and
uses a forward mode automatic differentiation algorithm to compute the derivatives of a scalar function.

The `CMakeLists.txt` file uses the command:

```cmake
find_package(lupnt)
```

to find the **lupnt** header files. The executable target `app` is then
linked against the imported target `lupnt::lupnt`:

```cmake
target_link_libraries(app lupnt::lupnt)
```

## Building and Executing the Application
To build the application, do:

```bash
cd cmake-project
mkdir build && cd build
cmake .. -DCMAKE_PREFIX_PATH=/path/to/lupnt/install/dir
make
```

To execute the application, do:

```bash
./app
```

Note: If **lupnt** has been installed system-wide, then the CMake argument
`CMAKE_PREFIX_PATH` should not be needed. Otherwise, you will need to specify
where **lupnt** is installed in your machine. For example:

```cmake
cmake .. -DCMAKE_PREFIX_PATH=$HOME/local
```

assuming directory `$HOME/local` is where **lupnt** was installed to, which should then contain the following directory:

```
$HOME/local/include/lupnt/
```

where the **lupnt** header files are located.

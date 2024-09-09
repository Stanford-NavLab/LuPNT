# LuPNT

[![Actions Status](https://github.com/Stanford-NavLab/LuPNT/workflows/MacOS/badge.svg)](https://github.com/Stanford-NavLab/LuPNT/actions)
[![Actions Status](https://github.com/Stanford-NavLab/LuPNT/workflows/Windows/badge.svg)](https://github.com/Stanford-NavLab/LuPNT/actions)
[![Actions Status](https://github.com/Stanford-NavLab/LuPNT/workflows/Ubuntu/badge.svg)](https://github.com/Stanford-NavLab/LuPNT/actions)
[![Actions Status](https://github.com/Stanford-NavLab/LuPNT/workflows/Style/badge.svg)](https://github.com/Stanford-NavLab/LuPNT/actions)
[![Actions Status](https://github.com/Stanford-NavLab/LuPNT/workflows/Install/badge.svg)](https://github.com/Stanford-NavLab/LuPNT/actions)
[![Actions Status](https://github.com/Stanford-NavLab/LuPNT/workflows/Python/badge.svg)](https://github.com/Stanford-NavLab/LuPNT/actions)
[![Actions Status](https://github.com/Stanford-NavLab/LuPNT/workflows/Examples/badge.svg)](https://github.com/Stanford-NavLab/LuPNT/actions)
[![Documentation Status](https://github.com/Stanford-NavLab/LuPNT/workflows/Docs/badge.svg)](https://github.com/Stanford-NavLab/LuPNT/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI Release](https://img.shields.io/pypi/v/pylupnt.svg)](https://pypi.org/project/pylupnt)
[![Python Versions](https://img.shields.io/pypi/pyversions/pylupnt)](https://pypi.org/project/pylupnt)
[![codecov](https://codecov.io/gh/Stanford-NavLab/LuPNT/branch/guillemc/graph/badge.svg)](https://codecov.io/gh/Stanford-NavLab/LuPNT)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Stanford-NavLab/LuPNT/development?labpath=examples%2Fpython%2Fnotebooks%2Fex_frozen_orbits.ipynb)
[![Open In Collab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1yhHImp3hB5P8dadLlv0CcYWQWLabNfuh)


<!-- [![Conda Release](https://img.shields.io/conda/v/conda-forge/pylupnt)](https://anaconda.org/conda-forge/pylupnt) -->

<p align="center">
  <img src="https://github.com/Stanford-NavLab/LuPNT/blob/guillemc/docs/_static/LuPNT_background.png?raw=true" width="auto" />
</p>

`LuPNT` is an open-source C++/Python library for Lunar Positioning, Navigation, and Timing (PNT) Research. This project is a product of the [Stanford NAV Lab](https://navlab.stanford.edu/).
If using this project in your own work please cite the following:

```bibtex
@inproceedings{IiyamaCasadesus2023,
  title = {LuPNT: Open-Souce Simulator for Lunar Positioning, Navigation, and Timing},
  author={Iiyama, Keidai and Casadesus Vila, Guillem and Gao, Grace},
  booktitle={Proceedings of the Institute of Navigation Gnss+ conference (ION Gnss+ 2023)},
  institution = {Stanford University},
  year = {2023},
  url = {https://github.com/Stanford-NavLab/LuPNT},
}
```

## Features

- [Modern CMake practices](https://pabloariasal.github.io/2018/02/19/its-time-to-do-cmake-right/)
- Clean separation of library and executable code
- Integrated test suite
- Continuous integration via [GitHub Actions](https://help.github.com/en/actions/)
- Code coverage via [codecov](https://codecov.io)
- Code formatting enforced by [clang-format](https://clang.llvm.org/docs/ClangFormat.html) and [cmake-format](https://github.com/cheshirekow/cmake_format) via [Format.cmake](https://github.com/Stanford-NavLab/Format.cmake)
- Reproducible dependency management via [CPM.cmake](https://github.com/Stanford-NavLab/CPM.cmake)
- Installable target with automatic versioning information and header generation via [PackageProject.cmake](https://github.com/Stanford-NavLab/PackageProject.cmake)
- Automatic [documentation](https://Stanford-NavLab.github.io/LuPNT) and deployment with [Doxygen](https://www.doxygen.nl) and [GitHub Pages](https://pages.github.com)
- Support for [sanitizer tools, and more](#additional-tools)

## Usage

To cleanly separate the library and subproject code, the outer `CMakeList.txt` only defines the library itself while the tests and other subprojects are self-contained in their own directories.
During development it is usually convenient to [build all subprojects at once](#build-everything-at-once).

### Dependencies

LuPNT requires [OpenMP](https://www.openmp.org) library for multiprocessing and a data directory.
The installation scripts for MacOS, Ubuntu, and Windows can be found under `scripts`.
Note that the data directory can be place anywhere as long as its path is correctly set.
Execute the scripts before building the library.

### Build and run the examples

Use the following command to build and run the executable target.

```bash
cmake -S examples/cpp -B build/examples
cmake --build build/examples -j4
./build_examples/examples/<example-name> --help
```

### Build and run test suite

Use the following commands from the project's root directory to run the test suite.

```bash
cmake -S test -B build/test
cmake --build build/test -j4
CTEST_OUTPUT_ON_FAILURE=1 cmake --build build/test --target test

# or simply call the executable:
./build/test/LuPNTTests
```

To collect code coverage information, run CMake with the `-DENABLE_TEST_COVERAGE=1` option.

### Run clang-format

Use the following commands from the project's root directory to check and fix C++ and CMake source style.
This requires _clang-format_, _cmake-format_ and _pyyaml_ to be installed on the current system.

```bash
cmake -S test -B build/test

# view changes
cmake --build build/test --target format

# apply changes
cmake --build build/test --target fix-format
```

See [Format.cmake](https://github.com/Stanford-NavLab/Format.cmake) for details.
These dependencies can be easily installed using pip.

```bash
pip install clang-format==14.0.6 cmake_format==0.6.11 pyyaml
```

### Build the documentation

The documentation is automatically built and [published](https://Stanford-NavLab.github.io/LuPNT) whenever a [GitHub Release](https://help.github.com/en/github/administering-a-repository/managing-releases-in-a-repository) is created.
To manually build documentation, call the following command.

```bash
cmake -S documentation -B build/doc
cmake --build build/doc --target GenerateDocs
# view the docs
open build/doc/doxygen/html/index.html
```

To build the documentation locally, you will need Doxygen, jinja2 and Pygments installed on your system.

### Build everything at once

The project also includes an `all` directory that allows building all targets at the same time.
This is useful during development, as it exposes all subprojects to your IDE and avoids redundant builds of the library.

```bash
cmake -S all -B build
cmake --build build -j4

# run tests
./build/test/LuPNTTests
# format code
cmake --build build --target fix-format
# run standalone
./build/standalone/LuPNT --help
# build docs
cmake --build build --target GenerateDocs
```

### Additional tools

The test and standalone subprojects include the [tools.cmake](cmake/tools.cmake) file which is used to import additional tools on-demand through CMake configuration arguments.
The following are currently supported.

#### Sanitizers

Sanitizers can be enabled by configuring CMake with `-DUSE_SANITIZER=<Address | Memory | MemoryWithOrigins | Undefined | Thread | Leak | 'Address;Undefined'>`.

#### Static Analyzers

Static Analyzers can be enabled by setting `-DUSE_STATIC_ANALYZER=<clang-tidy | iwyu | cppcheck>`, or a combination of those in quotation marks, separated by semicolons.
By default, analyzers will automatically find configuration files such as `.clang-format`.
Additional arguments can be passed to the analyzers by setting the `CLANG_TIDY_ARGS`, `IWYU_ARGS` or `CPPCHECK_ARGS` variables.

#### Ccache

Ccache can be enabled by configuring with `-DUSE_CCACHE=<ON | OFF>`.

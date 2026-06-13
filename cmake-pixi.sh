#!/usr/bin/env bash
# cmake-pixi.sh — cmake wrapper for VSCode CMake Tools
#
# Routes every cmake invocation through `pixi run` so that compilers,
# libraries, and env vars (CONDA_PREFIX, LUPNT_DATA_PATH, …) are always
# sourced from the pixi environment, regardless of which shell or conda
# environment was active when VSCode was launched.
#
# VSCode setting:  "cmake.cmakePath": "${workspaceFolder}/cmake-pixi.sh"
#
# VSCode calls this script with the working directory set to "/" rather than
# the workspace root, so we cd to the directory that contains this script
# (= the project root, where pixi.toml lives) before calling pixi run.
cd "$(dirname "${BASH_SOURCE[0]}")"
exec ~/.pixi/bin/pixi run cmake "$@"

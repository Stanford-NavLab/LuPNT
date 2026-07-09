# Registers a Jupyter kernel named "lupnt" for the active Python interpreter so that VSCode (and any
# other Jupyter front-end) can find the pixi environment in the kernel picker without manual setup.
#
# Runs once at CMake configure time; skips silently if the kernel spec already exists. Requires
# PYTHON_EXECUTABLE to be set (done by CondaAware.cmake when pixi/conda is active). Disable with
# -DLUPNT_REGISTER_KERNEL=OFF.

option(LUPNT_REGISTER_KERNEL "Register the lupnt Jupyter kernel during CMake configure" ON)

if(NOT LUPNT_REGISTER_KERNEL)
  return()
endif()

if(NOT DEFINED PYTHON_EXECUTABLE)
  message(STATUS "JupyterKernel: PYTHON_EXECUTABLE not set — skipping kernel registration")
  return()
endif()

set(_KERNEL_JSON "$ENV{HOME}/.local/share/jupyter/kernels/lupnt/kernel.json")

if(NOT EXISTS "${_KERNEL_JSON}")
  message(STATUS "JupyterKernel: registering kernel 'lupnt (pixi)' for ${PYTHON_EXECUTABLE}")
  execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -m ipykernel install --user --name lupnt --display-name
            "lupnt (pixi)"
    RESULT_VARIABLE _KERNEL_RESULT
    OUTPUT_VARIABLE _KERNEL_OUTPUT
    ERROR_VARIABLE _KERNEL_ERROR
  )

  if(NOT _KERNEL_RESULT EQUAL 0)
    message(
      WARNING
        "JupyterKernel: kernel registration failed (exit ${_KERNEL_RESULT}):\n${_KERNEL_ERROR}\n"
        "Run manually: python -m ipykernel install --user --name lupnt --display-name \"lupnt (pixi)\""
    )
    return()
  endif()
endif()

# VSCode (and `jupyter` invoked outside `pixi run`/`pixi shell`) launches the kernelspec's python
# binary directly, without pixi's [activation] env -- so PYTHONPATH (where pylupnt's pure-Python
# package lives) and the data-path env vars are missing, breaking `import pylupnt`. The kernelspec
# "env" field is merged into the kernel subprocess's environment by every Jupyter front-end
# regardless of how it's launched, so patch it here to mirror pixi's [activation] env. Re-applied on
# every configure (cheap) in case pixi.toml's activation env changes.
if(EXISTS "${_KERNEL_JSON}")
  file(READ "${_KERNEL_JSON}" _kernel_json)
  string(
    JSON
    _kernel_json
    SET
    "${_kernel_json}"
    "env"
    "{\"PYTHONPATH\": \"${LUPNT_REPO_ROOT}/python\", \"LUPNT_DATA_PATH\": \"${LUPNT_REPO_ROOT}/data/LuPNT_data\", \"PECSIMPY_BASE_PATH\": \"${LUPNT_REPO_ROOT}/data/LuPNT_data/plasma\", \"LUPNT_OUTPUT_PATH\": \"${LUPNT_REPO_ROOT}/output\"}"
  )
  file(WRITE "${_KERNEL_JSON}" "${_kernel_json}")
  message(
    STATUS "JupyterKernel: kernel 'lupnt' env configured (PYTHONPATH=${LUPNT_REPO_ROOT}/python)"
  )
endif()

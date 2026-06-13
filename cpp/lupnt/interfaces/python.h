#pragma once
#ifdef LUPNT_WITH_PYTHON
#  include <pybind11/eigen.h>
#  include <pybind11/embed.h>
#  include <pybind11/pybind11.h>

#  include <optional>

namespace py = pybind11;

namespace lupnt {

  class PythonInterpreter {
  public:
    /// @brief Get the process-wide embedded-Python interpreter, starting it
    /// on first use.
    ///
    /// Called (via this singleton) before importing any Python module needed
    /// by LuPNT's C++ side -- e.g. the module-level initializers of `py_pnt`,
    /// `py_pio`, `py_go`, `py_plot`, `py_plt` below, which back the plotting
    /// helpers used in `cpp/examples/plotting/`. If LuPNT is itself running
    /// inside an already-initialized Python process (`Py_IsInitialized()`),
    /// no new interpreter is started, avoiding a double-initialization crash
    /// when used from `pylupnt` bindings.
    ///
    /// @return Reference to the single process-wide `PythonInterpreter`
    static PythonInterpreter& GetInstance() {
      static PythonInterpreter instance;
      return instance;
    }

  private:
    std::optional<py::scoped_interpreter> python_;  // Manages Python interpreter

    PythonInterpreter() {
      // Only initialize if we're not already in a Python process
      if (!Py_IsInitialized()) {
        python_.emplace();
      }
    }
    ~PythonInterpreter() = default;

    PythonInterpreter(const PythonInterpreter&) = delete;
    PythonInterpreter& operator=(const PythonInterpreter&) = delete;
  };

  extern py::module_ py_pnt;
  extern py::module_ py_pio;
  extern py::module_ py_go;
  extern py::module_ py_plot;
  extern py::module_ py_plt;
}  // namespace lupnt
#endif  // LUPNT_WITH_PYTHON

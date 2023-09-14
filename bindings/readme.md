# bindings
Python bindings of the lupnt library

## How to add a new python bindings
Suppose you want to add the python bindings to "XXX.cpp" and "XXX.h"
1. Define `void init_xxx(py::module &m)` in `MainBindings.cpp`
2. Add `init_xxx(m)` in the `PYBIND11_MODULE` in `MainBindings.cpp`
3. Create `PyXXX.cpp` in the `bindings/` folder
4. In `PythonXXX.cpp`, define your bindings as follows. Refer to the [official document](https://pybind11.readthedocs.io/en/stable/) on how to write the bindings. 
    ```
    #include <pybind11/pybind11.h>
    #include <lupnt/path/to/XXX.h>

    namespace py = pybind11;
    using namespace LPT;

    void init_state(py::module &m){
        /* define your bindings here */
    }
    ```
5. Compile using CMAKE (nothing special needed, the lupnt project CMAKE will build the bindings as well)
    - If the compiler cannot find the python library, add the following to `.vscode/settings.json`
        ```
        {
            // cmake settings
            "cmake.configureArgs": [
                "-DPYTHON_INCLUDE_DIRS=path1string",
                "-DPYTHON_LIBRARIES=path2string"
                "-DBUILD_EXAMPLES=ON"
            ]
        }
        ```
6. Install the updated python library to local from the root directory of this project. Don't forget to create a new local environment and activate it (`python3 -m venv venv`, and `. venv/bin/activate` or `source .venv/bin/activate`) if you don't want to install it globally. 
    ```
    pip install .
    ```
7. Create a test under `test_python/`
8. Test your bindings with ```python3 -m pytest test_python```

## References and Examples
- [Functions](https://pybind11.readthedocs.io/en/stable/basics.html)
    - [keyward and default arguments](https://pybind11.readthedocs.io/en/stable/basics.html#keyword-args)
    - [with templates](https://pybind11.readthedocs.io/en/stable/advanced/functions.html#binding-functions-with-template-parameters)
- [Class](https://pybind11.readthedocs.io/en/stable/classes.html)
    - [overloading methods](https://pybind11.readthedocs.io/en/stable/classes.html#overloaded-methods)
        - You can do the similar for non-class functions
        - For overwriting constructors, you can do it as follows
            ```
                py::class_<State>(m, "State")
                    .def(py::init<const ad::VectorXreal &>())
                    .def(py::init<const ad::real, const ad::real, const ad::real>())
            ```
    - [operators](https://pybind11.readthedocs.io/en/stable/advanced/classes.html#operator-overloading)
- Defining Global Variables
    ```
    m.attr("PI") = py::float_(3.141592653589793);
    ```
#!/bin/zsh

cd build && cmake --build . --target python-package --config Release
cd lib/python_package && pip install -e .
echo "Pylupnt installed successfully! Version:"
python3 -c "import pylupnt; print(pylupnt.__version__)"
cd ../
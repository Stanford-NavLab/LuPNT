cmake --build build --config Debug --target lupnt pylupnt_pybind python-package -j 14 --
# cp -r build/lib/python_package/pylupnt .venv/lib/python3.12/site-packages
# cp build/third_party/AI-Toolbox/AI2olbox.dylib .venv/lib/python3.12/site-packages/AI2olbox.so
./scripts/stubgen.sh
cp build/lib/python_package/pylupnt/pylupnt_pybind.cpython-*-darwin.so python/pylupnt
echo "PyLuPNT installed"
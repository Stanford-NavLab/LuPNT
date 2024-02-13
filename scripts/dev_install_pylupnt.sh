/opt/homebrew/bin/cmake --build /Users/guillemcv/Development/NavLab/LuPNT/build --config Debug --target python-package -j 14
cp -r build/lib/python_package/pylupnt .venv/lib/python3.12/site-packages
cp build/third_party/AI-Toolbox/AIToolbox.dylib .venv/lib/python3.12/site-packages/AIToolbox.so
echo "PyLuPNT and AI-Toolbox installed in .venv"
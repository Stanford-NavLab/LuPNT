# 1. Dependencies
$ProgressPreference = 'SilentlyContinue'
Invoke-WebRequest -Uri "https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.4.3/hdf5-1.14.4-3-win-vs2022_cl.zip" -OutFile "hdf5.zip"
Expand-Archive -Path "hdf5.zip" -DestinationPath "C:\HDF5" -Force
Remove-Item "hdf5.zip"
Expand-Archive -Path "C:/HDF5/hdf5/HDF5-1.14.4-win64.zip" -DestinationPath "C:/HDF5"
[System.Environment]::SetEnvironmentVariable('HDF5_DIR', "C:\HDF5\HDF5-1.14.4-win64\CMake", 'User')

choco install boost-msvc-14.1

# 2. Download data
$ProgressPreference = 'SilentlyContinue'
Invoke-WebRequest -Uri "https://bit.ly/LuPNT_data" -OutFile "LuPNT_data.zip"
Expand-Archive -Path "LuPNT_data.zip" -DestinationPath "."
Remove-Item "LuPNT_data.zip"

# 3. Set data path for the current session
$env:LUPNT_DATA_PATH = "$PWD\LuPNT_data"

# 4. Set data path permanently (User level)
[System.Environment]::SetEnvironmentVariable('LUPNT_DATA_PATH', "$PWD\LuPNT_data", 'User')

# 5. Check data path
echo "LUPNT_DATA_PATH is set to: $env:LUPNT_DATA_PATH"

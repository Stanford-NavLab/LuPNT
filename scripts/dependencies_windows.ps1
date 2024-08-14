# 1. The compiler should provide OpenMP support
# (No explicit command to install OpenMP on Windows as it's compiler-dependent)

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

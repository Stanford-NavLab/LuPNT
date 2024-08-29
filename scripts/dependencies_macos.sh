#!/bin/bash

# 1. Dependencies
echo "Installing dependencies..."
brew install boost libomp hdf5

# 2. Download data
echo "Downloading data..."
curl -L https://bit.ly/LuPNT_data -o LuPNT_data.zip
unzip LuPNT_data.zip
rm LuPNT_data.zip

# 3. Set data path based on the shell in use
if [ -n "$ZSH_VERSION" ]; then
    SHELL_RC=~/.zshrc
    SHELL_NAME="zsh"
elif [ -n "$BASH_VERSION" ]; then
    SHELL_RC=~/.bashrc
    SHELL_NAME="bash"
else
    echo "Unsupported shell. Please use bash or zsh."
    exit 1
fi

echo "Setting data path in $SHELL_NAME configuration file..."
echo "export LUPNT_DATA_PATH='$(pwd)/LuPNT_data'" >> $SHELL_RC

# 4. Source the shell configuration file to apply the changes
echo "Applying changes by sourcing $SHELL_RC..."
source $SHELL_RC

# 5. Check data path
echo "LUPNT_DATA_PATH is set to: $LUPNT_DATA_PATH"

#!/usr/bin/env bash
set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if conda environment is activated
if [ -z "${CONDA_PREFIX:-}" ]; then
    echo -e "${RED}Error: No conda environment activated.${NC}"
    echo -e "${YELLOW}Please activate a conda environment first: conda activate <env_name>${NC}"
    exit 1
fi

BASEDIR=$(realpath "$(dirname "$0")/..")


# ****************************************************************
# Add python dir to PYTHONPATH by adding it to .zshrc or .bashrc
PYTHON_DIR=$(realpath "$BASEDIR/python")

# Set LUPNT environment variables
LUPNT_PATH_LINE="export LUPNT_PATH=\"$BASEDIR\""
LUPNT_DATA_PATH_LINE="export LUPNT_DATA_PATH=\"\${LUPNT_PATH}/data/LuPNT_data\""
LUPNT_OUTPUT_PATH_LINE="export LUPNT_OUTPUT_PATH=\"\${LUPNT_PATH}/output\""

# Also add CPM_SOURCE_CACHE to .zshrc or .bashrc
CPM_SOURCE_CACHE_LINE='export CPM_SOURCE_CACHE="$HOME/.cache/CPM"'

rc_files=()
[ -f "$HOME/.zshrc" ] && rc_files+=("$HOME/.zshrc")
[ -f "$HOME/.bashrc" ] && rc_files+=("$HOME/.bashrc")

if [ "${#rc_files[@]}" -gt 0 ]; then
    # Array of [description, grep_pattern, line_to_add]
    env_vars=(
        "python dir|export PYTHONPATH=.*$PYTHON_DIR|export PYTHONPATH=\$PYTHONPATH:$PYTHON_DIR"
        "CPM_SOURCE_CACHE|^export CPM_SOURCE_CACHE=|$CPM_SOURCE_CACHE_LINE"
        "LUPNT_PATH|^export LUPNT_PATH=|$LUPNT_PATH_LINE"
        "LUPNT_DATA_PATH|^export LUPNT_DATA_PATH=|$LUPNT_DATA_PATH_LINE"
        "LUPNT_OUTPUT_PATH|^export LUPNT_OUTPUT_PATH=|$LUPNT_OUTPUT_PATH_LINE"
    )
    for rc_file in "${rc_files[@]}"; do
        for entry in "${env_vars[@]}"; do
            IFS="|" read -r desc pattern line <<< "$entry"
            if ! grep -q "$pattern" "$rc_file"; then
                echo -e "${YELLOW}Adding $desc to $(basename "$rc_file")...${NC}"
                echo "$line" >> "$rc_file"
            else
                echo -e "${GREEN}Found $desc in ~/${rc_file##*/}${NC}"
            fi
        done
    done
else
    echo -e "${RED}No .zshrc or .bashrc found. Please add PYTHONPATH, CPM_SOURCE_CACHE, and LUPNT variables manually.${NC}"
fi

# ****************************************************************
# Create .vscode/settings.json if it doesn't exist or is invalid
SETTINGS_JSON="$BASEDIR/.vscode/settings.json"
mkdir -p "$(dirname "$SETTINGS_JSON")"

if [ ! -f "$SETTINGS_JSON" ]; then
    echo -e "${YELLOW}Creating $(realpath "$SETTINGS_JSON")${NC}"
    echo '{}' > "$SETTINGS_JSON"
    echo -e "${GREEN}Created .vscode/settings.json${NC}"
elif ! jq . "$SETTINGS_JSON" >/dev/null 2>&1; then
    echo -e "${YELLOW}Invalid JSON in $SETTINGS_JSON, reinitializing...${NC}"
    echo '{}' > "$SETTINGS_JSON"
    echo -e "${GREEN}Reinitialized .vscode/settings.json${NC}"
else
    echo -e "${GREEN}Found valid .vscode/settings.json${NC}"
fi

# Detect platform and set compilers accordingly
if [[ "$OSTYPE" == "darwin"* ]]; then
  cc_bin="$CONDA_PREFIX/bin/clang"
  cxx_bin="$CONDA_PREFIX/bin/clang++"
else
    # --- CHANGE THIS BLOCK FOR LINUX ---
    # Use the environment variables set by Conda's compiler metapackage
    cc_bin="${CC:-}"
    cxx_bin="${CXX:-}"

    # Fallback to the generic names if CC/CXX are not set (less reliable)
    if [ -z "${CC:-}" ]; then
        echo "Warning: CC variable not found. Using generic names. Ensure 'conda activate' was run."
        cc_bin="${CONDA_PREFIX}/bin/gcc"
        cxx_bin="${CONDA_PREFIX}/bin/g++"
    fi
    # -----------------------------------
fi

# Add or update cmake settings
if jq --arg cmakePath "$CONDA_PREFIX/bin/cmake" \
      --arg conda_prefix "$CONDA_PREFIX" \
      --arg cc "$cc_bin" \
      --arg cxx "$cxx_bin" \
      --arg path "$CONDA_PREFIX/bin:\${env:PATH}" \
      --arg install_rpath "$CONDA_PREFIX/lib" \
      --arg build_rpath "$CONDA_PREFIX/lib" \
      '.["cmake.cmakePath"]=$cmakePath |
       .["cmake.generator"]="Ninja" |
       .["cmake.environment"]["CONDA_PREFIX"]=$conda_prefix |
       .["cmake.environment"]["PATH"]=$path |
       .["cmake.configureSettings"]["CMAKE_C_COMPILER"]=$cc |
       .["cmake.configureSettings"]["CMAKE_CXX_COMPILER"]=$cxx |
       .["cmake.configureSettings"]["CMAKE_PREFIX_PATH"]=$conda_prefix |
       .["cmake.configureSettings"]["CMAKE_INSTALL_RPATH"]=$install_rpath |
       .["cmake.configureSettings"]["CMAKE_BUILD_RPATH"]=$build_rpath' \
      "$SETTINGS_JSON" > tmp.json; then
  mv tmp.json "$SETTINGS_JSON"
else
  echo -e "${RED}jq failed to update $SETTINGS_JSON. Please check the file for JSON validity.${NC}"
  exit 1
fi

# ****************************************************************
# Create or update .clangd file
CLANGD_FILE="$BASEDIR/.clangd"
echo -e "Creating/updating $(realpath "$CLANGD_FILE")"

cat > "$CLANGD_FILE" << EOF
CompileFlags:
  Remove: ["-fopenmp=libomp"]
  Add: ["-isystem", "$CONDA_PREFIX/include"]
EOF

echo -e "${GREEN}Created .clangd configuration${NC}"

# ****************************************************************
echo -e "${GREEN}Development environment installed${NC}"
echo -e "${YELLOW}Please restart your IDE to apply changes.${NC}"

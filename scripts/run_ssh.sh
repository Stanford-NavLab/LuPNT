#!/bin/zsh
setopt MONITOR
# Get current date and time
CURRENT_DATETIME=$(date +"%y%m%d-%H%M%S")

# Check if required arguments were provided
if [ $# -lt 2 ]; then
    echo "Please provide a Python file and output directory"
    echo "Usage: ./run_ssh.sh <python_file> <output_dir> <args...>"
    exit 1
fi

# Get the Python file and output directory
PYTHON_FILE=$1
OUTPUT_DIR=$2
shift 2
ARGS="$@"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Run the script and redirect output to a file with the current date and time
nohup python $PYTHON_FILE $ARGS > "$OUTPUT_DIR/$CURRENT_DATETIME.txt" & disown

# Tail the output file
tail -f "$OUTPUT_DIR/$CURRENT_DATETIME.txt"

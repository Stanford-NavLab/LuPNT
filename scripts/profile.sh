#!/bin/bash

# Set the output directory
OUTPUT_DIR="$(dirname "${BASH_SOURCE[0]}")/../profiling"

# Set the template
TEMPLATE="Time Profiler"

# Check if the target executable is provided
if [ $# -eq 0 ]; then
  echo "Error: Target executable not provided"
  exit 1
fi

# Make directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run xctrace
xctrace record --output "$OUTPUT_DIR" --template "$TEMPLATE" --launch "$1" --target-stdout /dev/stdout
